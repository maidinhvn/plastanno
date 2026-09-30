#!/usr/bin/env python3
"""A CDS ends at its first in-frame stop, and a real repair keeps its warning.

Run DIRECTLY and check $?. Exit 0 = all pass.

The terminal-stop step used to look only 3 codons ahead, within a length guard. On development
genomes that left most truncated 3' ends without a stop, and never trimmed a read-through past
a stop (ndhB). A first redesign repaired both, but a repaired ORF passes NCBI's validator, so
the warning disappeared from gene models still wrong at their start or splice sites. Now:

  1 trim      a read-through of <= 5 codons past the first stop is cut back to it   (marked)
  2           a CDS already ending at a stop is left alone
  3 old snap  a stop within 3 codons, within the old length guard, is taken silently,
              whatever lies upstream -- exactly what 3.0.1 did
  4 complete  with no internal stop, the 3' end is extended to the first stop within
              100 codons                                                              (marked)
  5           anything else is left alone

and finalize.submission_check keeps every marked CDS NEEDS_REVIEW, even when its ORF is valid.

Every case is built on a synthetic genome whose filler (GCC repeats) has no stop codon in any
frame on either strand, so each case knows its answer in advance.
"""
import copy
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.core.feature import Feature                                   # noqa: E402
from plastanno.core import coords as C                                       # noqa: E402
from plastanno.core.finalize import submission_check, SUBMIT, COMPLETED_3P   # noqa: E402
from plastanno.annotate import refine_splice as RS                           # noqa: E402
from plastanno.output.writers import write_report                           # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-74s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
L = 9000
G = list("GCC" * (L // 3))


def put(pos, s):
    G[pos:pos + len(s)] = list(s)


def cds(name, s, e, strand=1, exons=None, **kw):
    f = Feature(gene_name=name, gene_type="CDS", start=s, end=e, strand=strand, engine="AB",
                flag=kw.pop("flag", "HIGH"), **kw)
    f.exons = exons or [(s, e)]
    return f


BODY = "ATG" + "GCC" * 10                         # 11 codons, no stop
cases = {}
put(100, BODY + "TAA");                   cases["adj"] = cds("adj", 100, 133)         # old snap, 1 codon
put(250, BODY + "GCC" * 2 + "TGA");       cases["adj3"] = cds("adj3", 250, 283)       # old snap, 3 codons
put(400, BODY + "GCC" * 19 + "TAG");      cases["ext20"] = cds("ext20", 400, 433)     # stop at 490..493
put(900, BODY);                           cases["none"] = cds("none", 900, 933)       # no stop within 100 codons
put(1200, BODY + "TAA");                  cases["guard"] = cds("guard", 1200, 1233)   # 1 codon, old guard blocks
put(1800, BODY + "TAA" + "TTT");          cases["trim1"] = cds("trim1", 1800, 1839)
put(2100, BODY + "TAA" + "GCC" * 3 + "TAG"); cases["trim4"] = cds("trim4", 2100, 2148)
# an internal stop far upstream: the old snap still takes an adjacent stop (rps11) ...
put(2400, "ATG" + "GCC" * 3 + "TAA" + "GCC" * 10 + "TGA"); cases["far_adj"] = cds("far_adj", 2400, 2445)   # TGA at 2445..2448
# ... but nothing completes it from further away
put(2600, "ATG" + "GCC" * 3 + "TAA" + "GCC" * 10 + "GCC" * 6 + "TAG"); cases["far_far"] = cds("far_far", 2600, 2651)
put(2700, BODY + "TAA");                  cases["done"] = cds("done", 2700, 2736)
put(3000, BODY + "GC");                   cases["frame"] = cds("frame", 3000, 3035)
put(3300, "ATG" + "GCC" * 3); put(3400, "GCC" * 4 + "GCC" * 6 + "TAA")
cases["twoexon"] = cds("twoexon", 3300, 3412, exons=[(3300, 3312), (3400, 3412)])      # stop 6 codons on
put(3600, "ATG" + "GCC" * 6 + "TAA"); put(3700, "GCC")
cases["cutexon"] = cds("cutexon", 3600, 3703, exons=[(3600, 3624), (3700, 3703)])
full = BODY + "GCC" * 7 + "TAA"
put(4200 - (len(full) - len(BODY)), rc(full))
cases["mext"] = cds("mext", 4200, 4200 + len(BODY), strand=-1)                        # 8 codons
put(4800, rc(BODY + "TAG" + "GCC" * 2 + "TAA"))                                      # 3 codons past TAG
cases["mtrim"] = cds("mtrim", 4800, 4800 + len(BODY) + 12, strand=-1)
put(L - 33, BODY); put(0, "GCC" * 5 + "TAA")
cases["wrap"] = cds("wrap", L - 33, L)                                               # 6 codons, across the origin
put(6000, BODY + "GCC" * 3 + "TAA"); cases["pseudo"] = cds("pseudo", 6000, 6033, is_pseudogene=True)
SEQ = "".join(G)
CAT = {"guard": {"expected_len": 28}}             # 33 + 3 = 36 > 1.2 x 28


def run(f):
    g = copy.deepcopy(f)
    return RS._snap_terminal_stop(g, SEQ, CAT), g


marked = lambda f: any(str(n).startswith(COMPLETED_3P) for n in f.notes)
last = lambda f: C.extract(SEQ, f, L)[-3:]

print("--- 3: the old snap, silent ---")
ch, f = run(cases["adj"])
check("an adjacent stop is taken, without a mark (as 3.0.1 did)", (ch, f.end, last(f), marked(f)), (True, 136, "TAA", False))
ch, f = run(cases["adj3"])
check("a stop 3 codons on is taken, without a mark", (ch, f.end, marked(f)), (True, 292, False))
ch, f = run(cases["far_adj"])
check("with an internal stop far upstream, the adjacent stop is still taken (rps11)",
      (ch, f.end, last(f), marked(f)), (True, 2448, "TGA", False))

print("--- 4: completion, marked ---")
ch, f = run(cases["ext20"])
check("a stop 20 codons on is reached, and the CDS is marked", (ch, f.end, last(f), marked(f)), (True, 493, "TAG", True))
ch, f = run(cases["guard"])
check("a 1-codon snap the old guard blocks is completed instead, and marked",
      (ch, f.end, marked(f)), (True, 1236, True))
ch, f = run(cases["twoexon"])
check("two exons: the LAST exon is completed, the first untouched, marked",
      (ch, f.exons, marked(f)), (True, [(3300, 3312), (3400, 3433)], True))
ch, f = run(cases["mext"])
check("minus strand: completed toward lower coordinates, marked",
      (ch, f.start, last(f), marked(f)), (True, 4200 - 24, "TAA", True))
ch, f = run(cases["wrap"])
check("a stop across the origin is reached, the CDS now wraps, marked",
      (ch, f.end, sorted(f.exons), last(f), marked(f)), (True, 18, [(0, 18), (L - 33, L)], "TAA", True))
ch, f = run(cases["none"])
check("no stop within 100 codons: unchanged", (ch, f.end, f.notes), (False, 933, []))
ch, f = run(cases["far_far"])
check("an internal stop far upstream and no nearby stop: not completed", (ch, f.end, f.notes), (False, 2651, []))

print("--- 1: trim, marked ---")
ch, f = run(cases["trim1"])
check("one codon past a stop is trimmed back, marked", (ch, f.end, last(f), marked(f)), (True, 1836, "TAA", True))
ch, f = run(cases["trim4"])
check("four codons past a stop are trimmed back, marked", (ch, f.end, last(f), marked(f)), (True, 2136, "TAA", True))
ch, f = run(cases["mtrim"])
check("minus strand: trimmed from the lower end, marked", (ch, f.start, last(f), marked(f)), (True, 4809, "TAG", True))
ch, f = run(cases["cutexon"])
check("a read-through longer than the last exon is left alone", (ch, f.exons), (False, [(3600, 3624), (3700, 3703)]))

print("--- left alone ---")
ch, f = run(cases["done"])
check("a CDS already ending at its stop is untouched", (ch, f.notes), (False, []))
ch, f = run(cases["frame"])
check("a broken frame is untouched", ch, False)

print("--- refine_all ---")
ann = [copy.deepcopy(cases[k]) for k in ("ext20", "trim1", "pseudo", "adj")]
ann[-1].exon_strands = [1]                        # the trans-spliced marker
RS.refine_all(ann, SEQ, CAT)
by = {a.gene_name: a for a in ann}
check("refine_all completes and trims ordinary CDS", (by["ext20"].end, by["trim1"].end), (493, 1836))
check("refine_all leaves a pseudogene and a trans-spliced CDS alone",
      (by["pseudo"].end, by["adj"].end, by["pseudo"].notes, by["adj"].notes), (6033, 133, [], []))

print("--- the warning stays ---")
done = {k: run(cases[k])[1] for k in ("adj", "ext20", "trim4", "guard")}
flagged = {f.gene_name: p for f, p in submission_check(list(done.values()), SEQ)}
check("a silently snapped CDS with a valid ORF is not flagged", done["adj"].flag, "HIGH")
check("a completed CDS with a valid ORF stays NEEDS_REVIEW", (done["ext20"].flag, "ext20" in flagged), ("NEEDS_REVIEW", True))
check("  ... and its note says why", any("3' end was completed by the pipeline" in str(n) and "check its start" in str(n)
                                          for n in done["ext20"].notes if str(n).startswith(SUBMIT)), True)
check("a trimmed CDS with a valid ORF stays NEEDS_REVIEW", done["trim4"].flag, "NEEDS_REVIEW")
check("a guard-blocked snap, completed, stays NEEDS_REVIEW", done["guard"].flag, "NEEDS_REVIEW")
snap = [(f.flag, list(f.notes)) for f in done.values()]
submission_check(list(done.values()), SEQ)
check("a second pass changes no flag and no note", [(f.flag, list(f.notes)) for f in done.values()], snap)
with tempfile.TemporaryDirectory() as d:
    p = os.path.join(d, "x.report")
    write_report(list(done.values()), "TEST", L, {}, [], 1.0, p)
    rep = open(p).read()
    check("the report lists completed CDS apart from NCBI failures",
          ("No CDS would fail NCBI validation" in rep,
           "3 CDS pass only because the pipeline completed their 3' end" in rep), (True, True))

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
