#!/usr/bin/env python3
"""A first codon that is not a start codon moves to the best-supported start codon nearby,
and the CDS keeps its warning.

Run DIRECTLY and check $?. Exit 0 = all pass.

NCBI rejects a CDS whose first codon is neither a table-11 initiation codon nor written as an
RNA-edited start (SPECIAL_START_CODONS). On development genomes most of them sat one to three
codons from the reference start: rpl22 one codon upstream of its ATG, rpl2 on the GCG after
its edited ACG. The start-codon pass in refine_splice.refine_all:

  - touches only a first codon that is not a start codon (a valid start is never moved);
  - looks up to 10 codons upstream (never across an in-frame stop, round the origin) and
    downstream (inside the first exon, leaving at least 3 bp of it);
  - takes ATG or one of the gene's special start codons, never another table-11 start;
  - picks the candidate whose N-terminus the reference proteins support best; ties go to ATG,
    then to the nearest, then upstream;
  - does not move if the current start is better supported than the winner;

and finalize.submission_check keeps every moved CDS NEEDS_REVIEW.

Every case is built on a synthetic genome whose filler (GCC repeats) has no stop codon in any
frame on either strand, and a protein database written for the test, so each case knows its
answer in advance.
"""
import copy
import io
import os
import sys
import tempfile
from contextlib import redirect_stdout

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.core.feature import Feature                                   # noqa: E402
from plastanno.core import coords as C                                       # noqa: E402
from plastanno.core.finalize import submission_check, finalize_qc, SUBMIT, MOVED_5P   # noqa: E402
from plastanno.annotate import refine_splice as RS                           # noqa: E402
from plastanno.output.writers import write_report                           # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-78s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
L = 12000
G = list("GCC" * (L // 3))


def put(pos, s):
    G[pos:pos + len(s)] = list(s)


def cds(name, s, e, strand=1, exons=None, **kw):
    f = Feature(gene_name=name, gene_type="CDS", start=s, end=e, strand=strand, engine="AB",
                flag=kw.pop("flag", "HIGH"), **kw)
    f.exons = exons or [(s, e)]
    return f


# Reference proteins, one file per gene. N-termini the cases are built to match.
DB = tempfile.mkdtemp(prefix="start_rescue_db_")
REFS = {
    "tA": "MKLVAAAAAAAAAAAAAAAA",        # downstream ATG cases
    "tB": "MAANKLVAAAAAAAAAAAAA",        # upstream ATG cases
    "tC": "TKIRMSVRKICEAAAAAAAA",        # the current (invalid) start is the supported one
    "tD": "MDEWAAAAAAAAAAAAAAAA",        # homology beats distance
    "tE": "MAAAAAAAAAAAAAAAAAAA",
    "tF": "MWWWWWWWWWWWWWWWWWWW",        # supports no candidate: ties decide
    "tG": "MSSWAAAAAAAAAAAAAAAA",        # a GTG start is better supported, but is no candidate
    "rpl2": "MAIHLYKTSTPSAAAAAAAA",
}
for gene, prot in REFS.items():
    with open(os.path.join(DB, gene + ".fasta"), "w") as fh:
        for i in range(3):
            fh.write(">%s_%d\n%s\n" % (gene, i, prot))

TAIL = "GCC" * 12 + "TAA"
cases = {}
# plus strand, the ATG one codon downstream of an AAC start
put(100, "AAC" + "ATG" + "AAACTGGTT" + TAIL)
cases["down1"] = cds("tA", 100, 100 + 3 + 3 + 9 + len(TAIL))
# minus strand, the same
s = "AAC" + "ATG" + "AAACTGGTT" + TAIL
put(400, rc(s)); cases["mdown1"] = cds("tA", 400, 400 + len(s), strand=-1)
# plus strand, the ATG three codons upstream: ATG GCC GCC [AAC AAA CTG GTT ...
put(700, "ATG" + "GCCGCC" + "AAC" + "AAACTGGTT" + TAIL)
cases["up3"] = cds("tB", 709, 709 + 3 + 9 + len(TAIL))
# minus strand, the same
s = "ATG" + "GCCGCC" + "AAC" + "AAACTGGTT" + TAIL
put(1000, rc(s)); cases["mup3"] = cds("tB", 1000, 1000 + len(s) - 9, strand=-1)
# the ATG 11 codons upstream: out of the window
put(1300, "ATG" + "GCC" * 10 + "AAC" + "AAACTGGTT" + TAIL)
cases["far"] = cds("tB", 1300 + 33, 1300 + 33 + 3 + 9 + len(TAIL))
# an in-frame stop between the ATG and the start
put(1600, "ATG" + "TAA" + "GCC" + "AAC" + "AAACTGGTT" + TAIL)
cases["stop"] = cds("tB", 1609, 1609 + 3 + 9 + len(TAIL))
# valid starts are never moved, even with a better ATG next door
put(1900, "GTG" + "ATG" + "AAACTGGTT" + TAIL); cases["gtg"] = cds("tA", 1900, 1900 + 15 + len(TAIL))
put(2200, "ATA" + "ATG" + "AAACTGGTT" + TAIL); cases["ata"] = cds("tA", 2200, 2200 + 15 + len(TAIL))
# rpl2: ACG [GCG ATA CAT TTA TAC AAA ACT TCT ACT CCT TCT ...; ACG is rpl2's special start
RPL2 = "ACG" + "GCGATACATTTATACAAAACTTCTACTCCTTCT" + TAIL
put(2500, RPL2); cases["rpl2"] = cds("rpl2", 2503, 2500 + len(RPL2))
# ... ACG is not a candidate for a gene without it among its special starts
put(2800, RPL2); cases["notspecial"] = cds("tA", 2803, 2800 + len(RPL2))
# the current start is better supported than the ATG four codons on (a fern ACG):
# ACG AAA ATT CGT [ATG TCT GTT CGT AAA ATT TGT GAA reads TKIRMSVRKICE, the ATG MSVRKICE...
put(3100, "ACGAAAATTCGT" + "ATGTCTGTTCGTAAAATTTGTGAA" + TAIL)
cases["guard"] = cds("tC", 3100, 3100 + 36 + len(TAIL))
# homology beats distance: ATG +1 (MAAA..., unsupported) or ATG -4 (MDEW..., supported)
put(3400, "ATGGATGAATGG" + "AAC" + "ATG" + "GCC" * 3 + TAIL)
cases["hom"] = cds("tD", 3412, 3412 + 3 + 3 + 9 + len(TAIL))
# ties: ATG over a special codon; nearest; then upstream
put(3700, "ACG" + "GCG" + "ATG" + TAIL)                       # rpl2: ACG -1 and ATG +1, both unsupported
cases["tie_atg"] = cds("rpl2", 3703, 3706 + 3 + len(TAIL))
put(4000, "ATG" + "GCC" + "AAC" + "ATG" + TAIL)                # ATG -2 and ATG +1
cases["tie_near"] = cds("tF", 4006, 4009 + 3 + len(TAIL))
put(4300, "ATG" + "AAC" + "ATG" + TAIL)                        # ATG -1 and ATG +1
cases["tie_up"] = cds("tF", 4303, 4306 + 3 + len(TAIL))
# a better-supported GTG next door is not a candidate; the ATG further on is
put(4600, "AAC" + "GTG" + "AGTTCATGG" + "ATG" + TAIL)
cases["t11"] = cds("tG", 4600, 4600 + 15 + 3 + len(TAIL))
# two exons: the ATG in the first exon is taken, the second exon is untouched
put(5000, "AAC" + "ATG" + "AAA"); put(5200, "CTGGTT" + TAIL)
cases["twoexon"] = cds("tA", 5000, 5200 + 6 + len(TAIL), exons=[(5000, 5009), (5200, 5206 + len(TAIL))])
# ... down to exactly 3 bp of the first exon: AAC GCC [ATG | AAA ...
put(5500, "AAC" + "GCC" + "ATG"); put(5700, "AAACTGGTT" + TAIL)
cases["exon3bp"] = cds("tA", 5500, 5700 + 9 + len(TAIL), exons=[(5500, 5509), (5700, 5709 + len(TAIL))])
# ... but not a start codon the intron splits: AAC GCC AT | G AAA ... (2 bp would remain)
put(5800, "AAC" + "GCC" + "AT"); put(5900, "G" + "AAACTGGTT" + TAIL)
cases["splitcodon"] = cds("tA", 5800, 5900 + 10 + len(TAIL), exons=[(5800, 5808), (5900, 5910 + len(TAIL))])
# across the origin, plus strand: ATG at L-3, then [GCC | AAC ... from 0
put(L - 3, "ATG"); put(0, "GCC")
put(3, "AAC" + "AAACTGGTT" + TAIL)
cases["wrap"] = cds("tB", 3, 3 + 3 + 9 + len(TAIL))
# left alone: pseudogene, trans-spliced, frame, internal stop, no reference proteins
put(6000, "AAC" + "ATG" + "AAACTGGTT" + TAIL)
cases["pseudo"] = cds("tA", 6000, 6000 + 15 + len(TAIL), is_pseudogene=True)
cases["trans"] = cds("tA", 6000, 6000 + 15 + len(TAIL), exon_strands=[1])
put(6300, "AAC" + "ATG" + "AAACTGGTT" + TAIL + "G")
cases["frame"] = cds("tA", 6300, 6300 + 15 + len(TAIL) + 1)
put(6600, "AAC" + "ATG" + "AAATAACTGGTT" + TAIL)
cases["internal"] = cds("tA", 6600, 6600 + 18 + len(TAIL))
put(6900, "AAC" + "ATG" + "AAACTGGTT" + TAIL)
cases["nodb"] = cds("tZ", 6900, 6900 + 15 + len(TAIL))
SEQ = "".join(G)
# across the origin, minus strand, on a second genome (one origin each): the 5' end at L-3,
# its upstream codons at L-3..L and 0..3
G2 = list("GCC" * (L // 3))
s = "AAC" + "AAACTGGTT" + TAIL                                  # coding orientation
G2[L - 3 - len(s):L - 3] = list(rc(s))
G2[L - 3:L] = list(rc("GCC")); G2[0:3] = list(rc("ATG"))
SEQ2 = "".join(G2)
cases["mwrap"] = cds("tB", L - 3 - len(s), L - 3, strand=-1)


def run(f, db=DB, seq=None):
    g = copy.deepcopy(f)
    return RS._rescue_invalid_start(g, seq or SEQ, db, {}), g


moved = lambda f: any(str(n).startswith(MOVED_5P) for n in f.notes)
first = lambda f: C.extract(SEQ, f, L)[:3]
five = lambda f: f.start if f.strand == 1 else f.end

print("--- an invalid start moves, on either strand ---")
ch, f = run(cases["down1"])
check("plus: AAC -> the ATG 1 codon downstream, marked", (ch, f.start, f.end, first(f), moved(f)),
      (True, 103, cases["down1"].end, "ATG", True))
ch, f = run(cases["mdown1"])
check("minus: AAC -> the ATG 1 codon downstream, marked", (ch, f.start, f.end, first(f), moved(f)),
      (True, cases["mdown1"].start, cases["mdown1"].end - 3, "ATG", True))
ch, f = run(cases["up3"])
check("plus: AAC -> the ATG 3 codons upstream", (ch, f.start, first(f)), (True, 700, "ATG"))
ch, f = run(cases["mup3"])
check("minus: AAC -> the ATG 3 codons upstream", (ch, f.end, first(f)), (True, cases["mup3"].end + 9, "ATG"))
check("... the note names both codons and the distance",
      any("start moved from AAC to ATG, 3 codon(s) upstream" in str(n) for n in f.notes), True)
check("... and the rest of the CDS is unchanged",
      C.extract(SEQ, f, L)[9:], C.extract(SEQ, cases["mup3"], L))

print("--- where it may look ---")
check("an ATG 11 codons upstream is out of the window", run(cases["far"])[0], False)
check("an ATG beyond an in-frame stop is out of reach", run(cases["stop"])[0], False)
ch, f = run(cases["twoexon"])
check("two exons: the first is trimmed to its ATG, the second untouched",
      (ch, f.exons, first(f)), (True, [(5003, 5009), (5200, 5206 + len(TAIL))], "ATG"))
ch, f = run(cases["exon3bp"])
check("... down to exactly 3 bp of the first exon", (ch, f.exons), (True, [(5506, 5509), (5700, 5709 + len(TAIL))]))
check("... but not a start codon the intron splits (2 bp would remain)", run(cases["splitcodon"])[0], False)
ch, f = run(cases["wrap"])
check("across the origin (plus): the start moves to L-3 and the CDS wraps",
      (ch, f.start, f.end, sorted(f.exons), first(f)),
      (True, L - 3, cases["wrap"].end, [(0, cases["wrap"].end), (L - 3, L)], "ATG"))
check("... and reads ATG GCC then the old CDS", C.extract(SEQ, f, L), "ATGGCC" + C.extract(SEQ, cases["wrap"], L))
ch, f = run(cases["mwrap"], seq=SEQ2)
check("across the origin (minus): the 5' end moves to 3 and the CDS wraps",
      (ch, f.start, f.end, sorted(f.exons), C.extract(SEQ2, f, L)[:3]),
      (True, cases["mwrap"].start, 3, [(0, 3), (cases["mwrap"].start, L)], "ATG"))
check("... and reads ATG GCC then the old CDS", C.extract(SEQ2, f, L), "ATGGCC" + C.extract(SEQ2, cases["mwrap"], L))

print("--- which codons ---")
check("a GTG start is valid: never moved", run(cases["gtg"])[0], False)
check("an ATA start is valid: never moved", run(cases["ata"])[0], False)
ch, f = run(cases["rpl2"])
check("rpl2 on GCG -> its special ACG one codon upstream", (ch, f.start, first(f)), (True, 2500, "ACG"))
check("ACG is no candidate for a gene without it among its special starts", run(cases["notspecial"])[0], False)
ch, f = run(cases["t11"])
check("a better-supported GTG is no candidate; the ATG further on is taken", (ch, first(f), f.start), (True, "ATG", 4615))

print("--- which one ---")
ch, f = run(cases["hom"])
check("homology beats distance: the supported ATG 4 codons upstream", (ch, f.start, first(f)), (True, 3400, "ATG"))
ch, f = run(cases["tie_atg"])
check("equal support: ATG (+1) over the special ACG (-1)", (ch, f.start), (True, 3706))
ch, f = run(cases["tie_near"])
check("equal support: the nearer ATG (+1) over the farther (-2)", (ch, f.start), (True, 4009))
ch, f = run(cases["tie_up"])
check("equal support and distance: upstream (-1) over downstream (+1)", (ch, f.start), (True, 4300))
check("no move when the current start is better supported (a fern ACG)", run(cases["guard"])[0], False)

print("--- left alone ---")
check("a pseudogene is not the pass's business (refine_all skips it)", RS.refine_all(
    [copy.deepcopy(cases["pseudo"])], SEQ, {}, protein_db=DB), 0)
p = copy.deepcopy(cases["pseudo"]); t = copy.deepcopy(cases["trans"])
with redirect_stdout(io.StringIO()):
    RS.refine_all([p, t], SEQ, {}, protein_db=DB)
check("... nor a pseudogene or a trans-spliced CDS through refine_all", (moved(p), moved(t)), (False, False))
check("a CDS that is not a whole number of codons", run(cases["frame"])[0], False)
check("a CDS with an internal stop", run(cases["internal"])[0], False)
check("a gene without reference proteins", run(cases["nodb"])[0], False)
check("no protein database at all", run(cases["down1"], db=None)[0], False)

print("--- through refine_all, and the warning is kept ---")
anns = [copy.deepcopy(cases[k]) for k in ("down1", "rpl2", "gtg")]
out = io.StringIO()
with redirect_stdout(out):
    RS.refine_all(anns, SEQ, {}, protein_db=DB)
check("refine_all moves the two invalid starts and says so",
      ([moved(a) for a in anns], "Start codon moved: 2 CDS" in out.getvalue()), ([True, True, False], True))
quiet = [copy.deepcopy(cases["down1"])]
with redirect_stdout(io.StringIO()):
    RS.refine_all(quiet, SEQ, {})
check("refine_all without a protein database moves nothing", moved(quiet[0]), False)
with redirect_stdout(io.StringIO()):
    finalize_qc(anns, SEQ)
check("finalize_qc reads the moved rpl2 ACG as an RNA-edited start", anns[1].rna_edited_start, True)
bad = submission_check(anns, SEQ)
check("submission_check keeps both moved CDS NEEDS_REVIEW, the GTG one untouched",
      ([a.flag for a in anns], len(bad)), (["NEEDS_REVIEW", "NEEDS_REVIEW", "HIGH"], 2))
note = next(str(n) for n in anns[0].notes if str(n).startswith(SUBMIT))
check("... with a note that the pipeline moved the start",
      note, SUBMIT + "its start codon was moved by the pipeline (start moved from AAC to ATG, 1 codon(s) "
      "downstream; AAC is not a start codon): check its start and exon structure before "
      "submitting; flag lowered from HIGH")
check("... and the reason it returns", bad[0][1], ["start codon moved by the pipeline"])
again = submission_check(anns, SEQ)
check("idempotent: a second pass gives the same flags and notes",
      ([a.flag for a in anns], [n for n in anns[0].notes if str(n).startswith(SUBMIT)], len(again)),
      (["NEEDS_REVIEW", "NEEDS_REVIEW", "HIGH"], [note], 2))
with tempfile.TemporaryDirectory() as d:
    rep = os.path.join(d, "r.report")
    write_report(anns, "TEST", L, {}, [], 1.0, rep)
    txt = open(rep).read()
check("the report lists them apart, as passing only because the start was moved",
      "2 CDS pass only because the pipeline moved their start codon" in txt, True)

print("\n%d checks, %d failed" % (RUN[0], len(FAIL)))
sys.exit(1 if FAIL else 0)
