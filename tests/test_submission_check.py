#!/usr/bin/env python3
"""Every CDS that NCBI's validator would reject is flagged NEEDS_REVIEW, and nothing moves.

Run DIRECTLY and check $?. Exit 0 = all pass.

Most CDS that table2asn rejects (no valid start, no stop, internal stops) used to reach the
user flagged MEDIUM or HIGH. finalize.submission_check applies table2asn's three CDS checks for
genetic code 11 and flags what fails, after the rescue gate, changing flags and notes only.

Every case below is built on a synthetic genome whose filler (GCC repeats) has no stop codon
in any frame on either strand, so each case knows its answer in advance.
"""
import copy
import inspect
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.core.feature import Feature                                   # noqa: E402
from plastanno.core import finalize as FZ                                    # noqa: E402
from plastanno.core.finalize import submission_check, finalize_qc, SUBMIT    # noqa: E402
from plastanno.output.writers import write_report                           # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-70s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
GENOME = list("GCC" * 2000)                      # 6000 bp, no stop in any frame, either strand


def place(pos, coding, strand=1):
    """Write a coding sequence into the genome; return the 0-based half-open span."""
    s = coding if strand == 1 else rc(coding)
    GENOME[pos:pos + len(s)] = list(s)
    return pos, pos + len(s)


def cds(name, pos, coding, strand=1, **kw):
    s, e = place(pos, coding, strand)
    f = Feature(gene_name=name, gene_type="CDS", start=s, end=e, strand=strand, engine="AB",
                flag=kw.pop("flag", "HIGH"), **kw)
    f.exons = [(s, e)]
    return f


BODY = "GCC" * 10
F = {}
F["valid"] = cds("validA", 100, "ATG" + BODY + "TAA")
F["gtg"] = cds("gtgB", 300, "GTG" + BODY + "TAG")
F["bad_start"] = cds("badC", 500, "CCC" + BODY + "TAA")
F["acg_edited"] = cds("psbL", 700, "ACG" + BODY + "TGA", rna_edited_start=True)
F["acg_plain"] = cds("acgD", 900, "ACG" + BODY + "TGA")
F["no_stop"] = cds("nostopE", 1100, "ATG" + BODY + "GCC")
F["frame"] = cds("frameF", 1300, "ATG" + BODY + "TAAG")           # 3' end not a whole codon
F["internal"] = cds("intG", 1500, "ATG" + "GCC" * 3 + "TAA" + "GCC" * 3 + "TAA")
F["minus_ok"] = cds("minusH", 1700, "ATG" + BODY + "TAA", strand=-1)
F["minus_bad"] = cds("minusI", 1900, "CCC" + BODY + "TAA", strand=-1)
F["partial"] = cds("partJ", 2100, "CCC" + BODY + "GCC", orf_incomplete=True)
F["partial_int"] = cds("partK", 2300, "CCC" + "GCC" * 3 + "TAG" + "GCC" * 3 + "GC",
                       orf_incomplete=True)
F["pseudo"] = cds("pseuL", 2500, "CCC" + BODY + "GCC", is_pseudogene=True)
F["already"] = cds("revM", 2700, "CCC" + BODY + "TAA", flag="NEEDS_REVIEW")
# two exons: the intron holds stop codons, the spliced sequence none
e1 = "ATG" + "GCC" * 3 + "GC"
e2 = "C" + "GCC" * 3 + "TAA"
place(2900, e1); place(2914, "GTAAGTAGTGATAA" + "TAG" * 10 + "AG"); place(2960, e2)
F["two_exon"] = Feature(gene_name="twoN", gene_type="CDS", start=2900, end=2973, strand=1,
                        engine="AB", flag="HIGH")
F["two_exon"].exons = [(2900, 2914), (2960, 2973)]
trna = Feature(gene_name="trnX-UUU", gene_type="tRNA", start=3100, end=3172, strand=1,
               engine="B", flag="HIGH")
trna.exons = [(3100, 3172)]
place(3100, "TAA" * 24)
SEQ = "".join(GENOME)
ALL = list(F.values()) + [trna]
geometry = {id(f): (f.start, f.end, f.strand, list(f.exons)) for f in ALL}

print("--- one pass ---")
out = submission_check(ALL, SEQ)
flagged = {f.gene_name: probs for f, probs in out}
note = lambda f: next((str(n) for n in f.notes if str(n).startswith(SUBMIT)), "")

check("a valid ATG ... TAA CDS passes, flag kept", (F["valid"].flag, note(F["valid"])), ("HIGH", ""))
check("GTG is a table-11 start and passes", F["gtg"].gene_name in flagged, False)
check("a CCC start is flagged", flagged.get("badC"), ["no valid start codon (CCC)"])
check("  ... and becomes NEEDS_REVIEW", F["bad_start"].flag, "NEEDS_REVIEW")
check("  ... and the note keeps the flag it had",
      note(F["bad_start"]).endswith("flag lowered from HIGH"), True)
check("an ACG start written with the RNA-editing exception passes",
      F["acg_edited"].gene_name in flagged, False)
check("an ACG start without it is flagged", flagged.get("acgD"), ["no valid start codon (ACG)"])
check("a missing stop codon is flagged", flagged.get("nostopE"),
      ["no stop codon at the 3' end (GCC)"])
check("a 3' end that is not a whole codon is flagged as no stop",
      flagged.get("frameF"), ["no stop codon at the 3' end (AAG)"])
check("an internal stop codon is flagged", flagged.get("intG"), ["1 internal stop codon(s)"])
check("a valid minus-strand CDS passes (read on its own strand)",
      F["minus_ok"].gene_name in flagged, False)
check("a minus-strand CCC start is flagged", flagged.get("minusI"),
      ["no valid start codon (CCC)"])
check("a partial CDS is exempt from the start and stop checks",
      F["partial"].gene_name in flagged, False)
check("  ... but not from the internal-stop check", flagged.get("partK"),
      ["1 internal stop codon(s)"])
check("a pseudogene is skipped", (F["pseudo"].gene_name in flagged, F["pseudo"].flag),
      (False, "HIGH"))
check("a CDS already NEEDS_REVIEW is listed, its note claims no lowering",
      (F["already"].gene_name in flagged, "lowered" in note(F["already"])), (True, False))
check("a two-exon CDS is checked on its spliced sequence, not its span",
      F["two_exon"].gene_name in flagged, False)
check("a tRNA is never checked", (trna.gene_name in flagged, trna.flag), (False, "HIGH"))
check("no coordinate, strand or exon moved",
      all(geometry[id(f)] == (f.start, f.end, f.strand, list(f.exons)) for f in ALL), True)
check("confidence untouched", all(f.confidence == 0.0 for f in ALL), True)

print("--- idempotent ---")
snap = [(f.flag, list(f.notes)) for f in ALL]
out2 = submission_check(ALL, SEQ)
check("a second pass changes no flag and no note",
      [(f.flag, list(f.notes)) for f in ALL] == snap, True)
check("  ... and returns the same features",
      [f.gene_name for f, _ in out2] == [f.gene_name for f, _ in out], True)
check("no note is duplicated",
      all(sum(str(n).startswith(SUBMIT) for n in f.notes) <= 1 for f in ALL), True)

print("--- after finalize_qc, twice, as the pipeline would ---")
fresh = [copy.deepcopy(f) for f in ALL]
for f in fresh:
    f.flag = "HIGH" if f.gene_name != "revM" else "NEEDS_REVIEW"
    f.notes = []
passes = []
for _ in range(2):
    finalize_qc(fresh, SEQ)
    submission_check(fresh, SEQ)
    passes.append([(f.flag, [str(n) for n in f.notes]) for f in fresh])
check("finalize_qc + submission_check twice: identical flags and notes",
      passes[0] == passes[1], True)

print("--- where it runs ---")
src = inspect.getsource(sys.modules["plastanno.pipeline"]) if "plastanno.pipeline" in sys.modules \
    else open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "plastanno", "pipeline.py")).read()
i_revoke = src.find("revoke_implausible_rescues(annotations")
i_check = src.find("submission_check(annotations")
i_write = src.find("write_all(")
check("the pipeline calls submission_check", i_check > 0, True)
check("  ... after the rescue gate, so its flags cannot revoke a feature",
      0 < i_revoke < i_check, True)
check("  ... and before any file is written", 0 < i_check < i_write, True)

print("--- the report ---")
with tempfile.TemporaryDirectory() as d:
    p = os.path.join(d, "x.report")
    write_report(ALL, "TEST", len(SEQ), {}, [], 1.0, p)
    rep = open(p).read()
    check("the report has a Submission check section", "Submission check" in rep, True)
    check("  ... counting every failing CDS",
          "%d CDS would fail NCBI validation" % len(out) in rep, True)
    check("  ... and naming them", all(g in rep for g in flagged), True)
    clean = [copy.deepcopy(F["valid"])]
    clean[0].notes = []
    write_report(clean, "TEST", len(SEQ), {}, [], 1.0, p)
    check("a clean annotation says so", "No CDS would fail NCBI validation" in open(p).read(), True)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
