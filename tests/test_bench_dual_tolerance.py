#!/usr/bin/env python3
"""multi_genome_bench.py must report +/-0 beside +/-60, and must not disturb +/-60.

The published figure runs at `--tol 60`, where a call one base out of place scores
exactly like an exact one. On the development set that hid a defect affecting 63%
of tRNA loci behind a tRNA F1 of 93.6%, and then hid the fix as well: the whole
--exon-mode / --trna-mode programme moved ONE locus at +/-60.

So the aggregator now reports both. Two things have to hold:

  1. the +/-60 numbers are byte-for-byte what they were, so every existing
     bench_runs/*.log stays reproducible and comparable;
  2. an old result dict, written before +/-0 existed, still aggregates — a
     benchmark that crashes on its own history is not a benchmark.
"""
import importlib.util
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

spec = importlib.util.spec_from_file_location(
    "mgb", os.path.join(REPO, "scripts", "benchmark", "multi_genome_bench.py"))
mgb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mgb)

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-66s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


def has(label, text, needle, want=True):
    RUN[0] += 1
    got = needle in text
    ok = got == want
    print("  %-66s %s" % (label, "ok" if ok else "FAIL %r present=%s" % (needle, got)))
    if not ok:
        FAIL.append(label)


def result(acc, tp, fp, fn, tp0=None, fp0=None, fn0=None):
    r = {
        "acc": acc, "status": "OK",
        "tp": tp, "fp": fp, "fn": fn, "ref": tp + fn, "pred": tp + fp,
        "tp_by_type":   {"CDS": tp, "tRNA": 0, "rRNA": 0},
        "ref_by_type":  {"CDS": tp + fn, "tRNA": 0, "rRNA": 0},
        "pred_by_type": {"CDS": tp + fp, "tRNA": 0, "rRNA": 0},
    }
    if tp0 is not None:
        r.update({"tp_exact": tp0, "fp_exact": fp0, "fn_exact": fn0,
                  "tp_by_type_exact": {"CDS": tp0, "tRNA": 0, "rRNA": 0}})
    return r


print("--- 1. both tolerances are reported, with the arithmetic stated ---")
# 90 TP at +/-60, of which only 50 are exact; the other 40 become FN at +/-0.
new = [result("A", 90, 10, 10, 50, 50, 50)]
txt = mgb.aggregate(new, {})
has("the +/-60 line is still there", txt, "GLOBAL (micro, pooled): TP=90 FP=10 FN=10")
has("a +/-0 line appears", txt, "GLOBAL at +/-0    (exact): TP=50 FP=50 FN=50")
has("the boundary cost is named", txt, "boundary cost:")
has("it counts the inexact matches", txt, "(40 of 90 matches are not exact)")
has("no spurious 'predate' note when every genome has it", txt, "predate", want=False)

# F1 at 60 = 2*.9*.9/1.8 = 90.0 ; at 0 = 2*.5*.5/1.0 = 50.0 ; cost 40.0 points
has("the cost is the F1 difference, not invented", txt, "40.0 F1 points")

print("--- 2. +/-0 per-type is a SEPARATE block, so the +/-60 rows are untouched ---")
has("the +/-60 CDS row is exactly as it always was",
    txt, "  CDS   TP=90    ref=100   pred=100    Sens=90.0% Prec=90.0% F1=90.0%")
has("a separate +/-0 block exists", txt, "By gene type at +/-0 (exact boundaries):")
has("with the exact TP", txt, "CDS   TP=50")
has("and the share of matches that are exact", txt, "exact/matched 55.6%")

print("--- 3. an old result dict, with no +/-0 keys, still aggregates ---")
old = [result("A", 90, 10, 10)]
txt_old = mgb.aggregate(old, {})
has("the +/-60 block is identical to before",
    txt_old, "GLOBAL (micro, pooled): TP=90 FP=10 FN=10")
has("no +/-0 block is fabricated", txt_old, "GLOBAL at +/-0", want=False)
has("and the per-type row has no +/-0 tail", txt_old, "+/-0: TP=", want=False)

print("--- 4. the +/-60 text is byte-identical with and without the new keys ---")
def sixty_block(t):
    """Everything except the +/-0 additions: the two GLOBAL/note lines and the
    whole per-type block, header and trailing blank included."""
    out, skip = [], False
    for line in t.splitlines():
        if line.startswith("By gene type at +/-0"):
            skip = True
            if out and out[-1] == "":
                out.pop()                 # the blank line that introduced it
            continue
        if skip:
            if line.strip() == "":
                skip = False
                out.append(line)      # the blank belongs to the NEXT section
            continue
        if "+/-0" in line or "boundary cost" in line or "predate" in line:
            continue
        out.append(line)
    return "\n".join(out)
check("adding the +/-0 keys does not perturb one +/-60 character",
      sixty_block(txt), sixty_block(txt_old))

print("--- 5. a mixed batch says so instead of quietly averaging ---")
mixed = [result("A", 90, 10, 10, 50, 50, 50), result("B", 80, 20, 20)]
txt_mix = mgb.aggregate(mixed, {})
has("the note names how many genomes lack the measurement",
    txt_mix, "1 of 2 genomes predate")
has("and the +/-0 pooled figure uses only the genomes that have it",
    txt_mix, "GLOBAL at +/-0    (exact): TP=50 FP=50 FN=50")

print("--- 6. +/-0 really is stricter, on a real prediction/reference pair ---")
bm_path = os.path.join(REPO, "scripts", "benchmark", "benchmark_gene_by_gene.py")
s2 = importlib.util.spec_from_file_location("bm", bm_path)
bm = importlib.util.module_from_spec(s2)
s2.loader.exec_module(bm)

pred_gb = os.path.join(REPO, "benchmark_v3", "gate_b", "aragorn_e2e",
                       "E_aragorn", "NC_019616.1", "NC_019616.1.gb")
ref_gb = os.path.join(os.environ.get("PLASTANNO_DATA", "benchmark_data"),
                      "eval_v2_raw", "NC_019616.1.gb")
if os.path.exists(pred_gb) and os.path.exists(ref_gb):
    ref = bm.load_units(ref_gb)
    pred = bm.load_units(pred_gb)
    tp60 = len(bm.run(ref, pred, 60, 0.6)[0])
    tp0 = len(bm.run(ref, pred, 0, 0.6)[0])
    check("tol=0 accepts no more than tol=60", tp0 <= tp60, True)
    check("and strictly fewer on a genome with known inexact boundaries",
          tp0 < tp60, True)
    print("       (%d exact of %d matched at +/-60)" % (tp0, tp60))
else:
    print("  SKIPPED: fixture missing (%s)" % ("pred" if not os.path.exists(pred_gb) else "ref"))

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
