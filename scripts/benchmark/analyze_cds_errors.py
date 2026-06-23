#!/usr/bin/env python3
"""
CDS error decomposition across a finished multi_genome_bench run.

For every genome in a bench run dir, re-loads the reference and predicted
GenBank, reproduces the exact TP/FP/FN matching used by
benchmark_gene_by_gene.py (same name_match / coord_close / similarity), then
classifies CDS errors so we can see WHERE the CDS gap comes from:

  False negatives (missed reference CDS):
    true_miss   - no predicted CDS of the same name overlaps it          (detection failure)
    boundary    - a same-name pred overlaps but ends are >tol apart        (boundary error)
    low_sim     - a same-name pred overlaps, ends within tol, but seq sim < threshold
    naming      - only a DIFFERENT-name pred overlaps it                   (synonym / mis-call)

  False positives (spurious predicted CDS):
    boundary    - overlaps a same-name ref (the partner of a boundary FN)
    low_sim     - overlaps a same-name ref, ends within tol, seq sim < threshold
    naming      - overlaps a different-name ref
    spurious    - overlaps no ref CDS at all                               (pure false ORF)

Usage:
    python analyze_cds_errors.py <bench_run_dir> [--rawdata DIR] [--tol 60] [--sim 0.6]
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from benchmark_gene_by_gene import (
    load_units, name_match, coord_close, overlaps, similarity, run as bench_run,
)

DEFAULT_RAW = os.environ.get("PLASTANNO_DATA","/data06/users/vutrinh/Apiales_Plastomes_20260514")+"/rawdata"


def classify_fn(r, pred_cds, tol, sim):
    """Classify one reference-CDS false negative."""
    same_overlap = [p for p in pred_cds if name_match(r["name"], p["name"]) and overlaps(r, p)]
    if same_overlap:
        # pick the best-overlapping same-name pred to judge boundary vs seq
        p = max(same_overlap, key=lambda p: min(r["end"], p["end"]) - max(r["start"], p["start"]))
        if not coord_close(r, p, tol):
            d = abs(r["start"] - p["start"]) + abs(r["end"] - p["end"])
            return "boundary", d
        return "low_sim", round(similarity(r, p), 2)
    diff_overlap = [p for p in pred_cds if overlaps(r, p)]
    if diff_overlap:
        return "naming", diff_overlap[0]["name"]
    return "true_miss", None


def classify_fp(p, ref_cds, tol, sim):
    """Classify one predicted-CDS false positive."""
    same_overlap = [r for r in ref_cds if name_match(r["name"], p["name"]) and overlaps(r, p)]
    if same_overlap:
        r = max(same_overlap, key=lambda r: min(r["end"], p["end"]) - max(r["start"], p["start"]))
        if not coord_close(r, p, tol):
            return "boundary", None
        return "low_sim", None
    diff_overlap = [r for r in ref_cds if overlaps(r, p)]
    if diff_overlap:
        return "naming", diff_overlap[0]["name"]
    return "spurious", None


def process_one(args):
    """Worker: classify CDS errors for one genome. Returns partial dicts or None."""
    acc, run_dir, rawdata, tol, sim = args
    pred_gb = os.path.join(run_dir, acc, acc + ".gb")
    ref_gb = os.path.join(rawdata, acc + ".gb")
    if not (os.path.exists(pred_gb) and os.path.exists(ref_gb)):
        return None
    try:
        ref = load_units(ref_gb)
        pred = load_units(pred_gb)
    except Exception:
        return None

    tp, fp, fn, detect, mism = bench_run(ref, pred, tol, sim)
    ref_cds = [u for u in ref if u["type"] == "CDS"]
    pred_cds = [u for u in pred if u["type"] == "CDS"]

    fnk, fpk = Counter(), Counter()
    fng, fpg = defaultdict(Counter), defaultdict(Counter)
    bdeltas, lsims = [], []

    for r in (u for u in fn if u["type"] == "CDS"):
        kind, info = classify_fn(r, pred_cds, tol, sim)
        fnk[kind] += 1
        fng[r["name"]][kind] += 1
        if kind == "boundary":
            bdeltas.append(info)
        elif kind == "low_sim":
            lsims.append(info)
    for p in (u for u in fp if u["type"] == "CDS"):
        kind, info = classify_fp(p, ref_cds, tol, sim)
        fpk[kind] += 1
        fpg[p["name"]][kind] += 1
    return fnk, fpk, fng, fpg, bdeltas, lsims


def main():
    import multiprocessing as mp
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--rawdata", default=DEFAULT_RAW)
    ap.add_argument("--tol", type=int, default=60)
    ap.add_argument("--sim", type=float, default=0.6)
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()

    fn_kinds = Counter()
    fp_kinds = Counter()
    fn_gene = defaultdict(Counter)   # gene -> kind -> count
    fp_gene = defaultdict(Counter)
    boundary_deltas = []
    low_sim_vals = []
    n_genomes = 0

    accs = sorted(os.listdir(a.run_dir))
    jobs = [(acc, a.run_dir, a.rawdata, a.tol, a.sim) for acc in accs]
    with mp.Pool(a.workers) as pool:
        for res in pool.imap_unordered(process_one, jobs):
            if res is None:
                continue
            fnk, fpk, fng, fpg, bdeltas, lsims = res
            n_genomes += 1
            fn_kinds.update(fnk)
            fp_kinds.update(fpk)
            for g, c in fng.items():
                fn_gene[g].update(c)
            for g, c in fpg.items():
                fp_gene[g].update(c)
            boundary_deltas.extend(bdeltas)
            low_sim_vals.extend(lsims)

    print("=" * 64)
    print("CDS ERROR DECOMPOSITION  (%d genomes, tol=+/-%dbp, sim>=%.2f)"
          % (n_genomes, a.tol, a.sim))
    print("=" * 64)

    tot_fn = sum(fn_kinds.values())
    tot_fp = sum(fp_kinds.values())
    print("\nFALSE NEGATIVES (missed CDS): %d" % tot_fn)
    for k, v in fn_kinds.most_common():
        print("  %-10s %5d  (%4.1f%%)" % (k, v, 100 * v / tot_fn if tot_fn else 0))
    print("\nFALSE POSITIVES (spurious CDS): %d" % tot_fp)
    for k, v in fp_kinds.most_common():
        print("  %-10s %5d  (%4.1f%%)" % (k, v, 100 * v / tot_fp if tot_fp else 0))

    if boundary_deltas:
        boundary_deltas.sort()
        n = len(boundary_deltas)
        print("\nBoundary FN end-offset (|Δstart|+|Δend| bp): "
              "median=%d  p25=%d  p75=%d  max=%d"
              % (boundary_deltas[n // 2], boundary_deltas[n // 4],
                 boundary_deltas[3 * n // 4], boundary_deltas[-1]))

    print("\nTop genes by FALSE NEGATIVE count:")
    fn_rank = sorted(fn_gene.items(), key=lambda kv: -sum(kv[1].values()))[:20]
    for gene, kinds in fn_rank:
        tot = sum(kinds.values())
        detail = " ".join("%s=%d" % (k, v) for k, v in kinds.most_common())
        print("  %-12s %3d   (%s)" % (gene, tot, detail))

    print("\nTop genes by FALSE POSITIVE count:")
    fp_rank = sorted(fp_gene.items(), key=lambda kv: -sum(kv[1].values()))[:20]
    for gene, kinds in fp_rank:
        tot = sum(kinds.values())
        detail = " ".join("%s=%d" % (k, v) for k, v in kinds.most_common())
        print("  %-12s %3d   (%s)" % (gene, tot, detail))


if __name__ == "__main__":
    main()
