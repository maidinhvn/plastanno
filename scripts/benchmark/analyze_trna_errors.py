#!/usr/bin/env python3
"""
tRNA error decomposition across a finished multi_genome_bench run.

Mirrors analyze_cds_errors.py but for tRNA, classifying each error as a
detection failure, a boundary error, or an anticodon/amino-acid naming mismatch,
and tabulating the most common ref->pred confusion pairs (e.g. trnM-CAU vs
trnfM-CAU vs trnI-CAU, or the trnS family).

  False negatives (missed reference tRNA):
    true_miss   - no predicted tRNA overlaps it                       (detection)
    boundary    - a name-matching pred overlaps but ends >tol apart   (boundary)
    naming      - only a DIFFERENT-name pred overlaps it              (anticodon/aa)

  False positives (spurious predicted tRNA): symmetric (spurious = no overlap).

Usage:
    python analyze_trna_errors.py <bench_run_dir> [--rawdata DIR] [--tol 60]
"""
import argparse
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from benchmark_gene_by_gene import (
    load_units, name_match, coord_close, overlaps, run as bench_run,
)

DEFAULT_RAW = os.environ.get("PLASTANNO_DATA","/data06/users/vutrinh/Apiales_Plastomes_20260514")+"/rawdata"
_TOL = 60
_RUN = None
_RAW = None


def _classify_fn(r, pred_t, tol):
    same = [p for p in pred_t if name_match(r["name"], p["name"]) and overlaps(r, p)]
    if same:
        return "boundary", None        # name matches + overlaps but missed tol/sim
    diff = [p for p in pred_t if overlaps(r, p)]
    if diff:
        return "naming", diff[0]["name"]
    return "true_miss", None


def _classify_fp(p, ref_t, tol):
    same = [r for r in ref_t if name_match(r["name"], p["name"]) and overlaps(r, p)]
    if same:
        return "boundary", None
    diff = [r for r in ref_t if overlaps(r, p)]
    if diff:
        return "naming", diff[0]["name"]
    return "spurious", None


def worker(acc):
    pred_gb = os.path.join(_RUN, acc, acc + ".gb")
    ref_gb = os.path.join(_RAW, acc + ".gb")
    if not (os.path.exists(pred_gb) and os.path.exists(ref_gb)):
        return None
    try:
        ref = load_units(ref_gb)
        pred = load_units(pred_gb)
    except Exception:
        return None
    tp, fp, fn, detect, mism = bench_run(ref, pred, _TOL, 0.6)
    ref_t = [u for u in ref if u["type"] == "tRNA"]
    pred_t = [u for u in pred if u["type"] == "tRNA"]
    fnk, fpk, fng, fpg, conf = Counter(), Counter(), Counter(), Counter(), Counter()
    for r in (u for u in fn if u["type"] == "tRNA"):
        kind, info = _classify_fn(r, pred_t, _TOL)
        fnk[kind] += 1
        fng[r["name"]] += 1
        if kind == "naming":
            conf[(r["name"], info)] += 1
    for p in (u for u in fp if u["type"] == "tRNA"):
        kind, info = _classify_fp(p, ref_t, _TOL)
        fpk[kind] += 1
        fpg[p["name"]] += 1
    return fnk, fpk, fng, fpg, conf


def main():
    import multiprocessing as mp
    global _RUN, _RAW, _TOL
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--rawdata", default=DEFAULT_RAW)
    ap.add_argument("--tol", type=int, default=60)
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()
    _RUN, _RAW, _TOL = a.run_dir, a.rawdata, a.tol

    fn_kinds, fp_kinds = Counter(), Counter()
    fn_gene, fp_gene, confusions = Counter(), Counter(), Counter()
    n = 0
    accs = sorted(os.listdir(a.run_dir))
    with mp.Pool(a.workers) as pool:
        for res in pool.imap_unordered(worker, accs):
            if res is None:
                continue
            fnk, fpk, fng, fpg, conf = res
            n += 1
            fn_kinds.update(fnk); fp_kinds.update(fpk)
            fn_gene.update(fng); fp_gene.update(fpg)
            confusions.update(conf)

    tot_fn = sum(fn_kinds.values()); tot_fp = sum(fp_kinds.values())
    print("=" * 64)
    print("tRNA ERROR DECOMPOSITION  (%d genomes, tol=+/-%dbp)" % (n, a.tol))
    print("=" * 64)
    print("\nFALSE NEGATIVES (missed tRNA): %d" % tot_fn)
    for k, v in fn_kinds.most_common():
        print("  %-10s %5d  (%4.1f%%)" % (k, v, 100 * v / tot_fn if tot_fn else 0))
    print("\nFALSE POSITIVES (spurious tRNA): %d" % tot_fp)
    for k, v in fp_kinds.most_common():
        print("  %-10s %5d  (%4.1f%%)" % (k, v, 100 * v / tot_fp if tot_fp else 0))
    print("\nTop tRNA by FALSE NEGATIVE count:")
    for g, v in fn_gene.most_common(15):
        print("  %-14s %3d" % (g, v))
    print("\nTop tRNA by FALSE POSITIVE count:")
    for g, v in fp_gene.most_common(15):
        print("  %-14s %3d" % (g, v))
    print("\nTop anticodon/aa confusion pairs (ref -> pred, overlapping diff name):")
    for (rn, pn), v in confusions.most_common(15):
        print("  %-14s -> %-14s %3d" % (rn, pn, v))


if __name__ == "__main__":
    main()
