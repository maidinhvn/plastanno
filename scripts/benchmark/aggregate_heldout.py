#!/usr/bin/env python3
"""
Aggregate the 10 held-out chunk results into the final generalisation tables.

Read-only: consumes bench_runs/heldout/chunk_*/results.json produced by the
run-once held-out evaluation (scripts/run_heldout.sh). Reports:
  (a) pooled micro F1 over all genomes, by gene type + global;
  (b) per-batch global F1 mean +/- SD across the 10 chunks (stability);
  (c) per-genome macro F1;
  (d) F1 by structural mode;
  (e) F1 by major taxonomic family.
"""
import os, re, json, glob
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HELD = os.path.join(REPO, "bench_runs", "heldout")
CLS  = os.environ.get("PLASTANNO_DATA","benchmark_data")+"/classify_results.json"
FAM  = os.path.join(REPO, "splits", "heldout_families.tsv")
strip = lambda a: re.sub(r"\.\d+$", "", a)


def prf(tp, ref, pred):
    sn = tp / ref if ref else 0.0
    pr = tp / pred if pred else 0.0
    f1 = 2 * sn * pr / (sn + pr) if (sn + pr) else 0.0
    return sn * 100, pr * 100, f1 * 100


def main():
    chunks = sorted(glob.glob(os.path.join(HELD, "chunk_*", "results.json")))
    if not chunks:
        raise SystemExit("No chunk results under %s (run scripts/run_heldout.sh first)" % HELD)

    modes = {strip(k): v.get("mode", "UNKNOWN")
             for k, v in json.load(open(CLS)).items()}
    fam = {}
    if os.path.exists(FAM):
        for line in open(FAM):
            a, f = line.rstrip("\n").split("\t"); fam[a] = f

    allr, per_batch = [], []
    for cp in chunks:
        res = [r for r in json.load(open(cp))["results"] if r.get("status") == "OK"]
        allr += res
        tp = sum(r["tp"] for r in res); ref = sum(r["ref"] for r in res)
        pred = sum(r["pred"] for r in res)
        per_batch.append(prf(tp, ref, pred)[2])

    n = len(allr)
    print("=" * 64)
    print("HELD-OUT EVALUATION — %d genomes, %d chunks (FROZEN tool, run once)"
          % (n, len(chunks)))
    print("=" * 64)

    print("\nPooled (micro):")
    for t in ("CDS", "tRNA", "rRNA"):
        tp = sum(r["tp_by_type"].get(t, 0) for r in allr)
        ref = sum(r["ref_by_type"].get(t, 0) for r in allr)
        pred = sum(r["pred_by_type"].get(t, 0) for r in allr)
        sn, pr, f1 = prf(tp, ref, pred)
        print("  %-6s TP=%-6d ref=%-6d pred=%-6d  Sn=%.1f%%  Pr=%.1f%%  F1=%.1f%%"
              % (t, tp, ref, pred, sn, pr, f1))
    gtp = sum(r["tp"] for r in allr); gref = sum(r["ref"] for r in allr)
    gpred = sum(r["pred"] for r in allr)
    sn, pr, f1 = prf(gtp, gref, gpred)
    print("  GLOBAL TP=%-6d                     Sn=%.1f%%  Pr=%.1f%%  F1=%.1f%%"
          % (gtp, sn, pr, f1))

    mean = sum(per_batch) / len(per_batch)
    sd = (sum((x - mean) ** 2 for x in per_batch) / len(per_batch)) ** 0.5
    print("\nPer-batch global F1: mean=%.1f%%  SD=%.1f%%  (%d batches)"
          % (mean, sd, len(per_batch)))
    print("  batch values:", ["%.1f" % x for x in per_batch])

    pg = sorted(prf(r["tp"], r["ref"], r["pred"])[2] for r in allr)
    print("Per-genome macro F1: mean=%.1f%%  median=%.1f%%  min=%.1f%%  max=%.1f%%"
          % (sum(pg) / len(pg), pg[len(pg) // 2], pg[0], pg[-1]))

    print("\nBy structural mode:")
    bym = defaultdict(lambda: [0, 0, 0, 0])
    for r in allr:
        d = bym[modes.get(strip(r["acc"]), "UNKNOWN")]
        d[0] += r["tp"]; d[1] += r["ref"]; d[2] += r["pred"]; d[3] += 1
    for m, d in sorted(bym.items(), key=lambda x: -x[1][3]):
        print("  %-12s n=%-5d F1=%.1f%%" % (m, d[3], prf(d[0], d[1], d[2])[2]))

    print("\nTop 12 families by count:")
    byf = defaultdict(lambda: [0, 0, 0, 0])
    for r in allr:
        d = byf[fam.get(r["acc"], "?")]
        d[0] += r["tp"]; d[1] += r["ref"]; d[2] += r["pred"]; d[3] += 1
    for f, d in sorted(byf.items(), key=lambda x: -x[1][3])[:12]:
        print("  %-20s n=%-5d F1=%.1f%%" % (f, d[3], prf(d[0], d[1], d[2])[2]))

    failed = sum(1 for cp in chunks
                 for r in json.load(open(cp))["results"] if r.get("status") != "OK")
    if failed:
        print("\n[!] %d genomes did NOT return OK — inspect chunk logs." % failed)


if __name__ == "__main__":
    main()
