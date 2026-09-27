#!/usr/bin/env python3
"""
Multi-genome benchmark for Plastanno over a stratified sample of the DEV set.

For each sampled accession: run the pipeline into an isolated output directory
(parallel ACROSS genomes, threads=1 each, since Exonerate does not scale with
threads), then score the predicted .gb against the rawdata reference .gb using
the SAME gene-by-gene benchmark. Aggregates F1 globally (micro), per genome
(macro mean/std), per gene type, and per structure mode.

The HELD-OUT set is guarded: it can only be read with --set heldout --final, so
it is never touched during ordinary development.

Usage:
  python multi_genome_bench.py --n 120 --workers 16            # sample DEV
  python multi_genome_bench.py --n 12 --workers 8 --keep       # quick smoke test
"""
import argparse, json, os, re, random, subprocess, importlib.util, statistics, datetime
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

REPO   = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CLS    = os.environ.get("PLASTANNO_DATA","benchmark_data")+"/classify_results.json"
RAW    = os.environ.get("PLASTANNO_DATA","benchmark_data")+"/rawdata"
SPLITS = os.path.join(REPO, "splits")

strip = lambda a: re.sub(r"\.\d+$", "", a)


def _bm():
    p = os.path.join(REPO, "scripts/benchmark/benchmark_gene_by_gene.py")
    spec = importlib.util.spec_from_file_location("bm", p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def load_modes():
    cls = json.load(open(CLS))
    return {strip(k): v.get("mode", "UNKNOWN") for k, v in cls.items()}


def sample(set_name, n, seed, modes):
    accs = [l.strip() for l in open(os.path.join(SPLITS, set_name + "_set.txt")) if l.strip()]
    by = defaultdict(list)
    for a in accs:
        by[modes.get(strip(a), "UNKNOWN")].append(a)
    rng = random.Random(seed)
    total = len(accs)
    chosen = []
    for m in sorted(by):
        grp = sorted(by[m]); rng.shuffle(grp)
        k = max(1, round(n * len(grp) / total)) if grp else 0
        chosen += grp[:k]
    rng.shuffle(chosen)
    return chosen


def run_one(acc, workdir, bm):
    fasta  = os.path.join(RAW, acc + ".fasta")
    ref_gb = os.path.join(RAW, acc + ".gb")
    out    = os.path.join(workdir, acc)
    if not (os.path.exists(fasta) and os.path.exists(ref_gb)):
        return {"acc": acc, "status": "NO_INPUT"}
    try:
        r = subprocess.run(
            ["python3", "plastanno.py", "run", fasta, "-o", out, "-t", "1", "--no-plot"],
            capture_output=True, text=True, timeout=900, cwd=REPO,
        )
    except subprocess.TimeoutExpired:
        return {"acc": acc, "status": "TIMEOUT"}
    pred_gb = os.path.join(out, acc + ".gb")
    if not os.path.exists(pred_gb):
        return {"acc": acc, "status": "NO_OUTPUT", "err": (r.stderr or "")[-200:]}
    try:
        ref  = bm.load_units(ref_gb)
        pred = bm.load_units(pred_gb)
        # Scored at TWO tolerances. +/-60 is the published figure and its keys are
        # untouched, so every number in bench_runs/*.log stays reproducible.
        #
        # +/-0 is added because +/-60 cannot see boundary quality at all: a call
        # that is one base out scores identically to an exact one. Measured on the
        # development set, 63% of tRNA had inexact boundaries while tRNA F1 read
        # 93.6%, and the whole --exon-mode/--trna-mode programme moved ONE locus
        # at +/-60. A benchmark that cannot see the defect cannot see the fix
        # either, in either direction.
        tp, fp, fn, detect, mism = bm.run(ref, pred, 60, 0.6)
        tp0, fp0, fn0, _d0, _m0 = bm.run(ref, pred, 0, 0.6)
        tp_t = Counter(rr["type"] for _p, rr, _s in tp)
        tp_t0 = Counter(rr["type"] for _p, rr, _s in tp0)
        return {
            "acc": acc, "status": "OK",
            "tp": len(tp), "fp": len(fp), "fn": len(fn),
            "ref": len(ref), "pred": len(pred),
            "tp_by_type":  dict(tp_t),
            "ref_by_type": dict(Counter(u["type"] for u in ref)),
            "pred_by_type": dict(Counter(u["type"] for u in pred)),
            "tp_exact": len(tp0), "fp_exact": len(fp0), "fn_exact": len(fn0),
            "tp_by_type_exact": dict(tp_t0),
        }
    except Exception as e:
        return {"acc": acc, "status": "BENCH_ERROR", "err": str(e)[:200]}


def f1(tp, fp, fn):
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    f = 2 * sens * prec / (sens + prec) if (sens + prec) else 0.0
    return sens, prec, f


def aggregate(results, modes):
    ok = [r for r in results if r["status"] == "OK"]
    bad = [r for r in results if r["status"] != "OK"]
    lines = []
    lines.append("Ran: %d | OK: %d | failed: %d %s" % (
        len(results), len(ok), len(bad),
        dict(Counter(r["status"] for r in bad)) if bad else ""))
    if not ok:
        return "\n".join(lines)

    TP = sum(r["tp"] for r in ok); FP = sum(r["fp"] for r in ok); FN = sum(r["fn"] for r in ok)
    s, p, f = f1(TP, FP, FN)
    lines.append("")
    lines.append("GLOBAL (micro, pooled): TP=%d FP=%d FN=%d  Sens=%.1f%% Prec=%.1f%% F1=%.1f%%"
                 % (TP, FP, FN, s*100, p*100, f*100))
    per = [f1(r["tp"], r["fp"], r["fn"])[2] for r in ok]
    lines.append("Per-genome F1 (macro): mean=%.1f%% median=%.1f%% std=%.1f%% min=%.1f%% max=%.1f%%"
                 % (statistics.mean(per)*100, statistics.median(per)*100,
                    (statistics.pstdev(per)*100 if len(per) > 1 else 0.0),
                    min(per)*100, max(per)*100))

    # The same pooled figure at zero tolerance. Reported beside the headline, not
    # instead of it: the gap between the two IS the boundary-quality measurement.
    if any("tp_exact" in r for r in ok):
        oke = [r for r in ok if "tp_exact" in r]
        TP0 = sum(r["tp_exact"] for r in oke)
        FP0 = sum(r["fp_exact"] for r in oke)
        FN0 = sum(r["fn_exact"] for r in oke)
        s0, p0, f0 = f1(TP0, FP0, FN0)
        lines.append("GLOBAL at +/-0    (exact): TP=%d FP=%d FN=%d  Sens=%.1f%% Prec=%.1f%% F1=%.1f%%"
                     % (TP0, FP0, FN0, s0*100, p0*100, f0*100))
        lines.append("  boundary cost: %.1f F1 points are spent on inexact boundaries"
                     " (%d of %d matches are not exact)"
                     % ((f - f0)*100, TP - TP0, TP))
        if len(oke) != len(ok):
            lines.append("  NOTE: %d of %d genomes predate the +/-0 measurement and are"
                         " excluded from it" % (len(ok) - len(oke), len(ok)))

    lines.append("")
    lines.append("By gene type (pooled):")
    for t in ("CDS", "tRNA", "rRNA"):
        tp_t = sum(r["tp_by_type"].get(t, 0) for r in ok)
        rf_t = sum(r["ref_by_type"].get(t, 0) for r in ok)
        pr_t = sum(r["pred_by_type"].get(t, 0) for r in ok)
        fp_t = pr_t - tp_t; fn_t = rf_t - tp_t
        s, p, ff = f1(tp_t, fp_t, fn_t)
        lines.append("  %-5s TP=%-5d ref=%-5d pred=%-5d  Sens=%.1f%% Prec=%.1f%% F1=%.1f%%"
                     % (t, tp_t, rf_t, pr_t, s*100, p*100, ff*100))

    # A separate block rather than a tail on the rows above, so the +/-60 section
    # stays character-for-character what it has always been and an old log can be
    # diffed against a new one.
    if any("tp_by_type_exact" in r for r in ok):
        lines.append("")
        lines.append("By gene type at +/-0 (exact boundaries):")
        for t in ("CDS", "tRNA", "rRNA"):
            rf_t = sum(r["ref_by_type"].get(t, 0) for r in ok)
            pr_t = sum(r["pred_by_type"].get(t, 0) for r in ok)
            tp_t = sum(r["tp_by_type"].get(t, 0) for r in ok)
            tp_e = sum(r.get("tp_by_type_exact", {}).get(t, 0) for r in ok)
            _s, _p, ff = f1(tp_t, pr_t - tp_t, rf_t - tp_t)
            se, pe, fe = f1(tp_e, pr_t - tp_e, rf_t - tp_e)
            lines.append("  %-5s TP=%-5d  Sens=%.1f%% Prec=%.1f%% F1=%.1f%%"
                         "   boundary cost %.1f pts   exact/matched %.1f%%"
                         % (t, tp_e, se*100, pe*100, fe*100, (ff - fe)*100,
                            100*tp_e/tp_t if tp_t else 0.0))

    lines.append("")
    lines.append("By structure mode (pooled):")
    bym = defaultdict(list)
    for r in ok:
        bym[modes.get(strip(r["acc"]), "UNKNOWN")].append(r)
    for m in sorted(bym, key=lambda x: -len(bym[x])):
        g = bym[m]
        TPm = sum(r["tp"] for r in g); FPm = sum(r["fp"] for r in g); FNm = sum(r["fn"] for r in g)
        s, p, ff = f1(TPm, FPm, FNm)
        note = "  (n nhỏ — diễn giải thận trọng)" if len(g) < 10 else ""
        lines.append("  %-11s n=%-4d TP=%-5d FP=%-4d FN=%-4d  F1=%.1f%%%s"
                     % (m, len(g), TPm, FPm, FNm, ff*100, note))

    worst = sorted(ok, key=lambda r: f1(r["tp"], r["fp"], r["fn"])[2])[:10]
    lines.append("")
    lines.append("10 genome F1 thấp nhất:")
    for r in worst:
        _s, _p, ff = f1(r["tp"], r["fp"], r["fn"])
        lines.append("  %-14s %-11s F1=%.1f%%  TP=%d FP=%d FN=%d"
                     % (r["acc"], modes.get(strip(r["acc"]), "?"), ff*100, r["tp"], r["fp"], r["fn"]))
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=120)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--set", choices=["dev", "heldout"], default="dev")
    ap.add_argument("--final", action="store_true", help="required to touch the held-out set")
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--acc-file", default=None,
                    help="evaluate exactly the accessions listed in this file "
                         "(one per line; bypasses sampling). Used for held-out "
                         "chunk runs, e.g. splits/heldout_chunks/chunk_00.txt")
    ap.add_argument("--raw", default=None,
                    help="override the rawdata directory holding <acc>.fasta and "
                         "<acc>.gb (e.g. a held-out v2 eval folder)")
    a = ap.parse_args()

    if a.raw:
        globals()["RAW"] = a.raw
        print("RAW overridden:", RAW)

    if a.set == "heldout" and not a.final:
        raise SystemExit("REFUSED: held-out is frozen. Pass --final only for the final evaluation.")

    modes = load_modes()
    bm = _bm()
    if a.acc_file:
        accs = [l.strip() for l in open(a.acc_file) if l.strip()]
        print("Acc-file: %d genome tu %s" % (len(accs), a.acc_file))
    else:
        accs = sample(a.set, a.n, a.seed, modes)
    print("Mau: %d genome tu '%s' (seed=%d)" % (len(accs), a.set, a.seed))
    print("Phan bo mau theo mode:", dict(Counter(modes.get(strip(x), "UNKNOWN") for x in accs)))

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    workdir = a.workdir or os.path.join(REPO, "bench_runs", stamp)
    os.makedirs(workdir, exist_ok=True)
    print("Workdir:", workdir, "| workers:", a.workers)

    results = []
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(run_one, acc, workdir, bm): acc for acc in accs}
        done = 0
        for fut in as_completed(futs):
            results.append(fut.result())
            done += 1
            if done % 10 == 0 or done == len(accs):
                print("  ... %d/%d xong" % (done, len(accs)))

    report = aggregate(results, modes)
    print("\n" + "=" * 64)
    print("KET QUA BENCHMARK DA-GENOME (set=%s, n=%d)" % (a.set, len(accs)))
    print("=" * 64)
    print(report)

    json.dump({"args": vars(a), "results": results},
              open(os.path.join(workdir, "results.json"), "w"), indent=2)
    print("\nChi tiet: %s/results.json" % workdir)


if __name__ == "__main__":
    main()
