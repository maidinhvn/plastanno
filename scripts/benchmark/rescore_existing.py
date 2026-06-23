#!/usr/bin/env python3
"""Re-score a finished bench run from its saved .gb files (no pipeline re-run),
parallel across genomes. Usage: rescore_existing.py <dir> [workers]"""
import sys, os, json, importlib.util
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW  = os.environ.get("PLASTANNO_DATA","/data06/users/vutrinh/Apiales_Plastomes_20260514")+"/rawdata"

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

mg = load(os.path.join(REPO, "scripts/benchmark/multi_genome_bench.py"), "mg")
bm = mg._bm()

def score_one(acc, run_dir):
    ref_gb  = os.path.join(RAW, acc + ".gb")
    pred_gb = os.path.join(run_dir, acc, acc + ".gb")
    if not (os.path.exists(ref_gb) and os.path.exists(pred_gb)):
        return None
    try:
        ref  = bm.load_units(ref_gb)
        pred = bm.load_units(pred_gb)
        tp, fp, fn, detect, mism = bm.run(ref, pred, 60, 0.6)
    except Exception as e:
        return {"acc": acc, "status": "ERR", "err": str(e)[:120]}
    return {
        "acc": acc, "status": "OK",
        "tp": len(tp), "fp": len(fp), "fn": len(fn),
        "ref": len(ref), "pred": len(pred),
        "tp_by_type": dict(Counter(r["type"] for _p, r, _s in tp)),
        "ref_by_type": dict(Counter(u["type"] for u in ref)),
        "pred_by_type": dict(Counter(u["type"] for u in pred)),
    }

def main(run_dir, workers):
    modes = mg.load_modes()
    res = json.load(open(os.path.join(run_dir, "results.json")))
    accs = [r["acc"] for r in res["results"] if r.get("status") == "OK"]
    results = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(score_one, a, run_dir): a for a in accs}
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            if r: results.append(r)
            done += 1
            if done % 25 == 0 or done == len(accs):
                print("  ... %d/%d" % (done, len(accs)))
    ok = [r for r in results if r["status"] == "OK"]
    print("\nRESCORE (tu .gb co san, song song, dung name_match moi):")
    print(mg.aggregate(ok, modes))

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".",
         int(sys.argv[2]) if len(sys.argv) > 2 else 16)
