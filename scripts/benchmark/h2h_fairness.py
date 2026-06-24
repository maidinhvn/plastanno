import os, warnings, importlib.util, collections
warnings.filterwarnings("ignore")
from collections import Counter
from multiprocessing import Pool
spec = importlib.util.spec_from_file_location("bm", "scripts/benchmark/benchmark_gene_by_gene.py")
bm = importlib.util.module_from_spec(spec); spec.loader.exec_module(bm)
EVAL = os.environ.get("PLASTANNO_DATA","benchmark_data")+"/eval_v2_raw"; H = "bench_runs/h2h"
subset = [l.strip() for l in open("splits/h2h_subset.txt") if l.strip()]

def work(a):
    try:
        r = bm.load_units(f"{EVAL}/{a}.gb")
        g = bm.load_units(f"{H}/pga_out/{a}.gb")
        p = bm.load_units(f"{H}/plast_out/{a}/{a}.gb")
    except Exception:
        return None
    def cs(ref, pred):
        tp, fp, fn, _, _ = bm.run(ref, pred, 60, 0.6)
        return (Counter(rr['type'] for _x, rr, _y in tp),
                Counter(u['type'] for u in ref), Counter(u['type'] for u in pred))
    return (cs(r, g), cs(r, p),
            Counter(u['name'] for u in r if u['type'] == 'rRNA'),
            Counter(u['name'] for u in g if u['type'] == 'rRNA'))

if __name__ == "__main__":
    with Pool(24) as P:
        res = [x for x in P.map(work, subset) if x]
    print(f"paired sach={len(res)}")
    for idx, lab in ((0, "PGA"), (1, "Plastanno")):
        TP = Counter(); RF = Counter(); PD = Counter()
        for r in res:
            TP += r[idx][0]; RF += r[idx][1]; PD += r[idx][2]
        print(lab + ":")
        for t in ("CDS", "tRNA", "rRNA"):
            sn = 100*TP[t]/RF[t] if RF[t] else 0
            pr = 100*TP[t]/PD[t] if PD[t] else 0
            print(f"    {t:5s} Sn={sn:5.1f} Pr={pr:5.1f}  (TP={TP[t]} ref={RF[t]} pred={PD[t]})")
    rn = Counter(); pn = Counter()
    for r in res:
        rn += r[2]; pn += r[3]
    print(f"\nrRNA ten REF: {dict(rn)}")
    print(f"rRNA ten PGA: {dict(pn)}")
    tiers = collections.Counter()
    for l in open("splits/h2h_refmap.tsv"):
        pp = l.rstrip().split("\t")
        if len(pp) >= 3:
            tiers[pp[2]] += 1
    print(f"\nrefmap PGA reference tier: {dict(tiers)}")
