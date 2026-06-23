#!/usr/bin/env python3
"""
Diagnose WHY features of a given type fail to match, on the broad bench sample.
Same classification as the rRNA diagnostic, but for any type, and it also
collects the (ref_name, pred_name) pairs at the SAME locus (best-overlapping
pred) -- these are safe synonym candidates because position already agrees.
Usage: python diagnose_type.py <bench_run_dir> <TYPE>     (TYPE: tRNA|rRNA|CDS)
"""
import sys, os, json, re, importlib.util, difflib
from collections import Counter, defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW  = os.environ.get("PLASTANNO_DATA","/data06/users/vutrinh/Apiales_Plastomes_20260514")+"/rawdata"

def _bm():
    p = os.path.join(REPO, "scripts/benchmark/benchmark_gene_by_gene.py")
    spec = importlib.util.spec_from_file_location("bm", p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

def uf(u, *names, default=None):
    for nm in names:
        if isinstance(u, dict) and nm in u: return u[nm]
        if hasattr(u, nm): return getattr(u, nm)
    return default

def main(run_dir, TYPE):
    bm = _bm()
    norm = bm.norm_name
    def sim(a, b):
        if not a or not b: return 0.0
        return difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()
    def overlap(p, r):
        return max(0, min(uf(p,"end"), uf(r,"end")) - max(uf(p,"start"), uf(r,"start")))

    res = json.load(open(os.path.join(run_dir, "results.json")))
    accs = [r["acc"] for r in res["results"] if r.get("status") == "OK"]

    modes = Counter()
    mismatch_pairs = Counter()
    ref_names, pred_names = Counter(), Counter()
    examples = defaultdict(list)
    for acc in accs:
        ref_gb  = os.path.join(RAW, acc + ".gb")
        pred_gb = os.path.join(run_dir, acc, acc + ".gb")
        if not (os.path.exists(ref_gb) and os.path.exists(pred_gb)): continue
        try:
            ref  = bm.load_units(ref_gb)
            pred = bm.load_units(pred_gb)
        except Exception:
            continue
        ref_t  = [u for u in ref  if uf(u,"type")==TYPE]
        pred_t = [u for u in pred if uf(u,"type")==TYPE]
        for u in ref_t:  ref_names[uf(u,"name")]  += 1
        for u in pred_t: pred_names[uf(u,"name")] += 1

        for rr in ref_t:
            cands = [pr for pr in pred_t if overlap(pr, rr) > 0]
            if not cands:
                modes["NO_PRED_NEARBY"] += 1
                if len(examples["NO_PRED_NEARBY"])<10:
                    examples["NO_PRED_NEARBY"].append("%s %s @%d-%d" % (acc, uf(rr,"name"), uf(rr,"start"), uf(rr,"end")))
                continue
            pr = max(cands, key=lambda p: overlap(p, rr))
            name_ok = uf(pr,"name") == uf(rr,"name")
            ds, de = uf(pr,"start")-uf(rr,"start"), uf(pr,"end")-uf(rr,"end")
            pos_ok = abs(ds)<=60 and abs(de)<=60
            s = sim(uf(pr,"seq") or "", uf(rr,"seq") or "")
            seq_ok = s >= 0.6
            if name_ok and pos_ok and seq_ok:
                modes["MATCHED"] += 1
            elif not name_ok:
                modes["NAME_MISMATCH"] += 1
                mismatch_pairs[(uf(rr,"name"), uf(pr,"name"))] += 1
            elif not pos_ok:
                modes["BOUNDARY_OFF"] += 1
                if len(examples["BOUNDARY_OFF"])<10:
                    examples["BOUNDARY_OFF"].append("%s %s dstart=%+d dend=%+d (reflen=%d)" % (acc, uf(rr,"name"), ds, de, uf(rr,"end")-uf(rr,"start")))
            else:
                modes["SEQ_LOW"] += 1
                if len(examples["SEQ_LOW"])<10:
                    examples["SEQ_LOW"].append("%s %s sim=%.2f reflen=%d predlen=%d strand_r=%s strand_p=%s" % (acc, uf(rr,"name"), s, len(uf(rr,"seq") or ""), len(uf(pr,"seq") or ""), uf(rr,"strand"), uf(pr,"strand")))

    tot = sum(modes.values())
    print("=== %s chan doan (%d genome, %d %s tham chieu) ===" % (TYPE, len(accs), tot, TYPE))
    for m in ("MATCHED","NAME_MISMATCH","BOUNDARY_OFF","SEQ_LOW","NO_PRED_NEARBY"):
        print("  %-16s %5d  (%.1f%%)" % (m, modes[m], 100*modes[m]/max(1,tot)))

    print("\nTop 25 cap (ref_name -> pred_name) CUNG LOCUS nhung khac ten (ung vien synonym AN TOAN):")
    for (rn, pn), c in mismatch_pairs.most_common(25):
        print("    %-18s -> %-18s  x%d" % (rn, pn, c))

    print("\nTop ten %s THAM CHIEU (25):" % TYPE)
    for nm, c in ref_names.most_common(25): print("    %-18s %d" % (nm, c))
    print("Top ten %s DU DOAN (25):" % TYPE)
    for nm, c in pred_names.most_common(25): print("    %-18s %d" % (nm, c))

    for m in ("NO_PRED_NEARBY","BOUNDARY_OFF","SEQ_LOW"):
        if examples[m]:
            print("\n  vi du %s:" % m)
            for e in examples[m][:10]: print("    " + e)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "tRNA")
