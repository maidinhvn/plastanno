#!/usr/bin/env python3
"""
Diagnose WHY rRNA features fail to match, on the broad bench sample.
For each reference rRNA, find the best-overlapping predicted rRNA and classify:
  MATCHED        name+position(<=60bp)+seq(difflib>=0.6) all OK
  NAME_MISMATCH  a pred rRNA overlaps but the (normalised) name differs
  BOUNDARY_OFF   name OK but start/end off by > 60 bp
  SEQ_LOW        name+position OK but difflib < 0.6 (likely truncated DB / partial)
  NO_PRED_NEARBY no predicted rRNA overlaps this locus at all
Usage: python diagnose_rrna.py <bench_run_dir>
"""
import sys, os, json, re, importlib.util, difflib
from collections import Counter, defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW  = os.environ.get("PLASTANNO_DATA","benchmark_data")+"/rawdata"

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

def main(run_dir):
    bm = _bm()
    SYN = getattr(bm, "GENE_SYNONYMS", {}) or {}
    def norm(name):
        n = re.sub(r"\.\d+$", "", str(name)).lower()
        return SYN.get(n, n)
    def sim(a, b):
        if not a or not b: return 0.0
        return difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()
    def overlap(p, r):
        return max(0, min(uf(p,"end"), uf(r,"end")) - max(uf(p,"start"), uf(r,"start")))

    res = json.load(open(os.path.join(run_dir, "results.json")))
    accs = [r["acc"] for r in res["results"] if r.get("status") == "OK"]

    modes = Counter()
    ref_names, pred_names = Counter(), Counter()
    examples = defaultdict(list)
    bnd_deltas = []
    for acc in accs:
        ref_gb  = os.path.join(RAW, acc + ".gb")
        pred_gb = os.path.join(run_dir, acc, acc + ".gb")
        if not (os.path.exists(ref_gb) and os.path.exists(pred_gb)): continue
        try:
            ref  = bm.load_units(ref_gb)
            pred = bm.load_units(pred_gb)
        except Exception:
            continue
        ref_r  = [u for u in ref  if uf(u,"type")=="rRNA"]
        pred_r = [u for u in pred if uf(u,"type")=="rRNA"]
        for u in ref_r:  ref_names[norm(uf(u,"raw","name"))]  += 1
        for u in pred_r: pred_names[norm(uf(u,"raw","name"))] += 1

        for rr in ref_r:
            cands = [pr for pr in pred_r if overlap(pr, rr) > 0]
            if not cands:
                modes["NO_PRED_NEARBY"] += 1
                if len(examples["NO_PRED_NEARBY"])<6:
                    examples["NO_PRED_NEARBY"].append("%s %s @%d-%d" % (acc, norm(uf(rr,"raw","name")), uf(rr,"start"), uf(rr,"end")))
                continue
            pr = max(cands, key=lambda p: overlap(p, rr))
            name_ok = norm(uf(pr,"raw","name")) == norm(uf(rr,"raw","name"))
            ds, de = uf(pr,"start")-uf(rr,"start"), uf(pr,"end")-uf(rr,"end")
            pos_ok = abs(ds)<=60 and abs(de)<=60
            s = sim(uf(pr,"seq") or "", uf(rr,"seq") or "")
            seq_ok = s >= 0.6
            if name_ok and pos_ok and seq_ok:
                modes["MATCHED"] += 1
            elif not name_ok:
                modes["NAME_MISMATCH"] += 1
                if len(examples["NAME_MISMATCH"])<8:
                    examples["NAME_MISMATCH"].append("%s ref=%s vs pred=%s @%d" % (acc, norm(uf(rr,"raw","name")), norm(uf(pr,"raw","name")), uf(rr,"start")))
            elif not pos_ok:
                modes["BOUNDARY_OFF"] += 1; bnd_deltas.append((ds,de))
                if len(examples["BOUNDARY_OFF"])<8:
                    examples["BOUNDARY_OFF"].append("%s %s dstart=%+d dend=%+d (reflen=%d)" % (acc, norm(uf(rr,"raw","name")), ds, de, uf(rr,"end")-uf(rr,"start")))
            else:
                modes["SEQ_LOW"] += 1
                if len(examples["SEQ_LOW"])<8:
                    examples["SEQ_LOW"].append("%s %s sim=%.2f reflen=%d predlen=%d" % (acc, norm(uf(rr,"raw","name")), s, len(uf(rr,"seq") or ""), len(uf(pr,"seq") or "")))

    tot = sum(modes.values())
    print("=== rRNA chan doan (%d genome, %d rRNA tham chieu) ===" % (len(accs), tot))
    for m in ("MATCHED","NAME_MISMATCH","BOUNDARY_OFF","SEQ_LOW","NO_PRED_NEARBY"):
        print("  %-16s %5d  (%.1f%%)" % (m, modes[m], 100*modes[m]/max(1,tot)))

    print("\nTen rRNA THAM CHIEU:", dict(ref_names))
    print("Ten rRNA DU DOAN   :", dict(pred_names))
    if bnd_deltas:
        import statistics as st
        starts=[d[0] for d in bnd_deltas]; ends=[d[1] for d in bnd_deltas]
        print("\nBoundary deltas (khi BOUNDARY_OFF): dstart median=%.0f, dend median=%.0f" %
              (st.median(starts), st.median(ends)))
    for m in ("NAME_MISMATCH","BOUNDARY_OFF","SEQ_LOW","NO_PRED_NEARBY"):
        if examples[m]:
            print("\n  vi du %s:" % m)
            for e in examples[m]: print("    " + e)

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
