#!/usr/bin/env python3
"""
Calibrate a minimum-length quality threshold for CDS, on the BROAD DEV sample.
Reanalyses kept predictions of a finished bench run; for each predicted CDS,
computes spliced_len / expected_len and whether it is TP or FP.
Usage: python calibrate_threshold.py <bench_run_dir>
"""
import sys, os, json, re, importlib.util
from collections import Counter

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW  = os.environ.get("PLASTANNO_DATA","benchmark_data")+"/rawdata"
CAT  = os.path.join(REPO, "database", "gene_catalog.json")

def _bm():
    p = os.path.join(REPO, "scripts/benchmark/benchmark_gene_by_gene.py")
    spec = importlib.util.spec_from_file_location("bm", p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

def ufield(u, *names, default=None):
    for nm in names:
        if isinstance(u, dict) and nm in u:
            return u[nm]
        if hasattr(u, nm):
            return getattr(u, nm)
    return default

def unit_len(u):
    seq = ufield(u, "seq", "sequence")
    if seq is not None:
        return len(seq)
    s, e = ufield(u, "start"), ufield(u, "end")
    if s is not None and e is not None:
        return abs(e - s)
    return None

def unit_type(u): return ufield(u, "type", "gene_type", default="?")
def unit_name(u): return ufield(u, "raw", "name", "gene_name", default="?")

def build_ci_index(cat):
    idx = {}
    for k in cat:
        idx.setdefault(k.lower(), k)
    return idx

def get_expected(cat, ci_index, raw):
    """original case -> case-insensitive -> resolve synonym_of."""
    for key in (raw, ci_index.get(raw.lower())):
        if not key:
            continue
        e = cat.get(key, {})
        if e.get("expected_len"):
            return e["expected_len"]
        syn = e.get("synonym_of")
        if syn:
            for sk in (syn, ci_index.get(syn.lower())):
                if sk and cat.get(sk, {}).get("expected_len"):
                    return cat[sk]["expected_len"]
    return 0

def main(run_dir):
    bm = _bm()
    cat = json.load(open(CAT))
    ci_index = build_ci_index(cat)
    res = json.load(open(os.path.join(run_dir, "results.json")))
    accs = [r["acc"] for r in res["results"] if r.get("status") == "OK"]
    print("Phan tich %d genome tu %s" % (len(accs), run_dir))
    n_exp_genes = sum(1 for k, v in cat.items() if isinstance(v, dict) and v.get("expected_len"))
    print("Gene trong catalog co expected_len: %d" % n_exp_genes)

    rows = []
    no_exp = 0
    missing = Counter()
    for acc in accs:
        ref_gb  = os.path.join(RAW, acc + ".gb")
        pred_gb = os.path.join(run_dir, acc, acc + ".gb")
        if not (os.path.exists(ref_gb) and os.path.exists(pred_gb)):
            continue
        try:
            ref  = bm.load_units(ref_gb)
            pred = bm.load_units(pred_gb)
            tp, fp, fn, detect, mism = bm.run(ref, pred, 60, 0.6)
        except Exception:
            continue
        tp_ids = set(id(p) for p, r, s in tp)
        for u in pred:
            if unit_type(u) != "CDS":
                continue
            raw = re.sub(r"\.\d+$", "", str(unit_name(u)))
            exp = get_expected(cat, ci_index, raw)
            if not exp:
                no_exp += 1
                missing[raw] += 1
                continue
            L = unit_len(u)
            if L is None:
                continue
            rows.append((L / exp, id(u) in tp_ids))

    tp_r = [r for r, t in rows if t]
    fp_r = [r for r, t in rows if not t]
    print("\nCDS co expected_len: %d (TP=%d, FP=%d) | CDS thieu: %d" %
          (len(rows), len(tp_r), len(fp_r), no_exp))
    print("Top 15 ten gene THIEU expected_len (loi-case se la gene pho bien; chua-co se la orf/hiem):")
    for nm, c in missing.most_common(15):
        print("    %-16s %d" % (nm, c))

    buckets = [(0,0.1),(0.1,0.2),(0.2,0.3),(0.3,0.4),(0.4,0.5),(0.5,0.6),
               (0.6,0.7),(0.7,0.8),(0.8,0.9),(0.9,1.1),(1.1,1.3),(1.3,9e9)]
    print("\nPhan bo ratio = spliced_len / expected_len:")
    print("  %-12s %8s %8s" % ("bucket", "TP", "FP"))
    for lo, hi in buckets:
        nt = sum(1 for r in tp_r if lo <= r < hi)
        nf = sum(1 for r in fp_r if lo <= r < hi)
        lab = "<%.1f" % hi if lo == 0 else (">=%.1f" % lo if hi > 9e8 else "%.1f-%.1f" % (lo, hi))
        print("  %-12s %8d %8d  %s" % (lab, nt, nf, "#" * min(40, nf)))

    print("\nNeu loai CDS co ratio < theta (chi CDS co expected_len):")
    print("  %-8s %-16s %-16s" % ("theta", "FP loai bo", "TP mat oan"))
    for th in (0.2,0.3,0.4,0.5,0.6,0.7):
        fpd = sum(1 for r in fp_r if r < th)
        tpl = sum(1 for r in tp_r if r < th)
        print("  %-8.1f %-16s %-16s" % (th,
              "%d (%.0f%%)" % (fpd, 100*fpd/max(1,len(fp_r))),
              "%d (%.1f%%)" % (tpl, 100*tpl/max(1,len(tp_r)))))

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
