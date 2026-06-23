#!/usr/bin/env python3
"""
Split tRNA NAME_MISMATCH (same locus, different name) into:
  SAFE_SYNONYM   ref dropped the anticodon / used a known variant of the SAME
                 amino acid+anticodon as pred  -> benchmark under-counts
  WRONG_ANTICODON/WRONG_AA  pred named a genuinely DIFFERENT gene -> real error
Usage: python split_trna_mismatch.py <bench_run_dir>
"""
import sys, os, json, re, importlib.util
from collections import Counter

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

AA3 = {"ala":"a","arg":"r","asn":"n","asp":"d","cys":"c","gln":"q","glu":"e",
       "gly":"g","his":"h","ile":"i","leu":"l","lys":"k","met":"m","phe":"f",
       "pro":"p","ser":"s","thr":"t","trp":"w","tyr":"y","val":"v","sec":"u",
       "ifmet":"m","fmet":"m"}
def parse_trna(name):
    n = (name or "").lower().replace("_","-").replace("(","-").replace(")","").strip("-")
    m = re.match(r"^trn([a-z]{1,3})-?([acgu]{3})?$", n)
    if not m:
        m2 = re.match(r"^trn-?([a-z]+?)-?([acgu]{3})?$", n)
        if m2:
            aa, anti = m2.group(1), m2.group(2)
            return (AA3.get(aa, aa[:1] if aa else ""), anti or "")
        return "?", ""
    aa, anti = m.group(1), m.group(2)
    if aa in ("fm","fmet","ifmet"): aa1 = "m"
    elif len(aa) == 1: aa1 = aa
    else: aa1 = AA3.get(aa, aa[:1])
    return aa1, (anti or "")

def classify(ref_name, pred_name):
    ra, rc = parse_trna(ref_name); pa, pc = parse_trna(pred_name)
    if ra == "?" or pa == "?": return "UNPARSED"
    if ra != pa: return "WRONG_AA"
    if rc and pc and rc == pc: return "SAFE_SYNONYM"
    if (not rc) or (not pc): return "SAFE_SYNONYM"
    return "WRONG_ANTICODON"

def main(run_dir):
    bm = _bm()
    def overlap(p, r):
        return max(0, min(uf(p,"end"), uf(r,"end")) - max(uf(p,"start"), uf(r,"start")))
    res = json.load(open(os.path.join(run_dir, "results.json")))
    accs = [r["acc"] for r in res["results"] if r.get("status") == "OK"]
    cls = Counter(); safe_pairs = Counter(); wrong_pairs = Counter()
    for acc in accs:
        ref_gb  = os.path.join(RAW, acc + ".gb"); pred_gb = os.path.join(run_dir, acc, acc + ".gb")
        if not (os.path.exists(ref_gb) and os.path.exists(pred_gb)): continue
        try:
            ref  = bm.load_units(ref_gb); pred = bm.load_units(pred_gb)
        except Exception:
            continue
        ref_t  = [u for u in ref  if uf(u,"type")=="tRNA"]
        pred_t = [u for u in pred if uf(u,"type")=="tRNA"]
        for rr in ref_t:
            cands = [pr for pr in pred_t if overlap(pr, rr) > 0]
            if not cands: continue
            pr = max(cands, key=lambda p: overlap(p, rr))
            if uf(pr,"name") == uf(rr,"name"): continue
            c = classify(uf(rr,"name"), uf(pr,"name")); cls[c] += 1
            if c == "SAFE_SYNONYM": safe_pairs[(uf(rr,"name"), uf(pr,"name"))] += 1
            elif c in ("WRONG_AA","WRONG_ANTICODON"): wrong_pairs[(uf(rr,"name"), uf(pr,"name"))] += 1
    tot = sum(cls.values())
    print("=== Phan loai tRNA NAME_MISMATCH (%d cap cung locus khac ten) ===" % tot)
    for k in ("SAFE_SYNONYM","WRONG_ANTICODON","WRONG_AA","UNPARSED"):
        print("  %-16s %5d  (%.1f%%)" % (k, cls[k], 100*cls[k]/max(1,tot)))
    print("\nSAFE_SYNONYM — top cap (map an toan):")
    for (rn,pn),c in safe_pairs.most_common(20): print("    %-16s -> %-16s x%d" % (rn,pn,c))
    print("\nWRONG (anticodon/aa) — top cap (LOI PIPELINE THAT, KHONG map):")
    for (rn,pn),c in wrong_pairs.most_common(20): print("    %-16s vs %-16s x%d" % (rn,pn,c))

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
