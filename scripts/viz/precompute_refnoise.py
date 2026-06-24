#!/usr/bin/env python3
"""Precompute reference-noise statistics over the held-out reference GenBank set:
(1) how many reference plastomes annotate the inverted repeat at all, and
(2) how inconsistently genes are named across references.
Writes docs/figures/fig5_data.json."""
import json, os, glob, warnings
import re
from collections import Counter, defaultdict
from multiprocessing import Pool
warnings.filterwarnings("ignore")
from Bio import SeqIO

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Full RefSeq reference collection used to build the databases (the same 12,581
# plastomes the Background quantifies), so the figure matches the text claim.
RAW = os.environ.get("PLASTANNO_DATA","benchmark_data")+"/rawdata"
accs = [os.path.basename(p)[:-3] for p in glob.glob(f"{RAW}/*.gb")]

def scan(acc):
    p = f"{RAW}/{acc}.gb"
    if not os.path.exists(p):
        return None
    try:
        rec = SeqIO.read(p, "genbank")
    except Exception:
        return None
    # Clean, reproducible definition: the reference annotates the inverted repeat
    # as a standard repeat_region feature >= 10 kb (the plastid IR is ~20-30 kb).
    # Substring matching of free-text notes is deliberately avoided — it both
    # over- and under-counts (e.g. a 1-bp "junction of IRa and LSC" misc_feature).
    has_ir = any(f.type == "repeat_region" and
                 (int(f.location.end) - int(f.location.start)) >= 10000
                 for f in rec.features)
    names = {"rRNA": [], "tRNA": [], "CDS": []}
    for f in rec.features:
        if f.type in names:
            g = (f.qualifiers.get("gene", [""])[0] or f.qualifiers.get("product", [""])[0]).strip()
            if g:
                names[f.type].append(g.lower())
    return has_ir, names

if __name__ == "__main__":
    with Pool(24) as P:
        res = [r for r in P.map(scan, accs) if r]
    n = len(res)
    n_ir = sum(1 for h, _ in res if h)
    # naming heterogeneity: distinct spellings seen for rRNA / tRNA / CDS across refs
    variants = {}
    for t in ("rRNA", "tRNA", "CDS"):
        c = Counter()
        for _h, nm in res:
            c.update(set(nm[t]))
        variants[t] = c
    # canonical rRNA: collapse to rrnXX core and count spelling variants per core
    rrna_core = defaultdict(set)
    for name in variants["rRNA"]:
        m = re.search(r"(\d+\.?\d*)\s*s", name) or re.search(r"rrn\s*(\d+\.?\d*)", name)
        core = m.group(1) if m else "other"
        rrna_core[core].add(name)
    out = {
        "n_refs": n,
        "ir_annotated": n_ir,
        "ir_absent": n - n_ir,
        "ir_absent_pct": round(100*(n-n_ir)/n, 1),
        "rrna_total_spellings": len(variants["rRNA"]),
        "trna_total_spellings": len(variants["tRNA"]),
        "cds_total_spellings": len(variants["CDS"]),
        "rrna_top_spellings": variants["rRNA"].most_common(12),
        "rrna_core_variant_counts": {k: len(v) for k, v in sorted(rrna_core.items())},
    }
    json.dump(out, open(f"{REPO}/docs/figures/fig5_data.json", "w"), indent=2)
    print(f"n_refs={n} | IR absent={out['ir_absent']} ({out['ir_absent_pct']}%)")
    print(f"distinct spellings: rRNA={out['rrna_total_spellings']} tRNA={out['trna_total_spellings']} CDS={out['cds_total_spellings']}")
    print("rRNA core variant counts:", out["rrna_core_variant_counts"])
    print("top rRNA spellings:", out["rrna_top_spellings"][:8])
