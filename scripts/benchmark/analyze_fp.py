#!/usr/bin/env python3
"""
Diagnose Plastanno false positives by CAUSE, to prioritise fixes.

Re-runs the gene-by-gene matching, then classifies each FP into one of:
  boundary  : a same-name reference gene OVERLAPS it (multi-exon/intron coord or
              sequence mismatch -> appears as both FP and FN). Needs intron/boundary work.
  extra_copy: the gene WAS matched correctly elsewhere (a TP with the same name
              exists); this is a surplus duplicate. Needs copy-number reconciliation.
  short_HMM : a short CDS (<= MIN_CDS_BP) with no ref support; classic HMM noise.
  other     : misplaced / cross-hit, none of the above.

If the .gb stores per-feature provenance in /note (engine=..,  C=..), FP are also
cross-tabulated by engine and by confidence band.

Usage: python analyze_fp.py reference.gb predicted.gb [--tol 60] [--sim 0.6]
"""
import argparse, difflib, re
from collections import Counter, defaultdict
from Bio import SeqIO

GENE_SYNONYMS = {"psbn": "pbf1", "clpp1": "clpp"}
WANTED = {"CDS", "tRNA", "rRNA"}
MIN_CDS_BP = 150


def norm_name(n):
    n = (n or "").strip().lower().replace("_", "-")
    return GENE_SYNONYMS.get(n, n)


def load_units(path):
    rec = SeqIO.read(path, "genbank")
    out = []
    for f in rec.features:
        if f.type not in WANTED:
            continue
        name = (f.qualifiers.get("gene") or f.qualifiers.get("product") or [""])[0]
        try:
            seq = str(f.extract(rec.seq)).upper()
        except Exception:
            seq = ""
        note = " ".join(f.qualifiers.get("note", []))
        em = re.search(r"engine\s*=\s*(\S+)", note)
        cm = re.search(r"\bC\s*=\s*([0-9.]+)", note)
        out.append({"name": norm_name(name), "raw": name, "type": f.type,
                    "start": int(f.location.start), "end": int(f.location.end),
                    "seq": seq, "engine": em[1] if em else None,
                    "conf": float(cm[1]) if cm else None})
    return out


def close(a, b, tol):
    return abs(a["start"]-b["start"]) <= tol and abs(a["end"]-b["end"]) <= tol


def overlaps(a, b):
    return min(a["end"], b["end"]) > max(a["start"], b["start"])


def sim(a, b):
    if not a["seq"] or not b["seq"]:
        return 0.0
    return difflib.SequenceMatcher(None, a["seq"], b["seq"], autojunk=False).ratio()


def match(ref, pred, tol, s_thr):
    ru = [False]*len(ref); pu = [False]*len(pred); tp = []
    for pi, p in enumerate(pred):
        best = None
        for ri, r in enumerate(ref):
            if ru[ri] or r["name"] != p["name"] or not close(p, r, tol):
                continue
            s = sim(p, r)
            if s >= s_thr and (best is None or s > best[1]):
                best = (ri, s)
        if best:
            ru[best[0]] = True; pu[pi] = True; tp.append((p, ref[best[0]]))
    fp = [pred[i] for i in range(len(pred)) if not pu[i]]
    fn = [ref[i] for i in range(len(ref)) if not ru[i]]
    return tp, fp, fn


def classify(fp, fn, tp_names):
    cats = {}
    for p in fp:
        if any(r["name"] == p["name"] and overlaps(p, r) for r in fn):
            cats[id(p)] = "boundary"
        elif p["name"] in tp_names:
            cats[id(p)] = "extra_copy"
        elif p["type"] == "CDS" and (p["end"]-p["start"]) <= MIN_CDS_BP:
            cats[id(p)] = "short_HMM"
        else:
            cats[id(p)] = "other"
    return cats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference"); ap.add_argument("predicted")
    ap.add_argument("--tol", type=int, default=60); ap.add_argument("--sim", type=float, default=0.6)
    a = ap.parse_args()
    ref = load_units(a.reference); pred = load_units(a.predicted)
    tp, fp, fn = match(ref, pred, a.tol, a.sim)
    tp_names = set(p["name"] for p, r in tp)
    cats = classify(fp, fn, tp_names)

    print("TP=%d FP=%d FN=%d   (pred=%d ref=%d)" % (len(tp), len(fp), len(fn), len(pred), len(ref)))
    print("-"*54)
    print("FP by CAUSE:")
    cc = Counter(cats.values())
    for k in ("short_HMM", "extra_copy", "boundary", "other"):
        print("  %-11s %d" % (k, cc.get(k, 0)))
    print("-"*54)

    if any(p["engine"] for p in pred):
        print("FP by ENGINE:", dict(Counter(p["engine"] for p in fp)))
        band = lambda c: ("HIGH" if c is not None and c >= 0.8 else
                          "MED" if c is not None and c >= 0.5 else
                          "REVIEW" if c is not None else "??")
        print("FP by CONFIDENCE:", dict(Counter(band(p["conf"]) for p in fp)))
        ex = defaultdict(Counter)
        for p in fp:
            ex[cats[id(p)]][p["engine"]] += 1
        print("cause x engine:")
        for k in ("short_HMM", "extra_copy", "boundary", "other"):
            print("  %-11s %s" % (k, dict(ex[k])))
    else:
        print("(no engine/confidence in .gb /note -> showing cause only)")

    print("-"*54)
    print("'other' FP detail:")
    for p in sorted([x for x in fp if cats[id(x)] == "other"], key=lambda z: z["start"]):
        print("  %-12s %-5s %8d-%-8d eng=%s C=%s" %
              (p["raw"], p["type"], p["start"], p["end"], p["engine"], p["conf"]))


if __name__ == "__main__":
    main()
