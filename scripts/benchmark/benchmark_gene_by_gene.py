#!/usr/bin/env python3
"""
Gene-by-gene benchmark for Plastanno output vs a reference GenBank annotation.

Compares CDS / tRNA / rRNA features one-to-one and multiplicity-aware (so the two
IR copies of a gene count as two units). A predicted feature is a TRUE POSITIVE
when a reference feature of the same normalised gene name lies within +/-`tol` bp
at BOTH ends and the two spliced nucleotide sequences are at least `sim` similar
(difflib ratio). Unmatched predictions are false positives (FP); unmatched
references are false negatives (FN).

Also reports a looser DETECTION count (same name + any overlap, ignoring the
+/-tol window) and, as a diagnostic, FP/FN pairs that overlap under different
names (usually synonym issues, not real errors).

Usage:
    python benchmark_gene_by_gene.py reference.gb predicted.gb [--tol 60] [--sim 0.6]
"""
import argparse
import difflib
from collections import defaultdict
from Bio import SeqIO

# canonical synonyms (names are lower-cased + '_'->'-' before lookup)
GENE_SYNONYMS = {
    "psbn": "pbf1",
    "clpp1": "clpp",
    # rRNA naming variants -> canonical (16S, 23S, 4.5S, 5S). Names are already
    # lower-cased and '_'->'-' before lookup. Applied to BOTH ref and pred, so
    # this only normalises the reference's many NCBI spellings to the canonical
    # names the pipeline emits (rrn16/rrn23/rrn4/rrn5).
    "rrn16s": "rrn16", "rrs": "rrn16", "rns": "rrn16",
    "16s ribosomal rna": "rrn16", "16s rrna": "rrn16", "16s": "rrn16",
    "rrn23s": "rrn23", "rrl": "rrn23", "rnl": "rrn23",
    "23s ribosomal rna": "rrn23", "23s rrna": "rrn23", "23s": "rrn23",
    "rrn4.5s": "rrn4", "rrn4.5": "rrn4",
    "4.5s ribosomal rna": "rrn4", "4.5s rrna": "rrn4", "4.5s": "rrn4",
    "rrn5s": "rrn5", "rrf": "rrn5",
    "5s ribosomal rna": "rrn5", "5s rrna": "rrn5", "5s": "rrn5",
}
WANTED = {"CDS", "tRNA", "rRNA"}


def norm_name(name):
    n = (name or "").strip().lower().replace("_", "-")
    return GENE_SYNONYMS.get(n, n)


_AA3 = {"ala":"a","arg":"r","asn":"n","asp":"d","cys":"c","gln":"q","glu":"e",
        "gly":"g","his":"h","ile":"i","leu":"l","lys":"k","met":"m","phe":"f",
        "pro":"p","ser":"s","thr":"t","trp":"w","tyr":"y","val":"v","sec":"u",
        "ifmet":"m","fmet":"m","imet":"m"}

def _parse_trna(name):
    """(aa1, anticodon_RNA) for a tRNA name, else (None, ''). T->U normalised."""
    import re
    n = (name or "").lower().replace("_", "-").replace("(", "-").replace(")", "").replace(" ", "")
    n = re.sub(r"-+", "-", n).strip("-")
    m = re.match(r"^trna?-([a-z]{3})-?([acgut]{3})?$", n)   # tRNA-Leu (CAA) style
    if m and m.group(1) in _AA3:
        return _AA3[m.group(1)], (m.group(2) or "").replace("t", "u")
    m = re.match(r"^trn([a-z]{1,4})-?([acgut]{3})?$", n)    # trnL-CAA / trnL / trnfM-CAU
    if m:
        aa = m.group(1)
        aa1 = "m" if aa in ("fm","fmet","ifmet","imet") else (aa if len(aa)==1 else _AA3.get(aa))
        if aa1:
            return aa1, (m.group(2) or "").replace("t", "u")
    return None, ""

def name_match(a, b):
    """Names equal, OR (for tRNA) same amino acid and equal-or-omitted anticodon.
    Different anticodon or different amino acid never match (real mis-calls kept)."""
    if a == b:
        return True
    a1, c1 = _parse_trna(a); a2, c2 = _parse_trna(b)
    if a1 is None or a2 is None:
        return False
    if a1 != a2:
        return False
    if c1 and c2:
        return c1 == c2
    return True


def load_units(gb_path):
    rec = SeqIO.read(gb_path, "genbank")
    units = []
    for feat in rec.features:
        if feat.type not in WANTED:
            continue
        name = (feat.qualifiers.get("gene")
                or feat.qualifiers.get("product") or [""])[0]
        try:
            seq = str(feat.extract(rec.seq)).upper()
        except Exception:
            seq = ""
        units.append({
            "name": norm_name(name), "raw": name, "type": feat.type,
            "start": int(feat.location.start), "end": int(feat.location.end),
            "strand": feat.location.strand, "seq": seq,
        })
    return units


def coord_close(a, b, tol):
    return abs(a["start"] - b["start"]) <= tol and abs(a["end"] - b["end"]) <= tol


def overlaps(a, b):
    return min(a["end"], b["end"]) > max(a["start"], b["start"])


def similarity(a, b):
    if not a["seq"] or not b["seq"]:
        return 0.0
    return difflib.SequenceMatcher(None, a["seq"], b["seq"], autojunk=False).ratio()


def run(ref, pred, tol, sim):
    ref_used = [False] * len(ref)
    pred_used = [False] * len(pred)
    tp = []
    for pi, p in enumerate(pred):
        best = None
        for ri, r in enumerate(ref):
            if ref_used[ri] or not name_match(r["name"], p["name"]) or not coord_close(p, r, tol):
                continue
            s = similarity(p, r)
            if s >= sim and (best is None or s > best[1]):
                best = (ri, s)
        if best:
            ref_used[best[0]] = True
            pred_used[pi] = True
            tp.append((p, ref[best[0]], best[1]))

    fp = [pred[i] for i in range(len(pred)) if not pred_used[i]]
    fn = [ref[i] for i in range(len(ref)) if not ref_used[i]]

    # looser detection: same name + any overlap (ignore tol), independent matching
    ru = [False] * len(ref)
    detect = 0
    for p in pred:
        for ri, r in enumerate(ref):
            if not ru[ri] and name_match(r["name"], p["name"]) and overlaps(p, r):
                ru[ri] = True
                detect += 1
                break

    # diagnostic: FP/FN that overlap under different names
    mism = []
    used_fn = set()
    for p in fp:
        for j, r in enumerate(fn):
            if j in used_fn:
                continue
            if overlaps(p, r) and p["name"] != r["name"]:
                mism.append((p, r))
                used_fn.add(j)
                break
    return tp, fp, fn, detect, mism


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("predicted")
    ap.add_argument("--tol", type=int, default=60)
    ap.add_argument("--sim", type=float, default=0.6)
    a = ap.parse_args()

    ref = load_units(a.reference)
    pred = load_units(a.predicted)
    tp, fp, fn, detect, mism = run(ref, pred, a.tol, a.sim)

    nTP, nFP, nFN = len(tp), len(fp), len(fn)
    sens = nTP / (nTP + nFN) if (nTP + nFN) else 0.0
    prec = nTP / (nTP + nFP) if (nTP + nFP) else 0.0
    f1 = 2 * sens * prec / (sens + prec) if (sens + prec) else 0.0

    print("=" * 62)
    print("Gene-by-gene benchmark (tol=+/-%dbp, sim>=%.2f)" % (a.tol, a.sim))
    print("=" * 62)
    print("Reference units : %d" % len(ref))
    print("Predicted units : %d" % len(pred))
    print("-" * 62)
    print("TP = %d   FP = %d   FN = %d" % (nTP, nFP, nFN))
    print("Sensitivity = %.1f%%   Precision = %.1f%%   F1 = %.1f%%"
          % (sens * 100, prec * 100, f1 * 100))
    print("Detection (name+overlap, any boundary): %d / %d ref" % (detect, len(ref)))
    print("-" * 62)

    bt = defaultdict(lambda: [0, 0, 0])
    for p, r, s in tp:
        bt[r["type"]][0] += 1
    for r in ref:
        bt[r["type"]][1] += 1
    for p in pred:
        bt[p["type"]][2] += 1
    print("Per type  (TP / ref / pred):")
    for t in ("CDS", "tRNA", "rRNA"):
        print("  %-5s %3d / %3d / %3d" % (t, bt[t][0], bt[t][1], bt[t][2]))
    print("-" * 62)

    if fp:
        print("FALSE POSITIVES (%d):" % nFP)
        for p in sorted(fp, key=lambda x: x["start"]):
            print("  %-12s %-5s %8d-%-8d" % (p["raw"], p["type"], p["start"], p["end"]))
    if fn:
        print("FALSE NEGATIVES (%d):" % nFN)
        for r in sorted(fn, key=lambda x: x["start"]):
            print("  %-12s %-5s %8d-%-8d" % (r["raw"], r["type"], r["start"], r["end"]))
    if mism:
        print("DIAGNOSTIC — overlap under different name (%d, likely synonyms):" % len(mism))
        for p, r in mism:
            print("  pred %-12s vs ref %-12s @ %d" % (p["raw"], r["raw"], p["start"]))


if __name__ == "__main__":
    main()
