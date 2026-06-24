"""Reference-transfer splice refinement for cis multi-exon CDS.

Reference exon sequences end exactly at their splice junctions, so BLASTN of the
gene's genomic region (coding orientation) against the per-gene reference exon
panel locates every exon plus the start/stop. The outer ends (start codon, stop)
are taken from the alignment; each internal junction is then resolved by a small
constrained search for canonical GT (donor) / AG (acceptor) positions that keep
the spliced CDS in-frame and free of internal stops. Accepted only if the result
translates cleanly. Genes without a panel and trans-spliced genes (rps12) fall
back to merging artefact short introns.
"""
import os, json, subprocess, tempfile, itertools

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def _rc(s): return s.translate(_COMP)[::-1]

MIN_INTRON = 100
PAD        = 75
JW         = 15      # junction search window (bp)
from ..paths import config_dir
_HERE = os.path.dirname(__file__)
_BDB  = str(config_dir() / "boundary_db" / "exons")
_META = None

_TPL = None
def _meta():
    global _META
    if _META is None:
        try: _META = json.load(open(str(config_dir() / "boundary_db" / "meta.json")))
        except Exception: _META = {}
    return _META

def _tpl(gene):
    global _TPL
    if _TPL is None:
        try: _TPL = json.load(open(str(config_dir() / "exon_templates.json")))
        except Exception: _TPL = {}
    return (_TPL.get(gene) or {}).get("exon_lens")

def _dist(exons, L):
    return sum(abs(a - b) for a, b in zip(sorted(e - s for s, e in exons), sorted(L)))

def _prot(seq):
    seq = seq[: len(seq)//3*3]
    try:
        from Bio.Seq import Seq
        return str(Seq(seq).translate(table=11))
    except Exception:
        return ""

def _merge_short(ex):
    m = [list(ex[0])]
    for s, e in ex[1:]:
        if 0 <= s - m[-1][1] < MIN_INTRON: m[-1][1] = e
        else: m.append([s, e])
    return [(s, e) for s, e in m]

def _blast_exons(region, gene):
    if not os.path.exists(_BDB + ".nin"):
        return {}
    qf = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    qf.write(">q\n" + region + "\n"); qf.close()
    try:
        r = subprocess.run(
            ["blastn", "-task", "blastn-short", "-query", qf.name, "-db", _BDB,
             "-outfmt", "6 sseqid qstart qend sstart send slen bitscore", "-strand", "plus",
             "-word_size", "7", "-dust", "no", "-evalue", "1", "-max_target_seqs", "300"],
            capture_output=True, text=True, timeout=60)
    except Exception:
        os.unlink(qf.name); return {}
    os.unlink(qf.name)
    n = len(region); best = {}
    for line in r.stdout.splitlines():
        sid, qs, qe, ss, se, slen, bits = line.split("\t")
        g, ei, _ = sid.split("|")
        if g != gene: continue
        i = int(ei[1:]); qs, qe, ss, se, slen, bits = int(qs), int(qe), int(ss), int(se), int(slen), float(bits)
        cs = max(0, (qs - 1) - (ss - 1)); ce = min(n, qe + (slen - se))
        if i not in best or bits > best[i][2]:
            best[i] = (cs, ce, bits)
    return {i: (v[0], v[1]) for i, v in best.items()}

def _cands(C, pos, dinuc, at_end):
    # canonical GT/AG positions in the window (preferred, nearest first) ...
    n = len(C); gt = []
    for off in range(0, JW + 1):
        for p in ((pos,) if off == 0 else (pos + off, pos - off)):
            if 2 <= p <= n - 2 and (C[p-2:p] if at_end else C[p:p+2]) == dinuc:
                gt.append(p)
    # ... plus phase neighbours of the alignment boundary as a fallback, because
    # plastid group-II introns do not always end in a strict GT..AG.
    neigh = [pos + d for d in (0, -1, 1, -2, 2, -3, 3) if 2 <= pos + d <= n - 2]
    seen, out = set(), []
    for p in gt + neigh:
        if p not in seen:
            seen.add(p); out.append(p)
    return out[:8]

def _constrained(C, cex):
    """Resolve internal junctions to in-frame GT..AG; outer ends fixed. -> coding exons or None."""
    k = len(cex)
    Ds = [_cands(C, cex[j][1], "GT", False) for j in range(k - 1)]
    As = [_cands(C, cex[j+1][0], "AG", True) for j in range(k - 1)]
    if any(not d for d in Ds) or any(not a for a in As):
        return None
    best = None
    for cd in itertools.product(*Ds):
        for ca in itertools.product(*As):
            exons = []; ps = cex[0][0]; ok = True
            for j in range(k - 1):
                d, a = cd[j], ca[j]
                if d <= ps + 2 or a <= d + MIN_INTRON: ok = False; break
                exons.append((ps, d)); ps = a
            if not ok: continue
            exons.append((ps, cex[-1][1]))
            if any(e <= s for s, e in exons): continue
            if sum(e - s for s, e in exons) % 3: continue
            prot = _prot("".join(C[s:e] for s, e in exons))
            if not prot or prot[:-1].count("*"): continue
            shift = sum(abs(cd[j]-cex[j][1]) + abs(ca[j]-cex[j+1][0]) for j in range(k-1))
            if best is None or shift < best[0]:
                best = (shift, exons)
    return best[1] if best else None

def _refine_one(feat, g):
    ex = sorted(feat.exons)
    if len(ex) < 2: return False
    ex = _merge_short(ex)
    changed = (ex != sorted(feat.exons))

    M = _meta().get(feat.gene_name)
    if M and M.get("n_exons") == len(ex):
        strand = feat.strand
        r0 = max(0, ex[0][0] - PAD); r1 = min(len(g), ex[-1][1] + PAD)
        C = _rc(g[r0:r1]) if strand == -1 else g[r0:r1]
        hsp = _blast_exons(C, feat.gene_name)
        if len(hsp) == len(ex):
            cex = sorted(hsp[i] for i in sorted(hsp))
            ce = _constrained(C, cex)
            if ce:
                if strand == -1:
                    gex = sorted((r0 + (len(C) - e), r0 + (len(C) - s)) for s, e in ce)
                else:
                    gex = sorted((r0 + s, r0 + e) for s, e in ce)
                L = _tpl(feat.gene_name)
                closer = (L is None) or (_dist(gex, L) < _dist(ex, L))
                if gex != ex and closer and all(e > s for s, e in gex) and \
                   all(gex[i+1][0] > gex[i][1] for i in range(len(gex)-1)):
                    feat.exons = gex; feat.start = gex[0][0]; feat.end = gex[-1][1]
                    feat.has_intron = len(gex) > 1
                    return True
    if changed:
        feat.exons = ex; feat.start = ex[0][0]; feat.end = ex[-1][1]; feat.has_intron = len(ex) > 1
    return changed

def refine_all(annotations, genome_seq):
    n = 0
    for f in annotations:
        if getattr(f, "gene_type", None) != "CDS": continue
        if getattr(f, "exon_strands", None): continue
        if not getattr(f, "exons", None) or len(f.exons) < 2: continue
        try:
            if _refine_one(f, genome_seq): n += 1
        except Exception: pass
    return n
