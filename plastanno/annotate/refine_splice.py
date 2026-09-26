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
from ..core import coords as _coords
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
    """Refine the exon junctions of one multi-exon CDS.

    Everything happens in WINDOW coordinates: a window that follows the gene round
    the origin, in which the gene is contiguous and has exactly the number of exons
    it biologically has. That matters for more than the sequence slice. An earlier
    version worked on the genomic exon list, where an origin cutting THROUGH an
    exon splits it into two arcs — so `n_exons == len(exons)` failed, the
    refinement was skipped, and a wrapped gene took a different path from the same
    gene presented differently. Its exon-length distance to the template was wrong
    for the same reason: two arcs of one exon have the wrong lengths.

    Note the failure was conditional, not universal: a gene whose origin falls in
    an intron, or outside it, keeps its arc count and was refined normally. Only
    the cut-through-an-exon case diverged.
    """
    glen = len(g)
    if not feat.exons:
        return False
    win = _coords.open_window(g, feat.start, feat.end, glen, buffer=PAD)
    wlen = len(win.seq)

    # the exons, in window coordinates, with any arc split by the origin rejoined
    loc = []
    for s0, e0, _st in _coords.transcript_parts(feat, glen):
        ls = (s0 - win.start) % glen
        loc.append((ls, ls + (e0 - s0)))
    loc.sort()
    ex = []
    for iv in loc:
        if ex and iv[0] <= ex[-1][1]:
            ex[-1] = (ex[-1][0], max(ex[-1][1], iv[1]))
        else:
            ex.append(iv)
    if len(ex) < 2:
        return False
    merged = _merge_short(ex)
    changed = (merged != ex)
    ex = merged

    def to_genome(local_exons):
        arcs = []
        for a, b in local_exons:
            arcs.extend(win.arcs(a, b))
        return sorted(arcs)

    def write_back(local_exons):
        arcs = to_genome(local_exons)
        lo = min(a for a, _ in local_exons)
        hi = max(b for _, b in local_exons)
        feat.exons = arcs
        feat.start = win.to_genomic(lo)
        feat.end = win.interval(lo, hi)[1]
        feat.has_intron = len(local_exons) > 1

    M = _meta().get(feat.gene_name)
    if M and M.get("n_exons") == len(ex):
        strand = feat.strand
        C = _rc(win.seq) if strand == -1 else win.seq
        hsp = _blast_exons(C, feat.gene_name)
        if len(hsp) == len(ex):
            cex = sorted(hsp[i] for i in sorted(hsp))
            ce = _constrained(C, cex)
            if ce:
                if strand == -1:
                    gex = sorted((wlen - e, wlen - s) for s, e in ce)
                else:
                    gex = sorted((s, e) for s, e in ce)
                L = _tpl(feat.gene_name)
                closer = (L is None) or (_dist(gex, L) < _dist(ex, L))
                if gex != ex and closer and all(e > s for s, e in gex) and \
                   all(gex[i+1][0] > gex[i][1] for i in range(len(gex)-1)) and \
                   all(0 <= s and e <= wlen for s, e in gex):
                    write_back(gex)
                    return True
    if changed:
        write_back(ex)
    return changed


_STOP = {"TAA", "TAG", "TGA"}
_SNAP_MAX_CODONS = 3          # only look ≤ 9 bp downstream (RNA-editing-safe)

def _snap_terminal_stop(feat, g, gene_catalog):
    """Extend a CDS 3' end by up to _SNAP_MAX_CODONS codons to reach the nearest
    in-frame stop, when the spliced CDS does not already end in a stop.

    Deliberately CONSERVATIVE. Many plastid CDS (esp. multi-exon) are annotated a
    codon short of the terminal stop; snapping to the adjacent genomic stop fixes
    that. The window is tiny (≤ 9 bp) and the result must not overshoot the
    expected length, so genes whose functional stop is created by RNA editing
    (genomic read-through, e.g. some ndh/rpl2) are NOT force-extended to a distant
    downstream stop — if no genomic stop sits right at the terminus we leave it.
    """
    glen = len(g)
    # The coding sequence in TRANSCRIPTION order. Concatenating sorted(exons) puts
    # a wrapped gene's parts the wrong way round — the piece after the origin first
    # — so the frame, and therefore the terminal codon, were read from the wrong
    # place.
    seq = _coords.extract(g, feat, glen)
    if len(seq) < 6 or len(seq) % 3 != 0:      # need an intact frame to snap
        return False
    if seq[-3:].upper() in _STOP:               # already terminated
        return False

    cur_len = _coords.spliced_length(feat, glen)
    exp = (gene_catalog.get(feat.gene_name, {}) or {}).get("expected_len", 0)
    # The exons in transcription order; the LAST one carries the 3' end, wherever
    # that sits on the circle.
    parts = _coords.transcript_parts(feat, glen)
    if not parts:
        return False
    s0, e0, st0 = parts[-1]

    def codon_at(k):
        """The k-th codon downstream of the 3' end, read round the origin."""
        if st0 == 1:
            a = (e0 + 3 * k) % glen
            piece = g[a:a + 3] if a + 3 <= glen else g[a:] + g[:a + 3 - glen]
        else:
            b = (s0 - 3 * k) % glen
            a = (b - 3) % glen
            piece = g[a:b] if a < b else g[a:] + g[:b]
            piece = _rc(piece)
        return piece.upper()

    for k in range(_SNAP_MAX_CODONS):
        codon = codon_at(k)
        if len(codon) != 3:
            break
        if codon not in _STOP:
            continue
        grow = 3 * (k + 1)
        if exp and (cur_len + grow) / exp > 1.2:      # never overshoot the gene
            return False
        # Extend the terminal exon along the circle; it may now cross the origin
        # and become two arcs.
        if st0 == 1:
            ns, ne = s0, e0 + grow
        else:
            ns, ne = s0 - grow, e0
        arcs = _coords.region_arcs(ns % glen, (ne % glen) or glen, glen) \
            if (ns < 0 or ne > glen) else [(ns, ne)]
        exons = sorted([(a, b) for a, b, _ in parts[:-1]] + arcs)
        feat.exons = exons
        if st0 == 1:
            feat.end = (e0 + grow) % glen or glen
        else:
            feat.start = (s0 - grow) % glen
        return True
    return False


def refine_all(annotations, genome_seq, gene_catalog=None):
    gene_catalog = gene_catalog or {}
    n = 0
    for f in annotations:
        if getattr(f, "gene_type", None) != "CDS": continue
        if getattr(f, "exon_strands", None): continue
        if not getattr(f, "exons", None) or len(f.exons) < 2: continue
        try:
            if _refine_one(f, genome_seq): n += 1
        except Exception as exc:                      # noqa: BLE001
            # One malformed feature must not abort the whole refinement pass, but
            # a bug in _refine_one was previously invisible: every gene it threw
            # on was simply left unrefined with no trace.
            print("      WARNING: splice refinement failed for %s (%s: %s)"
                  % (getattr(f, "gene_name", "?"), type(exc).__name__, exc))
    # Terminal-stop snap — all CDS (single- and multi-exon), skip pseudogenes and
    # trans-spliced features. Runs after junction refinement so it snaps the final
    # exon end.
    for f in annotations:
        if getattr(f, "gene_type", None) != "CDS": continue
        if getattr(f, "exon_strands", None): continue
        if getattr(f, "is_pseudogene", False): continue
        try:
            _snap_terminal_stop(f, genome_seq, gene_catalog)
        except Exception as exc:                      # noqa: BLE001
            print("      WARNING: terminal-stop snap failed for %s (%s: %s)"
                  % (getattr(f, "gene_name", "?"), type(exc).__name__, exc))
    return n
