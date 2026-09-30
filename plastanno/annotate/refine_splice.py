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

# Why the exon panel could not be searched, once said. Without the panel every panel gene keeps
# its engine junctions, and that used to happen without a word: a database that BLAST cannot
# read (a version-5 rebuild missing its index files, a truncated install) made _blast_exons
# return nothing, and junction refinement was skipped for every gene while the run looked normal.
_PANEL_PROBLEM = []


def _panel_unavailable(why):
    if not _PANEL_PROBLEM:
        _PANEL_PROBLEM.append(why)
        print("      WARNING: splice-junction refinement is skipped for every panel gene: %s" % why)
    return {}


def _blast_exons(region, gene):
    if not os.path.exists(_BDB + ".nin"):
        return _panel_unavailable("the exon panel database is missing (%s.nin)" % _BDB)
    qf = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    qf.write(">q\n" + region + "\n"); qf.close()
    try:
        r = subprocess.run(
            ["blastn", "-task", "blastn-short", "-query", qf.name, "-db", _BDB,
             "-outfmt", "6 sseqid qstart qend sstart send slen bitscore", "-strand", "plus",
             "-word_size", "7", "-dust", "no", "-evalue", "1", "-max_target_seqs", "300"],
            capture_output=True, text=True, timeout=60)
    except Exception as exc:                          # noqa: BLE001
        os.unlink(qf.name)
        return _panel_unavailable("blastn could not be run (%s: %s)" % (type(exc).__name__, exc))
    os.unlink(qf.name)
    if r.returncode != 0:
        err = (r.stderr or "").strip().splitlines()
        return _panel_unavailable("BLAST could not search the exon panel (%s)"
                                  % (err[0] if err else "exit %d" % r.returncode))
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
# How far to look. On development genomes, the first in-frame stop after a CDS that lacked one
# was the reference's 3' end in almost every case, often well beyond 3 codons, and a read-through
# past the real stop was never more than a few codons.
_OLD_SNAP_CODONS = 3          # what 3.0.1 did silently, with its length guard
_EXTEND_MAX_CODONS = 100
_TRIM_MAX_CODONS = 5


def _snap_terminal_stop(feat, g, gene_catalog):
    """End a CDS at its first in-frame stop codon, and say so when that was a real repair.

    On the spliced CDS read in transcription order, in this order:

      1 trim      the first in-frame stop before the last codon lies within the last
                  _TRIM_MAX_CODONS codons: the CDS read through its stop; end it there. Marked.
      2           the last codon is a stop: nothing to do.
      3 old snap  a stop within _OLD_SNAP_CODONS codons that keeps the CDS within 1.2x the
                  catalog's expected length: extend to it silently. This is exactly what 3.0.1
                  did, whatever lies upstream (a gene whose reference carries an internal stop,
                  like some rps11, still gets its adjacent stop).
      4 complete  no internal stop at all: extend to the first in-frame stop within
                  _EXTEND_MAX_CODONS codons, read round the origin. Marked.
      5           otherwise nothing: an internal stop further upstream is a frame or pseudogene
                  problem, and a clean 3' end would only hide it.

    A marked CDS carries a note starting with finalize.COMPLETED_3P, and
    finalize.submission_check keeps it NEEDS_REVIEW even when its ORF is now valid. Completing
    a 3' end makes the ORF pass NCBI's validator, which would otherwise remove the warning from a
    gene model that is still wrong at its start or its splice sites. On development genomes that
    happened to most of the CDS repaired this way.
    """
    from ..core.finalize import COMPLETED_3P
    glen = len(g)
    # The coding sequence in TRANSCRIPTION order. Concatenating sorted(exons) puts
    # a wrapped gene's parts the wrong way round — the piece after the origin first
    # — so the frame, and therefore the terminal codon, were read from the wrong
    # place.
    seq = _coords.extract(g, feat, glen).upper()
    if len(seq) < 6 or len(seq) % 3 != 0:      # need an intact frame
        return False
    # The exons in transcription order; the LAST one carries the 3' end, wherever
    # that sits on the circle.
    parts = _coords.transcript_parts(feat, glen)
    if not parts:
        return False
    s0, e0, st0 = parts[-1]
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    internal = [i for i, c in enumerate(codons[:-1]) if c in _STOP]

    # 1 trim a short read-through past the first stop
    if internal and len(codons) - 1 - internal[0] <= _TRIM_MAX_CODONS:
        cut = len(codons) - 1 - internal[0]      # codons read past the first stop
        if 3 * cut >= e0 - s0:
            return False                         # the cut would leave the last exon
        if st0 == 1:
            last = (s0, e0 - 3 * cut)
            feat.end = e0 - 3 * cut
        else:
            last = (s0 + 3 * cut, e0)
            feat.start = s0 + 3 * cut
        feat.exons = sorted([(a, b) for a, b, _ in parts[:-1]] + [last])
        feat.notes.append(COMPLETED_3P + "3' end trimmed to the first in-frame stop "
                          "(%d codon(s) of read-through)" % cut)
        return True
    # 2 already terminated
    if codons[-1] in _STOP:
        return False

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

    def extend(k):
        grow = 3 * (k + 1)
        # Extend the terminal exon along the circle; it may now cross the origin
        # and become two arcs.
        if st0 == 1:
            ns, ne = s0, e0 + grow
        else:
            ns, ne = s0 - grow, e0
        arcs = _coords.region_arcs(ns % glen, (ne % glen) or glen, glen) \
            if (ns < 0 or ne > glen) else [(ns, ne)]
        feat.exons = sorted([(a, b) for a, b, _ in parts[:-1]] + arcs)
        if st0 == 1:
            feat.end = (e0 + grow) % glen or glen
        else:
            feat.start = (s0 - grow) % glen

    # 3 the old snap, unchanged: a stop within 3 codons, within the old length guard
    cur_len = len(seq)
    exp = (gene_catalog.get(feat.gene_name, {}) or {}).get("expected_len", 0)
    for k in range(_OLD_SNAP_CODONS):
        codon = codon_at(k)
        if len(codon) != 3:
            break
        if codon in _STOP:
            if not (exp and (cur_len + 3 * (k + 1)) / exp > 1.2):
                extend(k)
                return True
            break                                # guard blocked: leave it to step 4
    # 4 complete a truncated 3' end
    if internal:
        return False
    for k in range(_EXTEND_MAX_CODONS):
        codon = codon_at(k)
        if len(codon) != 3:
            break
        if codon in _STOP:
            extend(k)
            feat.notes.append(COMPLETED_3P + "3' end extended to the first in-frame stop "
                              "(%d codon(s))" % (k + 1))
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
    # Terminal stop — all CDS (single- and multi-exon), skip pseudogenes and
    # trans-spliced features: trim a short read-through, snap an adjacent stop, or
    # complete a truncated 3' end. Runs after junction refinement so it works on the
    # final exon end.
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
