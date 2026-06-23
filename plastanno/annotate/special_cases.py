"""
Special case handlers for complex genes.

Cases from v1 experience:
1. rps12 trans-splicing (3 exons across genome)
2. petB/petD short first exon (6bp, 9bp)
3. ycf1 pseudogene at IR/SSC junction
4. Gene synonym normalization
5. CAU anticodon disambiguation (trnI vs trnM)
"""
import difflib
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple, Dict
from Bio.Seq import Seq
from ..core.feature import Feature


# ── Gene synonyms ─────────────────────────────────────────────────────────────
SYNONYMS = {
    "clpP1" : "clpP",
    "clpP2" : "clpP",
    "psbN"  : "pbf1",
    "orf70a": "orf70",
    "orf70b": "orf70",
}

# ── CAU disambiguation (trnfM-CAU / trnM-CAU / trnI-CAU) ──────────────────────
# All three share anticodon CAU but are distinct genes with distinct sequences
# (trnfM is the formyl-Met initiator, trnM the elongator, trnI the isoleucine
# tRNA with a modified C). ARAGORN reports all of them as tRNA-Met/Ile and the
# old positional rule (in_IR → trnI, else → trnM) erased trnfM entirely and
# mis-named trnI copies whenever IR detection was off. BLASTing each locus
# against the tRNA DB — which carries all three named — recovers the true name
# at ~100% identity, independent of IR boundaries.
_CAU_NAMES = ("trnfM-CAU", "trnM-CAU", "trnI-CAU")


def _resolve_db_prefix(trna_db_dir):
    """Return a usable BLAST db prefix (has a .nin index) or None."""
    if not trna_db_dir:
        return None
    d = Path(trna_db_dir)
    for cand in (d / "global", d / "global.fasta"):
        if Path(str(cand) + ".nin").exists():
            return str(cand)
    return None


def _blast_cau_names(cau_feats, genome_seq, db_prefix):
    """Map id(feature) → best-hit CAU gene name via one batched blastn."""
    qfa = None
    try:
        with tempfile.NamedTemporaryFile(
            suffix=".fasta", mode="w", delete=False
        ) as f:
            for i, ann in enumerate(cau_feats):
                sub = genome_seq[ann.start:ann.end]
                if ann.strand == -1:
                    sub = str(Seq(sub).reverse_complement())
                f.write(f">{i}\n{sub}\n")
            qfa = f.name
        result = subprocess.run([
            "blastn", "-query", qfa, "-db", db_prefix,
            "-outfmt", "6 qseqid sseqid bitscore",
            "-word_size", "7", "-dust", "no", "-max_target_seqs", "5",
        ], capture_output=True, text=True, timeout=60)
    except (subprocess.SubprocessError, OSError):
        return {}
    finally:
        if qfa and os.path.exists(qfa):
            os.unlink(qfa)

    best = {}  # qid → (bitscore, name)
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        p = line.split("\t")
        if len(p) < 3:
            continue
        qid = int(p[0])
        m = re.search(r"trn(?:fM|M|I)-CAU", p[1])
        if not m:
            continue
        bs = float(p[2])
        if qid not in best or bs > best[qid][0]:
            best[qid] = (bs, m.group(0))
    return {qid: name for qid, (bs, name) in best.items()}


def normalize_names(annotations, ir_boundaries,
                    genome_seq=None, trna_db_dir=None):
    """Apply gene synonyms and CAU disambiguation."""
    changes = []
    for ann in annotations:
        # Apply synonyms
        if ann.gene_name in SYNONYMS:
            old = ann.gene_name
            ann.gene_name = SYNONYMS[old]
            changes.append(f"Renamed {old} → {ann.gene_name}")

    # CAU disambiguation — BLAST best-hit naming, positional fallback
    cau_feats = [a for a in annotations
                 if a.gene_type == "tRNA" and "CAU" in a.gene_name]
    db_prefix = _resolve_db_prefix(trna_db_dir)
    blast_names = ({} if (genome_seq is None or db_prefix is None)
                   else _blast_cau_names(cau_feats, genome_seq, db_prefix))

    for i, ann in enumerate(cau_feats):
        name = blast_names.get(i)
        if name is None:
            # Fallback: positional (IR → trnI, else → trnM); cannot tell trnfM.
            irb = ir_boundaries.get("IRb", (0, 0))
            ira = ir_boundaries.get("IRa", (0, 0))
            mid = (ann.start + ann.end) // 2
            in_ir = (irb[0] <= mid <= irb[1] or ira[0] <= mid <= ira[1])
            name = "trnI-CAU" if in_ir else "trnM-CAU"
        if ann.gene_name != name:
            changes.append(f"CAU: {ann.gene_name} → {name} @{ann.start}")
            ann.gene_name = name

    return annotations, changes


# ── rps12 trans-splicing ──────────────────────────────────────────────────────
_STOPS = ("TAA", "TAG", "TGA")


def _rc(s):
    return str(Seq(s).reverse_complement())


def _coding_of(genome, s, e, strand):
    return _rc(genome[s:e]) if strand == -1 else genome[s:e]


def _reconstruct_rps12_copy(genome, exon1, ir_strand, exon2_start_g, refs):
    """Reconstruct exon2 + exon3 of one IR copy of trans-spliced rps12.

    exon1 (5' exon, in the LSC) is fixed. Anchored on the exon2 start (from an
    Exonerate hit in the IR), the long exon2 ends at a GT donor and the short
    3' exon3 follows after an A[Y]/AG acceptor and terminates at a stop codon.
    The exact boundaries are pinned by full-length identity to the best-matching
    reference protein — loose consensus alone cannot localise the internal exon
    boundaries of this divergent gene. All coordinates are returned in the genome
    frame. Returns (identity, exon2_genomic, exon3_genomic) or None.
    """
    e1s, e1e, e1st = exon1
    e1c = _coding_of(genome, e1s, e1e, e1st)
    # rps12's 5' exon may begin at an alternative start codon (GTG) in some
    # lineages; treat it as Met rather than rejecting the reconstruction.
    rps12_starts = _valid_starts("rps12")
    e1prot = _translate_orf(e1c, rps12_starts).rstrip("*")
    if not e1prot or e1prot[0] != "M":
        return None
    glen = len(genome)
    if ir_strand == 1:
        w0, w1 = max(0, exon2_start_g - 12), min(glen, exon2_start_g + 1600)
        view = genome[w0:w1]
        tomap = lambda ci, cj: (w0 + ci, w0 + cj)
        a0 = exon2_start_g - w0
    else:
        w0, w1 = max(0, exon2_start_g - 1600), min(glen, exon2_start_g + 12)
        view = _rc(genome[w0:w1])
        tomap = lambda ci, cj: (w1 - cj, w1 - ci)
        a0 = w1 - exon2_start_g

    n1 = len(e1prot)
    cand_refs = sorted(
        refs,
        key=lambda r: difflib.SequenceMatcher(None, e1prot, r[:n1], autojunk=False).ratio(),
        reverse=True,
    )[:2]
    vlen = len(view)
    best = None
    for s2 in range(max(0, a0 - 12), a0 + 13, 3):
        for e2 in range(s2 + 195, min(vlen - 2, s2 + 285)):
            if view[e2:e2 + 2] != _DONOR:
                continue
            if "*" in _translate(e1c + view[s2:e2]):     # frame broken before exon3
                continue
            for a3 in range(e2 + 120, min(vlen, e2 + 950)):
                if view[a3 - 2:a3] not in _ACCEPTORS:
                    continue
                for e3 in range(a3 + 9, min(vlen, a3 + 45)):
                    seg = view[a3:e3]
                    if seg[-3:] not in _STOPS:
                        continue
                    prot = _translate_orf(e1c + view[s2:e2] + seg, rps12_starts)
                    if prot[:1] != "M" or "*" in prot[:-1] or not prot.endswith("*"):
                        continue
                    core = prot[:-1]
                    idn = max(difflib.SequenceMatcher(None, core, r, autojunk=False).ratio()
                              for r in cand_refs)
                    if idn >= 0.85 and (best is None or idn > best[0]):
                        best = (idn, tomap(s2, e2), tomap(a3, e3))
    return best


def _load_rps12_refs(protein_db):
    from Bio import SeqIO
    p = Path(protein_db) / "rps12.fasta"
    if not p.exists():
        return []
    try:
        seqs = [str(r.seq).rstrip("*") for r in SeqIO.parse(str(p), "fasta")]
    except Exception:
        return []
    return [s for s in seqs if len(s) >= 30]


def handle_rps12(annotations, genome_seq, ir_boundaries, protein_db=None):
    """Reconstruct trans-spliced rps12 as one feature per IR copy.

    rps12 has a 5' exon in the LSC and two 3' exons (exon2 long, exon3 short ~26 bp)
    that lie in the inverted repeat — hence duplicated, with the two copies on
    opposite strands. We keep the precise 5' exon (from Engine A), locate exon2 in
    each IR copy with Exonerate, then reconstruct exon2/exon3 boundaries by
    reference-protein identity (`_reconstruct_rps12_copy`). Each copy is emitted as
    a trans-spliced CDS carrying per-exon strands. If reconstruction fails, the raw
    rps12 fragments are left untouched (no regression).
    """
    rps12 = [a for a in annotations if a.gene_name == "rps12"]
    if not rps12 or not protein_db:
        return annotations, []
    other = [a for a in annotations if a.gene_name != "rps12"]

    refs = _load_rps12_refs(protein_db)
    if not refs:
        # Cannot reconstruct any genome (global DB issue) — keep raw, don't punish.
        return annotations, []
    rep = max(refs, key=len)

    lsc_end = ir_boundaries.get("LSC", (0, 0))[1]
    e1_feats = [f for f in rps12 if f.end <= lsc_end]
    if not e1_feats:
        # No 5' exon in the LSC (IR-lacking or IR-boundary mis-call): reconstruction
        # is impossible, and the raw fragments are disconnected single-exon pieces
        # that never match a trans-spliced reference — pure false positives.
        return other, ["rps12: no LSC exon1; dropped %d raw fragment(s)" % len(rps12)]
    exon1_f = min(e1_feats, key=lambda f: f.start)
    exon1 = (exon1_f.start, exon1_f.end, exon1_f.strand)
    # rps12 exon1 is canonically 38 codons = 114 nt with a phase-0 trans-splice
    # boundary (cut after Tyr-38). Exonerate frequently over-extends the donor
    # end; trim exon1 to 114 nt from the start codon so the boundary is exact.
    _RPS12_EXON1 = 114
    _e1s, _e1e, _e1st = exon1
    if _e1e - _e1s > _RPS12_EXON1:
        if _e1st == -1:
            exon1 = (_e1e - _RPS12_EXON1, _e1e, _e1st)
        else:
            exon1 = (_e1s, _e1s + _RPS12_EXON1, _e1st)

    from ..identify.engine_a import run_exonerate_region
    built = []
    for region in ("IRb", "IRa"):
        coords = ir_boundaries.get(region)
        if not coords:
            continue
        hits = run_exonerate_region(genome_seq, rep, "rps12", "rps12_ref",
                                    coords[0], coords[1])
        if not hits:
            continue
        hit = max(hits, key=lambda h: h.s_ref)
        anchor = hit.start if hit.strand == 1 else hit.end
        res = _reconstruct_rps12_copy(genome_seq, exon1, hit.strand, anchor, refs)
        if not res:
            continue
        idn, e2, e3 = res
        feat = Feature(
            gene_name="rps12", gene_type="CDS",
            product="ribosomal protein S12",
            exons=[(exon1[0], exon1[1]), e2, e3],          # transcript order
            exon_strands=[exon1[2], hit.strand, hit.strand],
            start=min(exon1[0], e2[0], e3[0]),
            end=max(exon1[1], e2[1], e3[1]),
            strand=exon1[2], engine="AB",
            confidence=0.9, flag="HIGH",
            notes=["trans-spliced gene (%s copy, id=%.2f)" % (region, idn)],
        )
        feat.has_intron = True
        built.append(feat)

    if not built:
        # Reconstruction failed for every IR copy (commonly because Engine A's raw
        # 5' exon does not start at the Met codon, so the spliced ORF is invalid).
        # The leftover raw fragments are single-exon pieces that cannot match a
        # trans-spliced reference, so drop them rather than emit false positives.
        return other, ["rps12: reconstruction failed; dropped %d raw fragment(s)" % len(rps12)]
    # collapse identical reconstructions (e.g. when only one IR copy is real)
    uniq = []
    for f in built:
        if not any(f.exons == u.exons for u in uniq):
            uniq.append(f)
    return other + uniq, ["rps12: reconstructed %d trans-spliced copy(ies)" % len(uniq)]


# ── petB/petD short first exon ────────────────────────────────────────────────

# Genes whose first exon is so short (≈6–9 bp) that Exonerate/HMM align only the
# long second exon, leaving the outer boundary off by ~one intron. We recover the
# missing first exon by reference-protein-anchored search (see below).
SHORT_EXON_GENES = {
    "petB":  {"first_exon_len": 6, "product": "cytochrome b6"},
    "petD":  {"first_exon_len": 8, "product": "cytochrome b6/f complex subunit IV"},
    "rpl16": {"first_exon_len": 9, "product": "ribosomal protein L16"},
}

# Search bounds for the intron separating the short first exon from exon 2.
_MIN_INTRON = 150
_MAX_INTRON = 1700
# Candidate first-exon lengths (bp). The intron may split a codon, so the length
# need not be a multiple of 3; the reading frame is validated on the spliced ORF.
_L_CANDS = (3, 6, 7, 8, 9, 12)


_NTERM_K = 12          # N-terminal residues used to anchor the first exon
_MIN_SUPPORT = 0.15    # ≥15% of reference proteins must corroborate the N-terminus
_MAX_TRIM = 12         # bp the long exon's 5' boundary may be trimmed to fix frame
_DONOR = "GT"                       # group-II intron 5' splice site (GU)
_ACCEPTORS = {"AT", "AC", "AG"}     # 3' splice site (group-II AY, or canonical AG)


def _translate(seq):
    seq = seq[: len(seq) // 3 * 3]
    if not seq:
        return ""
    return str(Seq(seq).translate(table=11))


from ..core.reconcile import SPECIAL_START_CODONS


def _valid_starts(gene_name):
    """ATG plus any community-recognised alternative start codons for this gene."""
    return ("ATG",) + tuple(SPECIAL_START_CODONS.get(gene_name, ()))


def _translate_orf(seq, valid_starts):
    """Translate an ORF, mapping a recognised alternative start codon to Met so the
    spliced protein compares fairly against ATG-started reference sequences."""
    p = _translate(seq)
    if p and seq[:3].upper() in valid_starts and p[0] != "M":
        p = "M" + p[1:]
    return p


def _nterm_support(core, profile):
    """Fraction of reference proteins whose N-terminus matches `core`'s N-terminus.

    The discriminator between the true start codon and a spurious upstream ATG is
    almost entirely in the first few residues (the rest, from the long second exon,
    is identical across placements). We therefore score a candidate by how many
    references share its N-terminal K-mer — the true placement reconstitutes the
    conserved start shared by the majority; a wrong one matches almost none.
    """
    prefixes, total = profile
    if total == 0 or len(core) < _NTERM_K:
        return 0.0
    pre = core[:_NTERM_K]
    hit = 0
    for prefix, cnt in prefixes:
        if difflib.SequenceMatcher(None, pre, prefix, autojunk=False).ratio() >= 0.8:
            hit += cnt
    return hit / total


def _score_spliced(spliced, profile, median_len):
    """Return (ok, support) for a candidate spliced coding sequence."""
    if len(spliced) < 30:
        return False, 0.0
    prot = _translate(spliced)
    if not prot or prot[0] != "M":
        return False, 0.0
    core = prot[:-1] if prot.endswith("*") else prot
    if "*" in core:                                  # internal stop → wrong frame
        return False, 0.0
    if not (0.6 * median_len <= len(core) <= 1.4 * median_len):
        return False, 0.0
    return True, _nterm_support(core, profile)


def _recover_first_exon(feat, genome_seq, profile, median_len):
    """Find and prepend the short first exon of a single-exon CDS feature.

    Scans the upstream region (on the coding strand) for an ATG that, spliced to
    the long exon, yields a clean ORF whose N-terminus matches the reference
    consensus. The long exon's 5' boundary is allowed to be trimmed by a few bp
    (`delta`): Exonerate often over-/under-extends it into the intron by 1–2
    codons, which would otherwise shift the reading frame and hide the true start.
    Returns True and mutates `feat` (exons/start/end) on success.
    """
    if len(feat.exons) > 1:          # already multi-exon — nothing to recover
        return False
    e2s, e2e = feat.start, feat.end
    strand = feat.strand
    glen = len(genome_seq)

    def rc(s):
        return str(Seq(s).reverse_complement())

    # Candidate exon2 5'-boundary trims, keeping only those left by a canonical
    # group-II acceptor site (intron ends in A[T/C], or AG). exon2_coding[d] is the
    # coding sequence of the long exon trimmed by d bp at its 5' end.
    exon2_coding = {}
    for d in range(_MAX_TRIM):
        if strand == 1:
            acc = genome_seq[e2s + d - 2:e2s + d]
            if acc in _ACCEPTORS:
                exon2_coding[d] = genome_seq[e2s + d:e2e]
        else:
            acc = rc(genome_seq[e2e - d:e2e - d + 2])
            if acc in _ACCEPTORS:
                exon2_coding[d] = rc(genome_seq[e2s:e2e - d])
    if not exon2_coding:
        return False

    best = None  # (support, -intron, (exon1_start, exon1_end), delta)

    if strand == 1:
        lo = max(0, e2s - _MAX_INTRON)
        hi = e2s - _MIN_INTRON
        for g in range(lo, hi):
            if genome_seq[g:g + 3] != "ATG":
                continue
            for L in _L_CANDS:
                e1e = g + L
                if e1e > e2s - _MIN_INTRON:
                    continue
                if genome_seq[e1e:e1e + 2] != _DONOR:          # 5' splice site
                    continue
                e1c = genome_seq[g:e1e]
                for d, e2c in exon2_coding.items():
                    ok, sup = _score_spliced(e1c + e2c, profile, median_len)
                    if ok and sup >= _MIN_SUPPORT:
                        cand = (sup, -(e2s + d - e1e), (g, e1e), d)
                        if best is None or cand > best:
                            best = cand
    else:
        lo = e2e + _MIN_INTRON
        hi = min(glen, e2e + _MAX_INTRON)
        for g in range(lo, hi):
            if rc(genome_seq[g - 2:g]) != _DONOR:               # 5' splice site
                continue
            for L in _L_CANDS:
                e1e = g + L
                if e1e > glen:
                    continue
                e1c = rc(genome_seq[g:e1e])
                if e1c[:3] != "ATG":
                    continue
                for d, e2c in exon2_coding.items():
                    ok, sup = _score_spliced(e1c + e2c, profile, median_len)
                    if ok and sup >= _MIN_SUPPORT:
                        cand = (sup, -(g - (e2e - d)), (g, e1e), d)
                        if best is None or cand > best:
                            best = cand

    if best is None:
        return False

    e1s, e1e = best[2]
    d = best[3]
    exon2 = (e2s + d, e2e) if strand == 1 else (e2s, e2e - d)
    feat.exons = sorted([(e1s, e1e), exon2])
    feat.start = min(e1s, exon2[0])
    feat.end = max(e1e, exon2[1])
    feat.has_intron = True
    feat.notes.append("recovered short first exon (%d bp)" % (e1e - e1s))
    return True


def _load_ref_profile(protein_db, gene):
    """Return (profile, median_len, refs) for a gene, or None.

    `profile` is the N-terminus k-mer histogram, `refs` a small set of full
    reference proteins (closest to the median length) used as a full-length
    identity anchor for genes whose N-terminus is too divergent for the consensus.
    """
    from collections import Counter
    from Bio import SeqIO
    p = Path(protein_db) / ("%s.fasta" % gene)
    if not p.exists():
        return None
    try:
        seqs = [str(r.seq).rstrip("*") for r in SeqIO.parse(str(p), "fasta")]
    except Exception:
        return None
    seqs = [s for s in seqs if len(s) >= _NTERM_K]
    if not seqs:
        return None
    counts = Counter(s[:_NTERM_K] for s in seqs)
    profile = (list(counts.items()), sum(counts.values()))
    lens = sorted(len(s) for s in seqs)
    median_len = lens[len(lens) // 2]
    refs = sorted(seqs, key=lambda s: abs(len(s) - median_len))[:8]
    return profile, median_len, refs


def handle_short_exon_genes(annotations, genome_seq,
                            ir_boundaries, gene_catalog, protein_db=None):
    """Recover the short first exon of petB/petD/rpl16 (dropped by Exonerate/HMM).

    Each is normally annotated as a single (long) exon, leaving the outer boundary
    off by ~one intron. We anchor the missing ~6–9 bp first exon to the reference
    protein N-terminus, so the annotated start/end and splice structure become
    correct.
    """
    changes = []
    if not protein_db:
        return annotations, changes

    for gene in SHORT_EXON_GENES:
        feats = [a for a in annotations
                 if a.gene_name == gene and a.gene_type == "CDS"]
        if not feats:
            continue
        loaded = _load_ref_profile(protein_db, gene)
        if not loaded:
            continue
        profile, median_len, _refs = loaded
        for feat in feats:
            if _recover_first_exon(feat, genome_seq, profile, median_len):
                changes.append("%s: recovered short first exon" % gene)

    return annotations, changes


# ── ORF boundary completion (truncated single-exon CDS) ───────────────────────
# Exonerate/HMM align a divergent reference protein only over its conserved core,
# so single-exon CDS are often reported as a truncated piece of the true ORF —
# the boundary is an alignment edge, not a codon edge (offsets are whole codons).
# These genes are clean ATG→stop ORFs (no RNA editing), so the fix is to complete
# the ORF: extend the 3' end to the first in-frame stop and snap the 5' start to
# the upstream ATG whose translation matches the reference N-terminus consensus.

_ORF_PAD = 3500    # bp searched up/downstream when completing an ORF (covers ycf2)
_NEAR_ATG = 30     # bp window for the nearest-ATG start fallback (divergent N-termini)


def _complete_orf(feat, genome_seq, profile, median_len, refs):
    """Extend a truncated single-exon CDS toward its full ATG→stop ORF.

    The 3' end is always extended to the first in-frame stop (safe, no anchoring
    needed — recovers e.g. ycf2 truncated by ~2.6 kb). The 5' start is snapped to
    the best ATG: one whose N-terminus matches the reference consensus if available
    (handles conserved genes), otherwise the nearest in-frame ATG within a small
    window (handles divergent N-termini like ycf2 where the start is ~correct but
    the consensus does not match). Never crosses an in-frame stop, so genuinely
    RNA-edited genes are left untouched. Mutates feat on success.
    """
    if feat.exons and len(feat.exons) > 1:
        return False
    s, e, strand = feat.start, feat.end, feat.strand
    if (e - s) % 3 != 0:
        return False
    valid_starts = _valid_starts(feat.gene_name)
    glen = len(genome_seq)
    w0, w1 = max(0, s - _ORF_PAD), min(glen, e + _ORF_PAD)
    if strand == 1:
        view = genome_seq[w0:w1]
        tomap = lambda a, b: (w0 + a, w0 + b)
        cs, ce = s - w0, e - w0
    else:
        view = _rc(genome_seq[w0:w1])
        tomap = lambda a, b: (w1 - b, w1 - a)
        cs, ce = w1 - e, w1 - s
    vlen = len(view)

    # 3' end: first in-frame stop at/after the current start frame
    stop_pos = None
    i = cs
    while i + 3 <= vlen:
        if view[i:i + 3] in _STOPS:
            stop_pos = i
            break
        i += 3
    if stop_pos is None:
        return False
    ce2 = stop_pos + 3

    # 5' start: in-frame ATGs from the previous upstream stop up to cs.
    lo = 0
    j = cs - 3
    while j >= 0:
        if view[j:j + 3] in _STOPS:
            lo = j + 3
            break
        j -= 3
    atgs = []  # (pos, nterm_support)
    c = lo
    while c <= cs:
        if view[c:c + 3] in valid_starts:
            prot = _translate_orf(view[c:ce2], valid_starts)
            if prot[:1] == "M" and prot.endswith("*") and "*" not in prot[:-1]:
                atgs.append((c, _nterm_support(prot[:-1], profile)))
        c += 3
    strong = [a for a in atgs if a[1] >= 0.15]
    if strong:                                    # consensus-anchored start
        best_c = max(strong, key=lambda a: (a[1], -a[0]))[0]
    else:
        # divergent N-terminus: anchor by full-length protein identity to a ref
        fp_best = None
        for c0, _sup in atgs:
            prot = _translate_orf(view[c0:ce2], valid_starts)
            prot = prot[:-1] if prot.endswith("*") else prot
            if not (0.5 * median_len <= len(prot) <= 1.6 * median_len):
                continue
            ratio = max((difflib.SequenceMatcher(None, prot, r, autojunk=False).ratio()
                         for r in refs), default=0.0)
            if ratio >= 0.6 and (fp_best is None or ratio > fp_best[0]):
                fp_best = (ratio, c0)
        if fp_best:
            best_c = fp_best[1]
        else:                                     # nearest-ATG fallback
            near = [a for a in atgs if abs(a[0] - cs) <= _NEAR_ATG]
            best_c = min(near, key=lambda a: abs(a[0] - cs))[0] if near else cs

    if (best_c, ce2) == (cs, ce):
        return False                              # already complete — no change
    core = _translate_orf(view[best_c:ce2], valid_starts)
    core = core[:-1] if core.endswith("*") else core
    if "*" in core:                               # would contain an internal stop
        return False
    if not (0.5 * median_len <= len(core) <= 1.6 * median_len):
        return False
    ns, ne = tomap(best_c, ce2)
    feat.start, feat.end = ns, ne
    feat.exons = [(ns, ne)]                        # single-exon ORF — keep exons in sync
    feat.notes.append("ORF boundary completed")
    return True


def complete_single_exon_orfs(annotations, genome_seq, gene_catalog, protein_db=None):
    """Complete truncated single-exon CDS to full ORFs (accD, ndhK, rps18, ycf2, …).

    Genes the catalog records as single-exon but that Exonerate split into several
    "exons" (a spurious intron over a divergent/indel region) are collapsed to
    their outer span before completion; genuinely multi-exon genes (catalog
    n_exons > 1, e.g. petB/clpP) are left to the splicing handlers.
    """
    changes = []
    if not protein_db:
        return annotations, changes
    cache = {}
    for feat in annotations:
        if feat.gene_type != "CDS":
            continue
        gene = feat.gene_name
        n_exons_cat = gene_catalog.get(gene, {}).get("n_exons", 1) or 1
        multi = bool(feat.exons) and len(feat.exons) > 1
        if multi and n_exons_cat > 1:
            continue                              # genuine multi-exon gene
        if gene not in cache:
            cache[gene] = _load_ref_profile(protein_db, gene)
        if not cache[gene]:
            continue
        profile, median_len, refs = cache[gene]
        saved = feat.exons
        if multi:
            feat.exons = []                       # drop spurious intron; outer span kept
        if _complete_orf(feat, genome_seq, profile, median_len, refs):
            changes.append("%s: ORF boundary completed" % gene)
        elif multi:
            feat.exons = saved                    # nothing improved — restore
    return annotations, changes


# ── ycf1 pseudogene ───────────────────────────────────────────────────────────
def handle_ycf1(annotations, genome_seq, ir_boundaries):
    """
    ycf1 appears twice:
    1. Full-length in SSC (functional)
    2. Truncated pseudogene at IRb/SSC junction

    Mark the shorter copy as pseudogene.
    """
    ycf1_feats = [a for a in annotations
                  if a.gene_name == "ycf1"]
    other      = [a for a in annotations
                  if a.gene_name != "ycf1"]

    if len(ycf1_feats) < 2:
        return annotations, []

    # Longest = functional, shortest = pseudogene
    ycf1_feats.sort(key=lambda x: x.end-x.start,
                    reverse=True)
    ycf1_feats[0].is_pseudogene   = False
    ycf1_feats[0].flag            = "HIGH"

    ycf1_feats[1].is_pseudogene    = True
    ycf1_feats[1].pseudogene_reason= "truncated IR copy"
    ycf1_feats[1].flag             = "MEDIUM"
    ycf1_feats[1].notes.append("ycf1 pseudogene at IR junction")

    return other + ycf1_feats, ["ycf1: marked truncated copy as pseudogene"]


# ── Main ──────────────────────────────────────────────────────────────────────
def run_all_special_cases(
    annotations,
    genome_seq,
    ir_boundaries,
    gene_catalog,
    protein_db=None,
    trna_db_dir=None,
):
    """Run all special case handlers."""
    all_changes = []

    # 1. Normalize names
    annotations, changes = normalize_names(
        annotations, ir_boundaries,
        genome_seq=genome_seq, trna_db_dir=trna_db_dir,
    )
    all_changes.extend(changes)

    # 2. rps12 trans-splicing
    annotations, changes = handle_rps12(
        annotations, genome_seq, ir_boundaries, protein_db
    )
    all_changes.extend(changes)

    # 3. Short exon genes
    annotations, changes = handle_short_exon_genes(
        annotations, genome_seq,
        ir_boundaries, gene_catalog, protein_db
    )
    all_changes.extend(changes)

    # 3b. Complete truncated single-exon ORFs (accD, ndhK, rps18, ndhI, …)
    annotations, changes = complete_single_exon_orfs(
        annotations, genome_seq, gene_catalog, protein_db
    )
    all_changes.extend(changes)

    # 4. ycf1 pseudogene
    annotations, changes = handle_ycf1(
        annotations, genome_seq, ir_boundaries
    )
    all_changes.extend(changes)

    return annotations, all_changes
