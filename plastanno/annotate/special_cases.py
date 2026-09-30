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
import hashlib
import tempfile
from pathlib import Path
from typing import List, Tuple, Dict
from Bio.Seq import Seq
from ..core.feature import Feature
from ..core import trna_identity as _TI
from ..core import coords as _coords
from ..core import ambiguity as _ambiguity


# ── Gene synonyms ─────────────────────────────────────────────────────────────
# Output names follow the MAJORITY usage in GenBank plastomes: psbN, ycf3, ycf4
# and clpP, not the newer pbf1, pafI, pafII and clpP1 (first seen 2014-15, used
# by 23-28% of records submitted in 2024; parked/arbiter_oracle/
# gene_name_timeline.txt). This table once mapped psbN -> pbf1, the minority name.
#
# Applied here, in step 6 only. Applying it to the engines' calls before
# reconciliation was tried and rejected (parked/gene_name_policy/RESULT.md): it
# is not a rename -- Engine A's pbf1 then pairs with Engine B's psbN, and the pafI
# and pafII profiles' calls compete inside the ycf3 and ycf4 clusters -- and on
# 30 dev genomes it moved coordinates in 8 and lost an exact psbN.
SYNONYMS = {
    "clpP1" : "clpP",
    "clpP2" : "clpP",
    "pbf1"  : "psbN",
    "pafI"  : "ycf3",
    "pafII" : "ycf4",
    "orf70a": "orf70",
    "orf70b": "orf70",
}


def apply_synonyms(features):
    """Rename in place to the SYNONYMS target; return how many were renamed.

    The name the engine used is kept in the feature's notes, so the provenance
    written to /note still says where the call came from.
    """
    n = 0
    for f in features:
        new = SYNONYMS.get(f.gene_name)
        if new:
            f.notes.append("named %s; the engine called it %s" % (new, f.gene_name))
            f.gene_name = new
            n += 1
    return n

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
                sub = _coords.extract(genome_seq, ann, len(genome_seq) or None)
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
    renamed = []
    for ann in annotations:
        old = ann.gene_name
        if apply_synonyms([ann]):
            renamed.append(ann)
            changes.append(f"Renamed {old} → {ann.gene_name}")

    # A renamed call overlapping a call that already carries the target name is a
    # second call of one gene -- a pafI profile hit beside the ycf3 call, or
    # Engine A's pbf1 beside Engine B's psbN, which never paired under two names.
    # Keep the call that already had the majority name.
    if renamed:
        _glen = len(genome_seq or "") or None
        rid = {id(a) for a in renamed}
        drop = set()
        for r in renamed:
            for o in annotations:
                if (id(o) not in rid and o.gene_name == r.gene_name
                        and o.gene_type == r.gene_type
                        and _coords.overlap_bp(r, o, _glen) > 0):
                    drop.add(id(r))
                    changes.append(f"Dropped a renamed {r.gene_name} call at "
                                   f"{r.start}-{r.end}: it overlaps an existing "
                                   f"{o.gene_name} call")
                    break
        if drop:
            annotations = [a for a in annotations if id(a) not in drop]

    # CAU disambiguation — BLAST best-hit naming, positional fallback
    cau_feats = [a for a in annotations
                 if a.gene_type == "tRNA" and "CAU" in a.gene_name]
    db_prefix = _resolve_db_prefix(trna_db_dir)
    blast_names = ({} if (genome_seq is None or db_prefix is None)
                   else _blast_cau_names(cau_feats, genome_seq, db_prefix))

    for i, ann in enumerate(cau_feats):
        name = blast_names.get(i)
        source = _TI.BLAST_TRNA_DB if name is not None else _TI.POSITIONAL_IR
        if name is None:
            # Fallback: positional (IR → trnI, else → trnM); cannot tell trnfM.
            _glen = len(genome_seq or "") or None
            mid = _coords.circular_midpoint(ann, _glen)
            in_ir = any(_coords.region_contains(ir_boundaries[k], mid, _glen)
                        for k in ("IRb", "IRa") if ir_boundaries.get(k))
            name = "trnI-CAU" if in_ir else "trnM-CAU"
        # Same single identity API. The positional fallback is NOT authoritative
        # on identity -- it cannot distinguish trnfM from trnM -- so a
        # disagreement there yields CONFLICT rather than a forced rename.
        canon, note = _TI.resolve(_TI.from_gene_name(ann.gene_name),
                                  _TI.from_gene_name(name), source,
                                  detector_symbol=ann.gene_name)
        if (_TI.gene_name(canon) != ann.gene_name
                or _TI.product(canon) != ann.product):
            changes.append("CAU: %s → %s @%d (%s)"
                           % (ann.gene_name, _TI.gene_name(canon), ann.start,
                              source))
            _TI.apply_trna_identity(ann, canon, note)

    return annotations, changes


# ── rps12 trans-splicing ──────────────────────────────────────────────────────
_STOPS = ("TAA", "TAG", "TGA")


def _rc(s):
    return str(Seq(s).reverse_complement())


def _coding_of(genome, s, e, strand, genome_len=None):
    """The coding sequence of [s, e) read FORWARD around the circle.

    genome[s:e] is empty when the interval crosses the origin, which silently
    turned rps12's 5' exon into an empty protein and aborted the reconstruction.
    """
    L = genome_len or len(genome)
    seq = genome[s:e] if s < e else genome[s:L] + genome[0:e]
    return _rc(seq) if strand == -1 else seq


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
    e1c = _coding_of(genome, e1s, e1e, e1st, len(genome))
    # rps12's 5' exon may begin at an alternative start codon (GTG) in some
    # lineages; treat it as Met rather than rejecting the reconstruction.
    rps12_starts = _valid_starts("rps12")
    e1prot = _translate_orf(e1c, rps12_starts).rstrip("*")
    if not e1prot or e1prot[0] != "M":
        return None
    glen = len(genome)
    # The search window is taken around the circle. Clamping it at the ends of the
    # sequence made the window short, or empty, whenever the origin fell inside the
    # inverted repeat — and the origin is inside IRb for any genome deposited with
    # a different rotation, which is not rare.
    if ir_strand == 1:
        win = _coords.open_window(genome, (exon2_start_g - 12) % glen,
                                  (exon2_start_g + 1600) % glen, glen)
        view = win.seq
        tomap = lambda ci, cj: win.interval(ci, cj)
    else:
        win = _coords.open_window(genome, (exon2_start_g - 1600) % glen,
                                  (exon2_start_g + 12) % glen, glen)
        view = _rc(win.seq)
        _wl = len(win.seq)
        tomap = lambda ci, cj: win.interval(_wl - cj, _wl - ci)
    a0 = 12                       # the anchor sits 12 bases into the window

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

    # Membership in the LSC, read forward around the circle. The test used to be
    # `f.end <= LSC_end`, which is right only when the LSC does not cross the
    # origin. In a genome deposited with a different origin the LSC is recorded as
    # start > end, every exon fails `end <= LSC_end`, and the branch below then
    # deletes rps12 outright — both copies, in every such genome.
    glen = len(genome_seq) or None
    lsc = ir_boundaries.get("LSC")
    if lsc:
        e1_feats = [f for f in rps12 if _coords.feature_in_region(f, lsc, glen)]
    else:
        e1_feats = list(rps12)
    if not e1_feats:
        # No 5' exon in the LSC (IR-lacking or IR-boundary mis-call): reconstruction
        # is impossible, and the raw fragments are disconnected single-exon pieces
        # that never match a trans-spliced reference — pure false positives.
        return other, ["rps12: no LSC exon1; dropped %d raw fragment(s)" % len(rps12)]

    # Among those, prefer a fragment that actually looks like exon1: its own
    # translation starts with Met.
    #
    # Position alone is not enough. When ir_detector finds no inverted repeat the
    # pipeline sets ir_boundaries = {"LSC": (0, genome_len)}, so EVERY fragment is
    # "in the LSC" and the tiebreak below degenerates to lowest coordinate. On
    # NC_037507.1 that chose the fragment at 92..851, which translates to TTTPKKPN
    # and is the 3' half of the gene, over the real exon1 at 66686..66800
    # (MPTIQQLI). _reconstruct_rps12_copy then failed its own Met check and both
    # fragments were discarded. Four of the 30 development genomes failed this way
    # and all four have no detected IR — conifers and several bryophytes have
    # genuinely lost one repeat.
    #
    # If nothing qualifies the list is left alone, so a genome whose exon1
    # Exonerate truncated is no worse off than before.
    _starts = _valid_starts("rps12")
    _met = [f for f in e1_feats
            if (_translate_orf(_coding_of(genome_seq, f.start, f.end, f.strand,
                                          len(genome_seq)), _starts).rstrip("*") or " ")[0] == "M"]
    if _met:
        e1_feats = _met
    # "first" means closest to the start of the LSC around the circle, not the
    # lowest coordinate: with a wrapping LSC the lowest coordinate is the far end.
    _lsc_s = lsc[0] if lsc else 0
    exon1_f = min(e1_feats,
                  key=lambda f: (f.start - _lsc_s) % (glen or (max(f.end, 1))))
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
    # Where to look for the 3' exons. Normally one search per IR copy.
    #
    # A plastome with no inverted repeat has neither IRb nor IRa, so this loop
    # used to skip every iteration, `built` stayed empty, and the function
    # reported "reconstruction failed" without ever having attempted one. Four of
    # the 30 development genomes failed exactly this way, all of them IR-lacking:
    # conifers and several bryophytes have genuinely lost one repeat. In such a
    # genome the 3' exons are not duplicated, so one search over the whole
    # sequence is both correct and sufficient.
    _regions = [(r, ir_boundaries[r]) for r in ("IRb", "IRa") if ir_boundaries.get(r)]
    if not _regions:
        _regions = [("whole", (0, glen or len(genome_seq)))]
    for region, coords in _regions:
        hits = run_exonerate_region(genome_seq, rep, "rps12", "rps12_ref",
                                    coords[0], coords[1], genome_len=glen)
        # The anchor must be the 3' block, never exon1 itself. Searching a whole
        # IR-lacking genome returns both, and on NC_039155.1 exon1 scored higher
        # (0.895 against 0.880), so the resolver below chose it and the
        # reconstruction could not succeed. In the IR case exon1 lies in the LSC
        # and no IR hit overlaps it, so this filter is a no-op there.
        _e1s, _e1e = exon1[0], exon1[1]
        hits = [h for h in hits if not (h.start < _e1e and h.end > _e1s)]
        if not hits:
            continue
        # The shared ambiguity policy: rank on alignment identity alone, collapse
        # hits describing the same span, and if several still rank equally take a
        # stable representative rather than whichever Exonerate emitted first.
        import hashlib as _hl
        hit, _hit_alts, _hit_amb, _hit_ord = _ambiguity.resolve(
            hits,
            evidence=lambda h: h.s_ref,
            digest=lambda h: _hl.sha256(
                _coords.extract(genome_seq, h, glen).encode()).hexdigest(),
            model=lambda h: (h.start, h.end, h.strand))
        anchor = hit.start if hit.strand == 1 else hit.end
        res = _reconstruct_rps12_copy(genome_seq, exon1, hit.strand, anchor, refs)
        if not res:
            continue
        idn, e2, e3 = res
        # Any of the three exons may cross the origin, in which case it occupies
        # two arcs; each arc is listed separately, in transcription order, carrying
        # the strand of the exon it came from. A wrapped exon written as a single
        # (start, end) pair with start > end would be measured as empty by
        # everything downstream.
        parts, pstr = [], []
        for (ps, pe), pst in (((exon1[0], exon1[1]), exon1[2]),
                              (e2, hit.strand), (e3, hit.strand)):
            arcs = _coords.region_arcs(ps, pe, glen)
            # An exon split by the origin is transcribed tail-first on the plus
            # strand and head-first on the minus, exactly as core.coords orders a
            # wrapped feature's parts.
            if len(arcs) > 1 and pst != -1:
                arcs = arcs[::-1]
            parts.extend(arcs)
            pstr.extend([pst] * len(arcs))
        feat = Feature(
            gene_name="rps12", gene_type="CDS",
            product="ribosomal protein S12",
            exons=parts,                                   # transcript order
            exon_strands=pstr,
            start=min(a for a, _ in parts),
            end=max(b for _, b in parts),
            strand=exon1[2], engine="AB",
            confidence=0.9, flag="HIGH",
            notes=["trans-spliced gene (%s copy, id=%.2f)" % (region, idn)],
        )
        feat.has_intron = True
        feat.is_trans_spliced = True
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


def _ref_protein_file(protein_db, gene):
    """The file of a gene's reference proteins, or None.

    Step 6 renames a gene to its majority name before anything else, but the protein database
    keeps the proteins under the name the engines search with: psbN's are in pbf1.fasta. A
    lookup by the new name alone found nothing, so from 3.0.1 on ORF completion skipped every
    psbN, and a psbN call four codons short of its stop stayed short. The new name is read when
    its file exists; otherwise the file of a name SYNONYMS maps to it.
    """
    p = Path(protein_db) / ("%s.fasta" % gene)
    if p.exists():
        return p
    for old in sorted(o for o, new in SYNONYMS.items() if new == gene):
        q = Path(protein_db) / ("%s.fasta" % old)
        if q.exists():
            return q
    return None


def _load_ref_profile(protein_db, gene):
    """Return (profile, median_len, refs) for a gene, or None.

    `profile` is the N-terminus k-mer histogram, `refs` a small set of full
    reference proteins (closest to the median length) used as a full-length
    identity anchor for genes whose N-terminus is too divergent for the consensus.
    """
    from collections import Counter
    from Bio import SeqIO
    p = _ref_protein_file(protein_db, gene)
    if p is None:
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


def _missing_nterm_aa(feat, genome_seq, refs):
    """How many residues are missing from the N-terminus of a single-exon CDS.

    The feature holds only the long 3' exon, so its translation starts partway
    into the protein. Aligning it to the closest reference gives that offset, and
    3x it is where the missing 5' exon should be looked for. Returns None when no
    reference corroborates the fragment.
    """
    seq = _coding_of(genome_seq, feat.start, feat.end, feat.strand)
    prot = _translate(seq[: len(seq) // 3 * 3]).rstrip("*")
    prot = prot.split("*")[0] if "*" in prot else prot
    if len(prot) < 30 or not refs:
        return None
    best = None
    for r in refs:
        m = difflib.SequenceMatcher(None, prot, r, autojunk=False)
        a, b, size = m.find_longest_match(0, len(prot), 0, len(r))
        if size >= 25 and (best is None or size > best[1]):
            best = (b - a, size)                  # offset of prot's start inside r
    if best is None or best[0] <= 0:
        return None
    return best[0]


def recover_missing_first_exon(annotations, genome_seq, gene_catalog, protein_db=None):
    """Recover the 5' exon of a multi-exon CDS that was annotated as a single exon.

    `handle_short_exon_genes` covers only petB/petD/rpl16, and its length search
    (`_L_CANDS`) tops out at 12 bp, so it can only ever find a *tiny* first exon.
    But the dominant boundary error is not that the first exon is short — it is
    that it is MISSING. On DEV, 71.7% of all boundary errors are at the 5' end
    alone, and the biggest offenders are genes whose first exon is nowhere near
    tiny: rpoC1 (435 bp, 17 genomes), atpF (145 bp, 13), rps16 (40 bp, 11),
    clpP (217 bp, 6). Exonerate reports only exon 2 and the annotated start is then
    off by exon1 + intron — a median 5' error of +1123 bp for rpoC1.

    This generalises the recovery: it fires whenever the catalog says the gene has
    at least two exons but we produced one, and it sizes the search from the data
    rather than a fixed list — the fragment is aligned to its reference proteins to
    learn how many N-terminal residues are missing, and only that neighbourhood is
    scanned for a start codon with a canonical donor site whose spliced product is
    a clean ORF matching the reference N-terminus.
    """
    changes = []
    if not protein_db:
        return annotations, changes
    cache = {}
    for feat in annotations:
        if feat.gene_type != "CDS" or feat.is_pseudogene or feat.exon_strands:
            continue
        if feat.exons and len(feat.exons) > 1:
            continue                                  # already spliced
        gene = feat.gene_name
        if (gene_catalog.get(gene, {}) or {}).get("n_exons", 1) < 2:
            continue                                  # single-exon gene — nothing to add
        if gene in SHORT_EXON_GENES:
            continue                                  # handled by the tiny-exon path
        if gene not in cache:
            cache[gene] = _load_ref_profile(protein_db, gene)
        if not cache[gene]:
            continue
        profile, median_len, refs = cache[gene]
        k = _missing_nterm_aa(feat, genome_seq, refs)
        if k is None or k < 2:
            continue
        if _recover_first_exon_sized(feat, genome_seq, profile, median_len, 3 * k):
            changes.append("%s: recovered missing 5' exon" % gene)
    return annotations, changes


_SIZE_SLACK = 15        # bp either side of the estimated first-exon length


def _recover_first_exon_sized(feat, genome_seq, profile, median_len, est_len):
    """Like `_recover_first_exon` but searching first-exon lengths around `est_len`
    instead of the fixed 3-12 bp list, so a 435 bp rpoC1 exon 1 is reachable."""
    e2s, e2e, strand = feat.start, feat.end, feat.strand
    glen = len(genome_seq)
    lengths = [L for L in range(max(3, est_len - _SIZE_SLACK), est_len + _SIZE_SLACK + 1)]

    exon2_coding = {}
    for d in range(_MAX_TRIM):
        if strand == 1:
            if genome_seq[e2s + d - 2:e2s + d] in _ACCEPTORS:
                exon2_coding[d] = genome_seq[e2s + d:e2e]
        else:
            if _rc(genome_seq[e2e - d:e2e - d + 2]) in _ACCEPTORS:
                exon2_coding[d] = _rc(genome_seq[e2s:e2e - d])
    if not exon2_coding:
        return False

    best = None
    if strand == 1:
        lo = max(0, e2s - _MAX_INTRON - est_len - _SIZE_SLACK)
        hi = e2s - _MIN_INTRON
        for g in range(lo, hi):
            if genome_seq[g:g + 3] != "ATG":
                continue
            for L in lengths:
                e1e = g + L
                if e1e > e2s - _MIN_INTRON:
                    continue
                if genome_seq[e1e:e1e + 2] != _DONOR:
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
        hi = min(glen, e2e + _MAX_INTRON + est_len + _SIZE_SLACK)
        for g in range(lo, hi):
            if _rc(genome_seq[g - 3:g]) != "ATG":
                continue
            for L in lengths:
                e1s = g - L
                if e1s < e2e + _MIN_INTRON:
                    continue
                if _rc(genome_seq[e1s - 2:e1s]) != _DONOR:
                    continue
                e1c = _rc(genome_seq[e1s:g])
                for d, e2c in exon2_coding.items():
                    ok, sup = _score_spliced(e1c + e2c, profile, median_len)
                    if ok and sup >= _MIN_SUPPORT:
                        cand = (sup, -(e1s - (e2e - d)), (e1s, g), d)
                        if best is None or cand > best:
                            best = cand
    if best is None:
        return False
    _sup, _gap, (a, b), d = best
    if strand == 1:
        feat.exons = [(a, b), (e2s + d, e2e)]
    else:
        feat.exons = [(e2s, e2e - d), (a, b)]
    feat.exons.sort()
    feat.start, feat.end = feat.exons[0][0], feat.exons[-1][1]
    feat.has_intron = True
    feat.notes.append("recovered missing 5' exon (%d bp)" % (b - a))
    return True


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
    glen = len(genome_seq)
    # Length around the circle, not e - s: a CDS that crosses the origin has
    # start > end, and the subtraction is a large negative number whose remainder
    # says nothing about the reading frame.
    span = _coords.region_length(s, e, glen)
    if span == 0 or span % 3 != 0:
        return False
    valid_starts = _valid_starts(feat.gene_name)
    # Search window. The fixed pad covers most genes, but a long gene whose
    # Exonerate hit is only a short conserved core needs a window wide enough to
    # still reach the real stop codon — ycf1 typically comes back as ~800 bp of a
    # ~5 kb ORF, and with the fixed 3.5 kb pad the stop lies outside the view, so
    # completion silently fails. Scale the pad to the gene's own coding length.
    pad = max(_ORF_PAD, int(median_len * 3 * 1.8))
    # The window is taken around the circle. It used to be genome_seq[max(0, s-pad):
    # min(glen, e+pad)], which is EMPTY for a CDS crossing the origin — completion
    # then silently failed and the gene was left without its terminal stop codon
    # (accD, 3 bp short, in every genome deposited with an origin inside it).
    win = _coords.open_window(genome_seq, s, e, glen, buffer=pad)
    wlen = len(win.seq)
    cs_w = (s - win.start) % glen                  # the CDS, in window coordinates
    ce_w = cs_w + span
    if strand == 1:
        view = win.seq
        tomap = lambda a, b: win.interval(a, b)
        cs, ce = cs_w, ce_w
    else:
        view = _rc(win.seq)
        tomap = lambda a, b: win.interval(wlen - b, wlen - a)
        cs, ce = wlen - ce_w, wlen - cs_w
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
            # two ATGs equidistant from the alignment edge, one upstream and one
            # downstream, tie here; prefer the upstream one rather than whichever
            # the list happened to hold first
            best_c = min(near, key=lambda a: (abs(a[0] - cs), a[0]))[0] if near else cs

    if (best_c, ce2) == (cs, ce):
        return False                              # already complete — no change
    # The nearest-ATG fallback above can return `cs`, the alignment edge, when no
    # start codon was found nearby. Completing an ORF onto a codon that cannot
    # initiate translation produces a gene model that is wrong at its 5' end, which
    # happened for 4.2% of completed CDS. Refuse instead.
    if view[best_c:best_c + 3].upper() not in valid_starts:
        return False
    core = _translate_orf(view[best_c:ce2], valid_starts)
    core = core[:-1] if core.endswith("*") else core
    if "*" in core:                               # would contain an internal stop
        return False
    if not (0.5 * median_len <= len(core) <= 1.6 * median_len):
        return False
    ns, ne = tomap(best_c, ce2)
    feat.start, feat.end = ns, ne
    # One exon on the circle is one or two ascending arcs; a completed ORF that now
    # crosses the origin must record both, or every later step measures it as empty.
    feat.exons = _coords.region_arcs(ns, ne, glen)
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


# ── frame-broken CDS rescue ───────────────────────────────────────────────────
_FS_WINDOW = 400      # bp searched either side when re-deriving a broken boundary


def rescue_frame_broken_cds(annotations, genome_seq, gene_catalog, protein_db=None):
    """Re-derive the boundary of a CDS whose spliced length is not a whole number
    of codons, by searching for a complete ORF around it.

    Exonerate's protein2genome models frameshifts: when the target carries a 1-2 bp
    indel the aligner keeps going across it, and the reported extent is then not a
    multiple of 3. Such a feature has no readable frame, so its translation picks up
    spurious internal stops — and both repair steps that would normally fix a
    boundary refuse to touch it, because `_complete_orf` gives up as soon as it sees
    an internal stop and the terminal-stop snap requires an intact frame. The broken
    extent is therefore written out verbatim.

    This is the fallback for exactly those features. Following PGA (which resolves a
    non-standard gene by hunting for a valid start and a valid stop rather than
    trusting the alignment), we scan all three frames of a window around the hit for
    complete ORFs — valid start codon, in-frame stop, no internal stop — and take
    the one that best matches the gene's reference proteins.

    Unlike PGA we never delete: PGA drops a gene it cannot resolve (5.4 genes per
    genome in its own warning logs, against the 0.37 broken CDS per genome this
    rescues), which costs real annotations. A feature we cannot re-derive keeps its
    coordinates untouched and is written as a partial CDS instead.

    Gated on `spliced % 3` so it can only ever see features that are already
    invalid; a well-formed CDS is never touched.
    """
    changes = []
    if not protein_db:
        return annotations, changes
    glen  = len(genome_seq)
    cache = {}
    for feat in annotations:
        if feat.gene_type != "CDS" or feat.is_pseudogene or feat.exon_strands:
            continue
        exons = feat.exons or [(feat.start, feat.end)]
        if sum(e - s for s, e in exons) % 3 == 0:
            continue                              # frame intact — leave alone
        gene = feat.gene_name
        # This rescue rewrites the feature as ONE continuous ORF, which would destroy
        # a genuine splice structure. Only single-exon genes are eligible; a
        # multi-exon gene with a broken frame is left for the splice handlers.
        if len(exons) > 1 or (gene_catalog.get(gene, {}) or {}).get("n_exons", 1) > 1:
            continue
        if gene not in cache:
            cache[gene] = _load_ref_profile(protein_db, gene)
        if not cache[gene]:
            continue
        _profile, median_len, refs = cache[gene]
        starts = _valid_starts(gene)

        w0, w1 = max(0, feat.start - _FS_WINDOW), min(glen, feat.end + _FS_WINDOW)
        if feat.strand == 1:
            view  = genome_seq[w0:w1]
            tomap = lambda a, b: (w0 + a, w0 + b)
        else:
            view  = _rc(genome_seq[w0:w1])
            tomap = lambda a, b: (w1 - b, w1 - a)
        vlen = len(view)

        best = None                               # (identity, length, (s, e))
        for frame in range(3):
            i = frame
            while i + 3 <= vlen:
                if view[i:i + 3].upper() in starts:
                    # walk to the first in-frame stop
                    j = i
                    while j + 3 <= vlen and view[j:j + 3].upper() not in _STOPS:
                        j += 3
                    if j + 3 <= vlen:             # a real terminator was reached
                        prot = _translate_orf(view[i:j + 3], starts)
                        core = prot[:-1] if prot.endswith("*") else prot
                        if (core and "*" not in core
                                and 0.5 * median_len <= len(core) <= 1.6 * median_len):
                            idn = max((difflib.SequenceMatcher(None, core, r,
                                                               autojunk=False).ratio()
                                       for r in refs), default=0.0)
                            if idn >= 0.5 and (best is None or idn > best[0]):
                                best = (idn, len(core), tomap(i, j + 3))
                i += 3
        if best is None:
            continue
        ns, ne = best[2]
        # The rescued ORF must describe the SAME locus, not a neighbour that merely
        # fell inside the search window.
        if min(ne, feat.end) - max(ns, feat.start) <= 0:
            continue
        feat.start, feat.end = ns, ne
        feat.exons = [(ns, ne)]
        feat.notes.append(
            "frame-broken alignment re-derived to a complete ORF (id=%.2f)" % best[0])
        changes.append("%s: frame-broken boundary re-derived" % gene)
    return annotations, changes


# ── ycf1 ──────────────────────────────────────────────────────────────────────
# A copy must reach this fraction of the catalog length to count as functional;
# below it the copy is either the genuine IR-junction remnant or an Exonerate core
# that could not be completed.
_YCF1_FULL_RATIO = 0.6
# Ignore anything shorter than this — stray micro-hits, not a reportable remnant.
_YCF1_MIN_FRAG   = 200
# Floor for the IR-derived copy. It can be lower than _YCF1_MIN_FRAG because that
# copy is not a standalone hit: it is anchored on a confirmed full-length ycf1 and
# on the measured IR boundaries. NCBI itself annotates ycf1 remnants as short as
# 81 bp. Keeping a floor at all guards against an IR boundary that is off by a few
# bases producing a spurious sliver.
_YCF1_MIN_DERIVED = 60


def _arcs_to_span(arcs, genome_len):
    """(start, end) read forward around the circle, for arcs that ARE contiguous.

    Two arcs meeting at the origin describe one interval crossing position 1, and
    the span is (the arc that starts later, the arc that ends earlier). Taking
    min(start) and max(end) instead returns (0, genome_len) — the whole genome —
    which is how a derived ycf1 copy of a few hundred bases could be recorded as
    spanning the entire plastome.

    Returns None when the arcs are NOT contiguous. Bridging a gap here would invent
    coordinates covering bases that were never mirrored: [(100,200),(400,500)]
    became (100,500), silently claiming the 200 bp between them.
    """
    arcs = sorted(arcs)
    if not arcs:
        return None
    if len(arcs) == 1:
        return arcs[0]
    if len(arcs) == 2 and arcs[0][0] == 0 and arcs[-1][1] == genome_len:
        return arcs[-1][0], arcs[0][1]          # one interval across the origin
    if all(b == arcs[i + 1][0] for i, (_, b) in enumerate(arcs[:-1])):
        return arcs[0][0], arcs[-1][1]          # abutting: genuinely one interval
    return None                                  # a real gap: not one interval


def _mirror_into_other_ir(span, ir_boundaries, genome_len=None):
    """Reflect a genomic interval from one inverted repeat into the other.

    IRb and IRa are reverse complements of each other, so a base k positions from
    the START of IRb sits k positions from the END of IRa. Returns the mirrored
    (start, end) in the opposite repeat, or None when the interval is not wholly
    inside either repeat. Pure arithmetic on the boundaries `ir_detector` already
    computed — no search, no extra tool call.
    """
    irb, ira = ir_boundaries.get("IRb"), ir_boundaries.get("IRa")
    if not irb or not ira:
        return None
    s, e = span
    # Containment on the circle: a repeat recorded with start > end crosses the
    # origin, and `irb[0] <= s and e <= irb[1]` is then false for every interval
    # inside it.
    if _within(span, irb, genome_len):
        src, dst = irb, ira
    elif _within(span, ira, genome_len):
        src, dst = ira, irb
    else:
        return None
    if genome_len:
        off_s = (s - src[0]) % genome_len
        off_e = (e - src[0]) % genome_len
        return ((dst[1] - off_e) % genome_len, (dst[1] - off_s) % genome_len)
    return (dst[1] - (e - src[0]), dst[1] - (s - src[0]))


def _within(span, region, genome_len):
    """Is the whole interval inside the region, reading forward around the circle?"""
    s, e = span
    if not genome_len:
        return region[0] <= s and e <= region[1]
    want = _coords.region_arcs(s, e, genome_len)
    have = _coords.region_arcs(region[0], region[1], genome_len)
    from plastanno.core.coords import _arcs_bp, _inter_bp
    return _arcs_bp(want) > 0 and _inter_bp(want, have) == _arcs_bp(want)


def handle_ycf1(annotations, genome_seq, ir_boundaries, gene_catalog=None):
    """Resolve the ycf1 copies of a plastome.

    ycf1 straddles the IRb/SSC junction, so a standard plastome carries two copies:
    the full-length functional gene in the SSC, and a truncated duplicate inside
    the IR that NCBI annotates as a `gene` with /pseudo and no CDS.

    ycf1 is the most divergent plastid ORF, so Exonerate routinely reports only a
    short conserved core of it (e.g. 807 bp of a 5040 bp gene) and the HMM engine
    often misses it entirely. `reconcile._select` therefore exempts ycf1 from its
    minimum-length filter — that filter used to delete the gene outright from
    23/111 DEV genomes that have it — and the length decision is made here instead,
    after `complete_single_exon_orfs` has had its chance to extend each core to a
    real ORF:

      * copies at (near) full length stay functional CDS;
      * a short copy inside the IR is the expected truncated duplicate and is kept
        as a pseudogene;
      * if no copy could be completed, the longest remnant is still reported as a
        pseudogene rather than dropped silently;
      * any other short copy is a stray core and is dropped, which is what the
        length filter used to do.
    """
    ycf1_feats = [a for a in annotations if a.gene_name == "ycf1"]
    if not ycf1_feats:
        return annotations, []
    other = [a for a in annotations if a.gene_name != "ycf1"]

    expected = (gene_catalog or {}).get("ycf1", {}).get("expected_len", 0)
    min_full = expected * _YCF1_FULL_RATIO

    def spliced(f):
        return _coords.spliced_length(f)

    # `_select` clusters overlapping copies of a gene, but it runs BEFORE ORF
    # completion. Two disjoint Exonerate cores of the same ycf1 are separate
    # clusters at that point and both survive; completion then extends them into
    # (nearly) the same full-length ORF, leaving a duplicate. Collapse copies that
    # overlap now, keeping the one closest to the catalog length.
    _glen = len(genome_seq or "") or None

    # Group by overlap, then keep the best of each group. This used to sort by
    # (start, end) and walk the list keeping the first non-overlapping candidate —
    # both the sort key and the traversal depend on where position 1 falls, so a
    # chain where A overlaps B and B overlaps C but A does not overlap C resolved
    # differently depending on which end of the chain was reached first. Grouping
    # is transitive and needs no order; overlap is measured on the circle.
    _n = len(ycf1_feats)
    _parent = list(range(_n))

    def _find(x):
        while _parent[x] != x:
            _parent[x] = _parent[_parent[x]]
            x = _parent[x]
        return x

    # Two candidates are the same copy only when each covers most of the other.
    # Grouping on ANY overlap is transitive, so one candidate straddling two real
    # loci would bridge them into a single group and delete a genuine copy — ycf1
    # spans the SSC/IR junction, where a long call really can touch both sides.
    # Reciprocal coverage cannot be bridged that way.
    _DUP_RECIPROCAL = 0.5

    def _same_copy(a, b):
        ov = _coords.region_overlap_bp((a.start, a.end), (b.start, b.end), _glen)
        if ov <= 0:
            return False
        la = _coords.region_length(a.start, a.end, _glen)
        lb = _coords.region_length(b.start, b.end, _glen)
        return (min(la, lb) > 0
                and ov / la >= _DUP_RECIPROCAL and ov / lb >= _DUP_RECIPROCAL)

    for _i in range(_n):
        for _j in range(_i + 1, _n):
            if _same_copy(ycf1_feats[_i], ycf1_feats[_j]):
                ra, rb = _find(_i), _find(_j)
                if ra != rb:
                    _parent[ra] = rb

    _groups = {}
    for _i in range(_n):
        _groups.setdefault(_find(_i), []).append(ycf1_feats[_i])

    def _seq_key(f):
        """Invariant last resort: the feature's own bases, not its coordinates."""
        try:
            return hashlib.sha256(
                _coords.extract(genome_seq, f, _glen).encode()).hexdigest()
        except Exception:
            return ""

    def _rank_evidence(f):
        """Biological evidence only — how close the model is to the expected
        length. The sequence digest is NOT part of this: it orders, it does not
        rank."""
        return -abs(spliced(f) - expected) if expected else spliced(f)

    deduped = []
    for members in _groups.values():
        primary, alts, amb, orderable = _ambiguity.resolve(
            members, evidence=_rank_evidence, digest=_seq_key,
            model=lambda f: (tuple(_coords.occupied_arcs(f, _glen)), f.strand))

        def _rec(a):
            return {"caller": a.engine or "?", "gene_name": a.gene_name,
                    "type": a.gene_type, "strand": a.strand,
                    "start": a.start, "end": a.end,
                    "exons": list(_coords.occupied_arcs(a, _glen)),
                    "wrapped": _coords.is_wrapped(a), "score": a.confidence,
                    "distance_bp": 0, "ambiguity_id": amb}

        if amb and not orderable:
            group = [primary] + list(alts)
            for f in group:
                f.ambiguity_id = amb
                f.flag = "NEEDS_REVIEW"
                f.alternatives.extend(_rec(o) for o in group if o is not f)
                f.notes.append(
                    "%d ycf1 models here are identical in sequence and rank equally "
                    "(%s); all are reported as one ambiguous locus" % (len(group), amb))
            deduped.extend(group)
            continue
        if amb:
            primary.ambiguity_id = amb
            primary.flag = "NEEDS_REVIEW"
            for a in alts:
                a.ambiguity_id = amb
                primary.alternatives.append(_rec(a))
            primary.notes.append(
                "%d ycf1 models here rank equally on the evidence (%s); one is "
                "reported and the rest recorded as alternatives" % (len(alts) + 1, amb))
        deduped.append(primary)
    deduped.sort(key=_seq_key)          # a stable order that no rotation changes
    n_dup = len(ycf1_feats) - len(deduped)
    ycf1_feats = deduped

    def in_ir(f):
        _g = len(genome_seq or "") or None
        for key in ("IRa", "IRb"):
            r = (ir_boundaries or {}).get(key)
            if r and _coords.region_overlap_bp((f.start, f.end), r, _g) > 0:
                return True
        return False

    def as_pseudo(f, why):
        f.is_pseudogene     = True
        f.pseudogene_reason = why
        f.flag              = "MEDIUM"

    full  = [f for f in ycf1_feats if spliced(f) >= min_full]
    short = [f for f in ycf1_feats if spliced(f) <  min_full]

    kept, changes = [], []
    for f in full:
        f.is_pseudogene = False
        kept.append(f)

    remnants = [f for f in short if spliced(f) >= _YCF1_MIN_FRAG]
    for f in remnants:
        if in_ir(f):
            as_pseudo(f, "truncated IR copy")
            kept.append(f)
    if not kept and remnants:
        # No copy reached full length (IR-lacking or heavily reduced plastome).
        # Report the best remnant as a pseudogene instead of annotating nothing.
        best = max(remnants, key=spliced)
        as_pseudo(best, "truncated copy; full-length ORF not recovered")
        kept.append(best)

    # Recover the IR-symmetric copy when the engines found only one. ycf1 straddles
    # the SSC/IR junction, so whatever part of the functional copy lies inside one
    # repeat is duplicated verbatim in the other — that duplicate IS the truncated
    # pseudogene. Deriving it from the IR geometry is exact and free, whereas
    # waiting for Exonerate to report a second hit on a gene this divergent failed
    # on 8 of the 125 DEV genomes whose reference does annotate the pseudogene.
    # The copy is marked as derived so it is never mistaken for an evidence-based
    # call, and it is emitted as a pseudogene (gene only, no CDS), so it cannot
    # affect CDS scoring.
    if (len(kept) == 1 and not kept[0].is_pseudogene
            and ir_boundaries and ir_boundaries.get("IRa") and ir_boundaries.get("IRb")):
        f = kept[0]
        _g = len(genome_seq or "") or None
        for key in ("IRb", "IRa"):
            rng = ir_boundaries[key]
            # The part of the functional copy that lies inside this repeat, taken
            # on the circle. max(start, lo)/min(end, hi) is the linear intersection
            # and is wrong the moment either the gene or the repeat wraps.
            inter = _coords._inter_bp
            f_arcs = _coords.region_arcs(f.start, f.end, _g)
            r_arcs = _coords.region_arcs(rng[0], rng[1], _g)
            shared = [(max(a, c), min(b, d)) for a, b in f_arcs for c, d in r_arcs
                      if min(b, d) > max(a, c)]
            if not shared:
                continue
            if sum(b - a for a, b in shared) < _YCF1_MIN_DERIVED:
                continue
            # Each shared arc is mirrored on its own. Collapsing them first with
            # min(start)/max(end) turns two arcs lying either side of the origin
            # into one interval covering nearly the whole genome, and the mirror of
            # that is meaningless.
            mirrored = []
            for _a, _b in shared:
                _m = _mirror_into_other_ir((_a, _b), ir_boundaries, _g)
                if _m:
                    mirrored.extend(_coords.region_arcs(_m[0], _m[1], _g))
            if not mirrored or inter(sorted(mirrored), f_arcs) > 0:
                continue
            mirrored = sorted(mirrored)
            m = _arcs_to_span(mirrored, _g)
            if m is None:
                # the mirrored pieces do not form one interval; emitting a span
                # over the gap would claim bases that were never mirrored
                continue
            g = Feature(
                gene_name="ycf1", gene_type="CDS",
                product=f.product or "hypothetical chloroplast RF1",
                start=m[0], end=m[1], strand=-f.strand,
                exons=mirrored,
                engine=f.engine, confidence=f.confidence,
                notes=["derived from the %s copy by inverted-repeat symmetry" % key],
            )
            as_pseudo(g, "truncated IR copy")
            kept.append(g)
            changes.append("ycf1: derived truncated %s copy from IR symmetry"
                           % ("IRa" if key == "IRb" else "IRb"))
            break

    n_pseudo  = sum(1 for f in kept if f.is_pseudogene)
    n_dropped = len(ycf1_feats) - len(kept) + n_dup
    if kept:
        changes.append("ycf1: %d functional copy(ies), %d pseudogene(s)"
                       % (len(kept) - n_pseudo, n_pseudo))
    if n_dropped:
        changes.append("ycf1: dropped %d unrecoverable fragment(s)" % n_dropped)

    return other + kept, changes


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

    # 3a-bis. Recover a missing 5' exon (rpoC1/atpF/rps16/clpP annotated as one exon)
    annotations, changes = recover_missing_first_exon(
        annotations, genome_seq, gene_catalog, protein_db
    )
    all_changes.extend(changes)

    # 3b. Complete truncated single-exon ORFs (accD, ndhK, rps18, ndhI, …)
    annotations, changes = complete_single_exon_orfs(
        annotations, genome_seq, gene_catalog, protein_db
    )
    all_changes.extend(changes)

    # 3c. Re-derive CDS whose alignment carried a frameshift (no readable frame)
    annotations, changes = rescue_frame_broken_cds(
        annotations, genome_seq, gene_catalog, protein_db
    )
    all_changes.extend(changes)

    # 4. ycf1 full-length / IR-junction pseudogene
    annotations, changes = handle_ycf1(
        annotations, genome_seq, ir_boundaries, gene_catalog
    )
    all_changes.extend(changes)

    return annotations, all_changes
