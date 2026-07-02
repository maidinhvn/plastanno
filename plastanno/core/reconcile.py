"""
Reconciliation layer — heart of the hybrid pipeline.

Takes features from Engine A (reference) and Engine B (model) and produces a
unified, scored, provenanced annotation set.

Scoring rewrite: confidence is a weighted average over the signals that are
ACTUALLY AVAILABLE for each feature, renormalised to sum to 1. The old code used
fixed weights and zeroed the missing signals, which capped any single-engine
feature at 0.4, and additionally hard-coded NEEDS_REVIEW for every B-only and
many A-only features. Both issues are removed here.

Fragment collapse: after binning/scoring, same-name features that overlap on the
same strand are collapsed to a single best representative, so a gene detected as
several overlapping fragments is reported once. Spatially separate copies (the
two IR copies of a gene) do not overlap and are preserved.
"""
from typing import List, Dict, Tuple
from .feature import Feature


# ── Matching ──────────────────────────────────────────────────────────────────

def jaccard_overlap(a: Feature, b: Feature) -> float:
    """Compute Jaccard overlap on genomic coordinates."""
    intersection = max(0, min(a.end, b.end) - max(a.start, b.start))
    union        = max(a.end, b.end) - min(a.start, b.start)
    return intersection / union if union > 0 else 0.0


def match_features(
    engine_a: List[Feature],
    engine_b: List[Feature],
    min_overlap: float = 0.5,
) -> Dict[str, List]:
    """Bin features into AB / A_only / B_only / conflict."""
    bins   = {"AB": [], "A_only": [], "B_only": [], "conflict": []}
    used_b = set()

    for fa in engine_a:
        best_match, best_overlap = None, 0.0
        for i, fb in enumerate(engine_b):
            if i in used_b or fb.gene_name != fa.gene_name:
                continue
            ovl = jaccard_overlap(fa, fb)
            if ovl > best_overlap:
                best_overlap, best_match = ovl, (i, fb)

        if best_match is None:
            bins["A_only"].append(fa)
        elif best_overlap >= min_overlap:
            used_b.add(best_match[0])
            bins["AB"].append((fa, best_match[1], best_overlap))
        else:
            used_b.add(best_match[0])
            bins["conflict"].append((fa, best_match[1]))

    for i, fb in enumerate(engine_b):
        if i not in used_b:
            bins["B_only"].append(fb)
    return bins


# ── ORF validation ────────────────────────────────────────────────────────────

# Plastid genes with community-recognised non-ATG (alternative) start codons.
# These are real, documented starts (often created/edited at the RNA level), so a
# genomic ORF beginning with the listed codon is treated as a valid M-start rather
# than rejected/truncated. Frequencies measured on the reference set (n≈200):
# rps19 GTG 60%, ndhD ACG ~40%, psbL ACG 24%, rpl2 ACG 20%, ycf1 GTG/ACG 14%,
# rps12 GTG (trans-spliced 5' exon). Kept gene-specific so non-listed genes stay
# strict ATG-only — a blanket alt-start rule would invent spurious ORFs.
SPECIAL_START_CODONS = {
    "psbL":  ["ACG"],
    "ndhD":  ["ACG", "GTG"],
    "rps19": ["GTG"],
    "rps12": ["GTG"],
    "rpl2":  ["ACG"],
    "ycf1":  ["GTG", "ACG"],
}

def validate_orf(feat: Feature, genome_seq: str, gene_catalog: dict) -> float:
    """Compute S_orf in [0,1]. Non-CDS features get a structural prior of 0.8."""
    if feat.gene_type != "CDS":
        return 0.8

    if feat.exons:                       # spliced coding sequence (multi-exon genes)
        seq = "".join(genome_seq[s:e] for s, e in sorted(feat.exons))
    else:
        seq = genome_seq[feat.start:feat.end]
    if feat.strand == -1:
        from Bio.Seq import Seq
        seq = str(Seq(seq).reverse_complement())
    if len(seq) < 30:
        return 0.0

    score = 0.0
    valid_starts = ["ATG"] + SPECIAL_START_CODONS.get(feat.gene_name, [])
    if seq[:3].upper() in valid_starts:
        score += 0.25
    if seq[-3:].upper() in ["TAA", "TAG", "TGA"]:
        score += 0.25
    internal = sum(1 for i in range(3, len(seq) - 3, 3)
                   if seq[i:i + 3].upper() in ["TAA", "TAG", "TGA"])
    if internal == 0:
        score += 0.25
    exp_len = gene_catalog.get(feat.gene_name, {}).get("expected_len", 0)
    if exp_len > 0:
        if 0.7 <= len(seq) / exp_len <= 1.3:
            score += 0.25
    else:
        score += 0.25
    return score


# ── Confidence-signal availability ──────────────────────────────────────────────

def _effective_detection(feat: Feature) -> float:
    """Detection score used as the 'model' signal for non-CDS features."""
    det = max(feat.s_model, feat.s_ref)
    if det <= 0 and "ARAGORN" in feat.notes:
        det = 0.7  # de novo structural prediction
    return det


# ── Fragment collapse ──────────────────────────────────────────────────────────

_ENGINE_RANK = {"AB": 3, "AB_conflict": 3, "A": 2, "B": 1, "": 0}

def _quality(f, gene_catalog):
    """Rank a feature within a same-gene cluster."""
    span = sum(e - s for s, e in f.exons) if f.exons else (f.end - f.start)  # SPLICED length
    exp = gene_catalog.get(f.gene_name, {}).get("expected_len", 0)
    closeness = (1.0 - min(1.0, abs(span - exp) / exp)) if exp > 0 else 0.0
    # length-closeness dominates (best signal for "complete gene" vs "fragment"),
    # then cross-confirmation, then completeness (span), then confidence.
    return (round(closeness, 3), _ENGINE_RANK.get(f.engine, 0), span, round(f.confidence, 3))


def _select(features, gene_catalog, ir_boundaries=None):
    """
    Unified locus selection — replaces fragment-collapse and HMM-crosshit removal.

    Cluster features that describe the same locus, then keep the single best one
    (by _quality) from each cluster. Two features join the same cluster only when
    they share a gene type AND either:
      * same gene name and they overlap at all          -> duplicate fragments, or
      * different names but one nearly covers the other at a similar size
        (reciprocal overlap > 0.5 and length ratio > 0.5) -> a paralog cross-hit
        (e.g. a psbD profile landing on psbA).

    Different gene types never join, so matK inside the trnK intron, or a CDS that
    merely abuts a tRNA, are both kept. Spatially separate copies (the two IR
    copies of a gene) do not overlap and both survive. A long spurious ORF that
    covers a much shorter real gene is not joined to it (length-ratio guard), so
    the real gene is never removed. rps12 (trans-spliced) is passed through.
    """
    SKIP = {"rps12"}
    passthrough = [f for f in features if f.gene_name in SKIP]
    work = [f for f in features if f.gene_name not in SKIP]
    n = len(work)

    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    def ov_frac(a, b):
        inter = max(0, min(a.end, b.end) - max(a.start, b.start))
        shorter = min(a.end - a.start, b.end - b.start)
        return inter / shorter if shorter > 0 else 0.0
    def len_ratio(a, b):
        la, lb = a.end - a.start, b.end - b.start
        return min(la, lb) / max(la, lb) if max(la, lb) > 0 else 0.0

    for i in range(n):
        for j in range(i + 1, n):
            a, b = work[i], work[j]
            if a.gene_type != b.gene_type:
                continue
            of = ov_frac(a, b)
            if of <= 0:
                continue
            if a.gene_name == b.gene_name:
                union(i, j)
            elif of > 0.5 and len_ratio(a, b) > 0.5:
                union(i, j)

    clusters = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(work[i])

    kept = []
    for members in clusters.values():
        best = max(members, key=lambda f: _quality(f, gene_catalog))
        if len(members) > 1:
            best.notes.append("selected over %d overlapping candidate(s)" % (len(members) - 1))
        kept.append(best)

    # Minimum-length quality threshold: drop CDS far shorter than expected.
    # Calibrated on the broad DEV sample — below ratio 0.6, candidates are ~99%
    # false positives (1561 FP vs 18 TP across 123 genomes). Applied to CDS only,
    # and only when the expected length is known; the spliced length is used so
    # multi-exon genes are judged on coding length, not outer span.
    MIN_LEN_RATIO = 0.6
    final = []
    for f in kept:
        if f.gene_type == "CDS":
            exp = gene_catalog.get(f.gene_name, {}).get("expected_len", 0)
            if exp:
                spliced = sum(e - s for s, e in f.exons) if f.exons else (f.end - f.start)
                if spliced / exp < MIN_LEN_RATIO:
                    continue
        final.append(f)

    # NOTE: a single-copy guard that collapsed same-name LSC/SSC CDS to one copy
    # was tried and REVERTED — it cost 128 true positives on DEV (CDS Sn
    # 91.7→90.5) because the catalog "region" reflects the typical position only:
    # genes near the LSC/IR or SSC/IR junction (rps19, rps15, ndhD, rpl22, …) are
    # legitimately duplicated when the inverted repeat expands, so collapsing by
    # catalog region deletes real second copies. Catalog region ≠ always single
    # copy; do not reintroduce without per-genome IR-boundary awareness.
    #
    # Confidence-gated spurious-copy drop (NOT region-based, so the reverted guard's
    # failure mode does not apply): when a gene already has a CONFIDENT copy
    # (HIGH/MEDIUM, i.e. flag != NEEDS_REVIEW), drop any *additional* copy that is
    # both single-engine AND NEEDS_REVIEW. These are stray short exonerate/HMM hits
    # that cleared the length filter (e.g. a 72 bp psbM fragment 15 kb from the real
    # gene). A real IR-expanded second copy sits in a perfect repeat and is
    # confirmed by both engines (AB / HIGH), so it is never single-engine
    # NEEDS_REVIEW and is protected; a gene whose only copy is weak keeps it.
    confident = {f.gene_name for f in final
                 if f.gene_type == "CDS" and f.flag != "NEEDS_REVIEW"}

    # IR-aware second discriminator (only when the genome's IR boundaries are
    # known): also drop a single-engine copy that is WEAKER than a dual-engine (AB)
    # copy of the same gene AND lies OUTSIDE the inverted repeat. A genuine
    # IR-expanded duplicate sits inside the IR, so it is protected; a spurious
    # exonerate/HMM hit elsewhere (e.g. a MEDIUM-confidence engine-A psbM copy in
    # the LSC, which the NEEDS_REVIEW gate above does not catch) is removed. Gated
    # on ir_boundaries so behaviour is unchanged when IR is unknown; the earlier
    # NON-IR-aware form of this rule cost 41 real IR-junction TP on DEV, which the
    # in-IR guard recovers.
    ab_best = {}
    for f in final:
        if f.gene_type == "CDS" and f.engine == "AB":
            ab_best[f.gene_name] = max(ab_best.get(f.gene_name, 0.0), f.confidence)

    def _in_ir(f):
        if not ir_boundaries:
            return False
        for key in ("IRa", "IRb"):
            r = ir_boundaries.get(key)
            if r and min(f.end, r[1]) - max(f.start, r[0]) > 0:
                return True
        return False

    def _spurious_extra(f):
        if f.gene_type != "CDS" or f.engine not in ("A", "B"):
            return False
        # (a) stray low-confidence copy when a confident copy exists
        if f.flag == "NEEDS_REVIEW" and f.gene_name in confident:
            return True
        # (b) single-engine copy weaker than an AB copy, outside the IR
        if (ir_boundaries and f.gene_name in ab_best
                and f.confidence < ab_best[f.gene_name] and not _in_ir(f)):
            return True
        return False

    final = [f for f in final if not _spurious_extra(f)]
    return passthrough + final


def reconcile(
    engine_a: List[Feature],
    engine_b: List[Feature],
    genome_seq: str,
    gene_catalog: dict,
    ir_boundaries: dict = None,   # enables the IR-aware spurious-copy guard in _select
    weights: Tuple = None,   # kept for signature compatibility; unused
) -> List[Feature]:
    """Merge Engine A and B features with provenance and a calibrated confidence."""
    bins = match_features(engine_a, engine_b)
    results = []

    # AB: both engines found the gene
    for fa, fb, overlap in bins["AB"]:
        # Coordinate donor is normally Engine A (exonerate gives precise CDS
        # boundaries). Exception: when A is only a FRAGMENT of a CDS that B
        # captured in full, trusting A's span makes the merged feature fail the
        # _select length filter and the whole copy is lost (this is what dropped
        # ycf2's ~6.8 kb IRa copy — exonerate breaks the large gene into pieces
        # while the HMM finds it whole). In that case take B's fuller coordinates.
        donor = fa
        if fa.gene_type == "CDS":
            exp = gene_catalog.get(fa.gene_name, {}).get("expected_len", 0)
            if exp:
                la = sum(e - s for s, e in fa.exons) if fa.exons else (fa.end - fa.start)
                lb = sum(e - s for s, e in fb.exons) if fb.exons else (fb.end - fb.start)
                if la / exp < 0.6 and lb > la:
                    donor = fb
        feat = Feature(
            gene_name=fa.gene_name, gene_type=fa.gene_type,
            product=donor.product or fa.product,
            start=donor.start, end=donor.end, strand=donor.strand,
            exons=donor.exons, protein=donor.protein, engine="AB",
            s_overlap=overlap, s_ref=fa.s_ref, s_model=fb.s_model,
        )
        feat.s_orf = validate_orf(feat, genome_seq, gene_catalog)
        if feat.gene_type == "CDS":
            avail = {"overlap", "ref", "model", "orf"}
        else:
            feat.s_model = max(feat.s_model, _effective_detection(fb))
            avail = {"overlap", "model", "orf"}
        feat.compute_confidence(avail)
        feat.notes.append("confirmed by both engines")
        results.append(feat)

    # A_only: reference transfer only
    for fa in bins["A_only"]:
        fa.engine = "A"
        fa.s_overlap = 0.0
        fa.s_orf = validate_orf(fa, genome_seq, gene_catalog)
        avail = {"ref", "orf"}
        fa.compute_confidence(avail)
        fa.notes.append("reference-only (model did not confirm)")
        results.append(fa)

    # B_only: model only (HMM CDS, or tRNA/rRNA)
    for fb in bins["B_only"]:
        fb.engine = "B"
        fb.s_overlap = 0.0
        fb.s_orf = validate_orf(fb, genome_seq, gene_catalog)
        if fb.gene_type != "CDS":
            fb.s_model = _effective_detection(fb)
        avail = {"model", "orf"}
        fb.compute_confidence(avail)
        fb.notes.append("model-only (no reference transfer)")
        results.append(fb)

    # conflict: same gene name, different position
    ir_genes = {"rpl2", "rpl23", "ndhB", "rps7", "ycf2",
                "trnI-CAU", "trnI-GAU", "trnA-UGC"}
    for fa, fb in bins["conflict"]:
        if fa.gene_name in ir_genes:
            for f in (fa, fb):
                f.engine = "AB"
                f.s_orf = validate_orf(f, genome_seq, gene_catalog)
                if f.gene_type == "CDS":
                    avail = {"ref", "orf"} if f.s_ref > 0 else {"model", "orf"}
                else:
                    f.s_model = _effective_detection(f)
                    avail = {"model", "orf"}
                f.compute_confidence(avail)
                f.notes.append("IR copy")
                results.append(f)
        else:
            fa.s_orf = validate_orf(fa, genome_seq, gene_catalog)
            fb.s_orf = validate_orf(fb, genome_seq, gene_catalog)
            winner = fa if fa.s_orf >= fb.s_orf else fb
            winner.engine = "A" if winner is fa else "B"
            winner.s_overlap = 0.0
            if winner.gene_type == "CDS":
                avail = {"ref", "orf"} if winner.s_ref > 0 else {"model", "orf"}
            else:
                winner.s_model = _effective_detection(winner)
                avail = {"model", "orf"}
            winner.compute_confidence(avail)
            winner.notes.append("conflict resolved by ORF validity")
            results.append(winner)

    return _select(results, gene_catalog, ir_boundaries)
