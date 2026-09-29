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
import copy
import hashlib
from collections import Counter
from typing import List, Dict, Tuple
from .feature import Feature
from .ir_genes import LEGACY_IR_CONFLICT_GENES  # legacy path only;
# the pooled path decides copy number on geometry, never on a gene-name list
from . import coords as _coords
from . import ambiguity as _ambiguity


# ── Matching ──────────────────────────────────────────────────────────────────

def _feature_key(f):
    """An ordering that depends only on what the feature IS, never on the order the
    engine happened to emit it in. Sorting by this before building the assignment
    matrix is what makes the pairing permutation-invariant: with an arbitrary input
    order, two equally good candidates are separated only by whichever scipy's
    solver reaches first, and re-ordering the engine output silently re-annotates
    the genome."""
    return (f.gene_name, f.gene_type, f.strand, f.start, f.end,
            tuple(sorted(f.exons or ())), f.engine or "",
            round(f.s_ref or 0.0, 6), round(f.s_model or 0.0, 6))


def _assign(weights, rows, cols):
    """One-to-one assignment over {(i, j): weight}, on two objectives in order:
    first as MANY pairs as possible, then the highest total weight among those.

    The order matters and getting it wrong loses genes. A single candidate scoring
    1.0 and two candidates scoring 0.5 each have the same total weight, so a plain
    maximum-weight assignment may take the one and leave two features unpaired —
    each then reported separately, as if the engines had never agreed. Offering
    every candidate a bonus larger than any achievable total weight makes an extra
    pair worth more than any amount of weight, and weight decides only between
    assignments of equal size.

    `rows`/`cols` fix the matrix order, so the result is a function of the weights
    alone and not of the order the engines emitted their candidates.
    """
    if not weights:
        return {}
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    ri = sorted({i for i, _ in weights}, key=rows.index)
    cj = sorted({j for _, j in weights}, key=cols.index)
    rpos = {i: k for k, i in enumerate(ri)}
    cpos = {j: k for k, j in enumerate(cj)}
    # Weights are in [0, 1] and at most min(rows, cols) pairs can be made, so this
    # bonus strictly exceeds the total weight of any assignment: one more pair
    # always wins, whatever the weights.
    bonus = min(len(ri), len(cj)) + 1.0
    cost = np.zeros((len(ri), len(cj)))
    for (i, j), w in weights.items():
        cost[rpos[i], cpos[j]] = -(bonus + w)
    out = {}
    for r, c in zip(*linear_sum_assignment(cost)):
        i, j = ri[r], cj[c]
        if (i, j) in weights:
            out[i] = j
    return out


def match_features(
    engine_a: List[Feature],
    engine_b: List[Feature],
    min_overlap: float = 0.5,
    genome_len: int = None,
) -> Dict[str, List]:
    """Bin features into AB / A_only / B_only / conflict.

    Two measurements, asking two different questions.

    *Locus containment* — how much of the smaller footprint the two calls share
    (coords.locus_overlap_fraction) — decides whether a pair describes the same gene
    copy. The engines routinely disagree about the intron: Exonerate splits a gene
    into exons, the HMM path calls one span across it, and either may reach only
    part of the gene. Measured on exon arcs, such a pair scores as a poor match and
    the same copy is reported twice, once per engine.

    The *exon-model Jaccard* is the assignment weight, and the value stored as
    s_overlap. Containment cannot be the weight: a fragment lying wholly inside a
    longer call contains at 1.000, so with several Engine A fragments of one gene
    every candidate ties and the assignment picks between them arbitrarily — on the
    development genome that handed Engine B's rps18 call to a partial instead of the
    exact-match copy. And s_overlap must stay the Jaccard, or a contained fragment
    would score 1.0 for agreement and inflate the confidence of exactly the features
    that deserve it least.

    Pairing runs in two stages, and the order matters. Stage 1 assigns only the
    pairs that clear min_overlap, so a strong pair can never lose its partner to a
    weak one that happens to have a higher exon Jaccard. Stage 2 then pairs what is
    left over, and everything it produces is a conflict. Running one assignment over
    both and classifying afterwards — which is what this did — let a sub-threshold
    candidate take a partner away from a valid one.

    A pair must agree on gene_name, gene_type and strand. The last two were once
    unchecked, so the two inverted-repeat copies of a gene, which lie on opposite
    strands, could merge into one "confirmed by both engines" feature.

    Same locus, but exon models sharing no base at all, is not confirmation. An
    Engine B call that lands inside an Engine A intron contains at 1.000 and agrees
    on nothing. Such a pair is barred from stage 1 ENTIRELY, not merely
    reclassified afterwards: while it was still admitted and only re-binned at the
    end, it could win the assignment and take a partner away from a pair that would
    have been real, and the feature it displaced was then reported as a
    single-engine call. Stage 1 admits a pair only when the containment clears the
    threshold AND the exon models share at least one base.
    """
    bins = {"AB": [], "A_only": [], "B_only": [], "conflict": []}
    if not engine_a or not engine_b:
        return {"AB": [], "A_only": list(engine_a), "B_only": list(engine_b),
                "conflict": []}

    rows = sorted(range(len(engine_a)), key=lambda i: _feature_key(engine_a[i]))
    cols = sorted(range(len(engine_b)), key=lambda j: _feature_key(engine_b[j]))

    pairs = {}                       # (i, j) -> (locus containment, exon Jaccard)
    for i in rows:
        fa = engine_a[i]
        for j in cols:
            fb = engine_b[j]
            if fb.gene_name != fa.gene_name:
                continue
            if fb.gene_type != fa.gene_type or fb.strand != fa.strand:
                continue
            same_locus = _coords.locus_overlap_fraction(fa, fb, genome_len)
            if same_locus > 0:
                pairs[(i, j)] = (same_locus, _coords.jaccard(fa, fb, genome_len))

    # stage 1 — only pairs that are the same locus AND agree on some coding base
    strong = {k: v[1] for k, v in pairs.items()
              if v[0] >= min_overlap and v[1] > 0.0}
    matched = _assign(strong, rows, cols)

    # stage 2 — everything left over, whether it fell short on containment or
    # shares no base with its neighbour. All of it is a conflict; pairing it here
    # only records WHICH other call the engines disagreed with.
    used_a, used_b = set(matched), set(matched.values())
    weak = {k: v[1] for k, v in pairs.items()
            if k not in strong and k[0] not in used_a and k[1] not in used_b}
    for i, j in _assign(weak, rows, cols).items():
        matched[i] = j
        used_a.add(i)
        used_b.add(j)

    paired_b = set(matched.values())
    for i, fa in enumerate(engine_a):
        j = matched.get(i)
        if j is None:
            bins["A_only"].append(fa)
            continue
        same_locus, exon_ovl = pairs[(i, j)]
        if same_locus < min_overlap:
            bins["conflict"].append(
                (fa, engine_b[j],
                 "locus overlap %.2f is below the %.2f threshold"
                 % (same_locus, min_overlap)))
        elif exon_ovl <= 0.0:
            bins["conflict"].append(
                (fa, engine_b[j],
                 "same locus (containment %.2f) but the exon models share no base"
                 % same_locus))
        else:
            bins["AB"].append((fa, engine_b[j], exon_ovl))

    for j, fb in enumerate(engine_b):
        if j not in paired_b:
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

def _internal_stops(seq: str) -> int:
    """In-frame stop codons strictly inside a coding sequence (first and last
    codons excluded). The one definition validate_orf and the rescue gate share."""
    return sum(1 for i in range(3, len(seq) - 3, 3)
               if seq[i:i + 3].upper() in ("TAA", "TAG", "TGA"))


def validate_orf(feat: Feature, genome_seq: str, gene_catalog: dict) -> float:
    """Compute S_orf in [0,1]. Non-CDS features get a structural prior of 0.8."""
    if feat.gene_type != "CDS":
        return 0.8

    # The coding sequence, in transcription order and coding orientation. Slicing
    # genome_seq[start:end] directly returns "" for a feature that crosses the
    # origin, and concatenating exons in ascending order before reverse-
    # complementing gets a wrapped multi-exon gene's parts in the wrong order.
    seq = _coords.extract(genome_seq, feat, len(genome_seq) or None)
    if len(seq) < 30:
        return 0.0

    score = 0.0
    valid_starts = ["ATG"] + SPECIAL_START_CODONS.get(feat.gene_name, [])
    if seq[:3].upper() in valid_starts:
        score += 0.25
    if seq[-3:].upper() in ["TAA", "TAG", "TGA"]:
        score += 0.25
    if _internal_stops(seq) == 0:
        score += 0.25
    exp_len = gene_catalog.get(feat.gene_name, {}).get("expected_len", 0)
    if exp_len > 0:
        if 0.7 <= len(seq) / exp_len <= 1.3:
            score += 0.25
    else:
        score += 0.25
    return score


# ── Confidence-signal availability ──────────────────────────────────────────────

def _from_aragorn(feat: Feature) -> bool:
    """Did ARAGORN call this feature?

    The test used to be `"ARAGORN" in feat.notes`, which is list membership, not a
    substring search: it matched only the bare note "ARAGORN" and missed both other
    forms the caller writes — "ARAGORN intron i(off,len)" for a tRNA with an intron
    and "ARAGORN (crosses the origin)" for one spanning position 1. Those features
    lost the structural prior below and were scored as if nothing had detected them,
    which put an origin-crossing trnH-GUG at confidence 0.40 / NEEDS_REVIEW while
    the same gene in a rotated copy of the same genome scored 0.75 / MEDIUM.
    """
    return any(str(n).startswith("ARAGORN") for n in feat.notes)


def _effective_detection(feat: Feature) -> float:
    """Detection score used as the 'model' signal for non-CDS features."""
    det = max(feat.s_model, feat.s_ref)
    if det <= 0 and _from_aragorn(feat):
        det = 0.7  # de novo structural prediction
    return det


# ── Fragment collapse ──────────────────────────────────────────────────────────

_ENGINE_RANK = {"AB": 3, "AB_conflict": 3, "A": 2, "B": 1, "": 0}

def _quality(f, gene_catalog):
    """Rank a feature within a same-gene cluster."""
    span = _coords.spliced_length(f)  # SPLICED length
    exp = gene_catalog.get(f.gene_name, {}).get("expected_len", 0)
    closeness = (1.0 - min(1.0, abs(span - exp) / exp)) if exp > 0 else 0.0
    # length-closeness dominates (best signal for "complete gene" vs "fragment"),
    # then cross-confirmation, then completeness (span), then confidence.
    return (round(closeness, 3), _ENGINE_RANK.get(f.engine, 0),
            span, round(f.confidence, 3))


_KEY_FLANK = 200


def _invariant_key(f, genome_seq, genome_len, flank=_KEY_FLANK):
    """A last-resort ordering key that a rotation or a strand flip cannot change.

    Coordinates are not usable here — they are precisely what the transform
    changes, so ordering by start or end is deterministic within one presentation
    of a genome and arbitrary between two. The feature's own coding sequence is
    invariant: rotating the genome moves the feature, reverse-complementing it
    moves and flips it, and in both cases the same bases come back.

    The sequence ALONE is not enough to separate two candidates in a repeat, where
    the same bases occur twice. The flanking context is included, in coding
    orientation, because it travels with the feature under both transforms while
    still differing between two positions that merely look alike.
    """
    if not genome_seq:
        return ""
    L = genome_len or len(genome_seq)
    try:
        core = _coords.extract(genome_seq, f, L)
        if flank and L:
            up = _slice_circular(genome_seq, f.start - flank, f.start, L)
            dn = _slice_circular(genome_seq, f.end, f.end + flank, L)
            if f.strand == -1:
                from Bio.Seq import Seq
                up, dn = (str(Seq(dn).reverse_complement()),
                          str(Seq(up).reverse_complement()))
            core = "%s|%s|%s" % (up, core, dn)
        return hashlib.sha256(core.encode()).hexdigest()
    except Exception:
        return ""


def _slice_circular(seq, start, end, genome_len):
    a, b = start % genome_len, end % genome_len
    return seq[a:b] if a < b else seq[a:] + seq[:b]


def _select(features, gene_catalog, ir_boundaries=None, genome_len=None,
            genome_seq=None):
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

    # Overlap and length are computed over the intervals a feature occupies, so a
    # feature crossing the origin — whose start exceeds its end — is measured
    # correctly instead of yielding a negative span.
    def ov_frac(a, b):
        return _coords.overlap_fraction(a, b, genome_len)
    def locus_ov_frac(a, b):
        """Same-gene fragments, compared on the footprint. Two calls of one gene
        that each caught a different exon share no exon arc at all, so on exon arcs
        they never cluster and the gene is emitted twice."""
        return _coords.locus_overlap_fraction(a, b, genome_len)
    def len_ratio(a, b):
        la = _coords.spliced_length(a, genome_len)
        lb = _coords.spliced_length(b, genome_len)
        return min(la, lb) / max(la, lb) if max(la, lb) > 0 else 0.0

    for i in range(n):
        for j in range(i + 1, n):
            a, b = work[i], work[j]
            if a.gene_type != b.gene_type:
                continue
            if a.gene_name == b.gene_name:
                # same gene: footprint, so the intron is not a gap between copies
                if locus_ov_frac(a, b) > 0:
                    union(i, j)
                continue
            # different genes: exon arcs only. Filling trnK-UUU's intron would make
            # it cover matK, and the cross-hit rule would then delete matK.
            of = ov_frac(a, b)
            if of > 0.5 and len_ratio(a, b) > 0.5:
                union(i, j)

    clusters = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(work[i])

    kept = []
    for members in clusters.values():
        # A tie in quality must not be resolved by the order the engines emitted
        # their candidates: that order reverses when the genome is presented on the
        # other strand. Nor by coordinates, which the transform changes. What is
        # left is the sequence itself.
        # Rank by biological evidence ONLY; the digest is a canonical ordering
        # device with no biological meaning and must not decide which model wins.
        primary, alts, amb, orderable = _ambiguity.resolve(
            members,
            evidence=lambda f: _quality(f, gene_catalog),
            digest=lambda f: _invariant_key(f, genome_seq, genome_len),
            model=lambda f: (tuple(_coords.occupied_arcs(f, genome_len)), f.strand))

        def _alt_record(a, ref):
            return {"caller": a.engine or "?", "gene_name": a.gene_name,
                    "type": a.gene_type, "strand": a.strand,
                    "start": a.start, "end": a.end,
                    "exons": list(_coords.occupied_arcs(a, genome_len)),
                    "wrapped": _coords.is_wrapped(a),
                    "score": round(a.confidence or 0.0, 4),
                    "distance_bp": (_coords.boundary_distance(a, ref, genome_len)
                                    if genome_len else 0),
                    "ambiguity_id": amb}

        if amb and not orderable:
            # No semantic primary exists: two distinct models share a sequence and
            # a context, so nothing presentation-independent puts one first.
            # Emitting only one would emit coordinates that depend on the order the
            # engines happened to produce their candidates. Every member is
            # emitted, all with the same ambiguity_id, so this is still ONE review
            # locus however many models it has.
            group = [primary] + list(alts)
            for f in group:
                f.ambiguity_id = amb
                f.flag = "NEEDS_REVIEW"
                f.alternatives.extend(_alt_record(o, f) for o in group if o is not f)
                f.notes.append(
                    "%d models here are identical in sequence and context and rank "
                    "equally (%s); no evidence chooses between them, so all are "
                    "reported as one ambiguous locus" % (len(group), amb))
            kept.extend(group)
            continue

        best = primary
        if amb:
            best.ambiguity_id = amb
            best.flag = "NEEDS_REVIEW"
            for a in alts:
                a.ambiguity_id = amb
                best.alternatives.append(_alt_record(a, best))
            best.notes.append(
                "%d distinct models here rank equally on the evidence (%s); the one "
                "reported is a stable representative chosen by sequence, the others "
                "are recorded as alternatives" % (len(alts) + 1, amb))
        if len(members) > 1:
            best.notes.append("selected over %d overlapping candidate(s)" % (len(members) - 1))
        kept.append(best)

    # Minimum-length quality threshold: drop CDS far shorter than expected.
    # Calibrated on the broad DEV sample — below ratio 0.6, candidates are ~99%
    # false positives (1561 FP vs 18 TP across 123 genomes). Applied to CDS only,
    # and only when the expected length is known; the spliced length is used so
    # multi-exon genes are judged on coding length, not outer span.
    # ycf1 is exempt. It is the most divergent plastid ORF: Exonerate usually
    # reports only a short conserved core of the real ~5 kb gene, and its second
    # copy at the IRb/SSC junction is *genuinely* truncated. Applying the filter
    # here deleted ycf1 entirely from 23 of the 111 DEV genomes that carry it. The
    # length decision is deferred to annotate.special_cases.handle_ycf1, which runs
    # after ORF completion has had a chance to extend the core to the full gene.
    MIN_LEN_RATIO = 0.6
    LEN_FILTER_EXEMPT = {"ycf1"}
    final = []
    for f in kept:
        if f.gene_type == "CDS" and f.gene_name not in LEN_FILTER_EXEMPT:
            exp = gene_catalog.get(f.gene_name, {}).get("expected_len", 0)
            if exp:
                spliced = _coords.spliced_length(f, genome_len)
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


#: what "the model passes the ORF check" means, numerically. validate_orf scores
#: start codon, terminal stop, absence of internal stops, and length-vs-expected
#: at 0.25 each, so 1.0 is clean on all four.
POOLED_ORF_PASS = 1.0


def reconcile_pooled(
    engine_a: List[Feature],
    engine_b: List[Feature],
    genome_seq: str,
    gene_catalog: dict,
    ir_boundaries: dict = None,
    orf_pass: float = POOLED_ORF_PASS,
) -> List[Feature]:
    """Candidate pooling with a B-primary / A-rescue selection rule.

    The H3 ablation found the cross-engine integration in `reconcile()` -- AB
    matching, the coordinate-donor rule, the agreement boost and conflict
    resolution -- bought no detectable improvement in strict CDS coordinate
    concordance over pooling the same candidates through the same selector. What
    DID help was Engine B's coordinates: on the development set Engine B alone
    had the best strict CDS precision (0.8541) of any arm.

    So this path keeps the candidates and drops the integration:

      * both engines found the locus -> keep ENGINE B's coordinates when its
        model passes the ORF check, or is at least as good as Engine A's;
        substitute Engine A only when B does not pass and A is better;
      * one engine found it -> keep it (this is the rescue that lifts inventory);
      * same name, too far apart to be one locus -> for an IR-duplicated gene
        that is the two COPIES and both are kept; otherwise ORF validity decides;
      * agreement between engines earns NO score boost. Two engines agreeing is
        not evidence the coordinates are right, and H3 found the boost bought
        nothing.

    Every feature is therefore scored on single-engine signals only. The score is
    a CANDIDATE SELECTION score: it ranks candidates for this locus and nothing
    more. It is deliberately not called a confidence -- H4 showed the old scalar
    was quantised and stale, and renaming it would relabel the same defective
    quantity.

    NOTE the IR set: this uses the complete IR_DUPLICATED_GENES, not the
    eight-gene subset `reconcile()` is frozen on. See core.ir_genes.
    """
    fallbacks = {}
    out = pooled_candidates(engine_a, engine_b, genome_seq, gene_catalog, orf_pass,
                            fallbacks=fallbacks)

    for f in out:
        _score_pooled(f, genome_seq, gene_catalog)

    survivors = _select_pooled(out, gene_catalog, ir_boundaries, genome_seq)
    return survivors + _rescue(out, survivors, fallbacks, genome_seq,
                               gene_catalog, ir_boundaries)


def _score_pooled(f, genome_seq, gene_catalog):
    """Single-engine scoring, identical for every pooled candidate."""
    f.s_overlap = 0.0                          # no cross-engine agreement term
    f.s_orf = validate_orf(f, genome_seq, gene_catalog)
    if f.gene_type != "CDS":
        f.s_model = _effective_detection(f)
        avail = {"model", "orf"}
    else:
        avail = {"ref", "orf"} if f.s_ref > 0 else {"model", "orf"}
    f.compute_confidence(avail)


def _select_pooled(feats, gene_catalog, ir_boundaries, genome_seq):
    try:
        return _select(feats, gene_catalog, ir_boundaries, len(genome_seq))
    except TypeError:
        return _select(feats, gene_catalog, ir_boundaries)


def _feature_key(f):
    """What makes two emitted features the same: type, name, strand, parts."""
    parts = tuple(tuple(x) for x in f.exons) if f.exons else ((f.start, f.end),)
    return (f.gene_type, f.gene_name, f.strand, parts)


def _rescue(out, survivors, fallbacks, genome_seq, gene_catalog, ir_boundaries):
    """Recover a locus whose preferred candidate did not survive `_select`.

    `pooled_candidates` keeps one candidate of each AB or conflict pair. `_select`
    applies its length filter AFTER the per-cluster pick, so when the kept
    candidate is a fragment it is deleted with nothing left to fall back to, and
    the locus vanishes although a complete call existed. Measured on 30
    development genomes: 6 of 2222 discarded Engine A calls left their gene
    absent from the output; on another 60, ndhA was absent from 21.

    An earlier fix put the discarded candidate into `_select`'s input. On the 14
    development genomes of its sample it rescued 5 loci and MOVED 41, failing on
    every one, because any added candidate
    changes what the selector decides — a new rank key reorders unrelated
    clusters, and union-find can bridge two clusters through the newcomer.

    So nothing is added to `_select`'s input, and its survivors are exactly what
    they were. Only afterwards is a fallback considered, and only if the candidate
    it lost to is gone. Each is then PROBED: `_select` itself is run on deep copies
    of the survivors plus the candidate, and the candidate is accepted only if the
    probe returns every survivor unchanged and keeps the candidate. The judge is
    the real selector — its ycf1 exemption, its rps12 passthrough, its weak-copy
    and IR rules all apply without being re-implemented here — and the copies mean
    a probe cannot touch a real survivor's notes.

    This guarantees the output of reconcile_pooled is a superset of what it was.
    Whether later steps keep that property is a separate, empirical question.
    """
    if not fallbacks:
        return []
    alive = {id(f) for f in survivors}
    cands = [fallbacks[id(k)] for k in out
             if id(k) in fallbacks and id(k) not in alive]
    if not cands:
        return []
    for c in cands:
        _score_pooled(c, genome_seq, gene_catalog)
    # rescue candidates are resolved against each other by the same selector
    cands = _select_pooled(cands, gene_catalog, ir_boundaries, genome_seq)
    cands.sort(key=lambda f: (f.gene_name or "", f.start, f.end))

    accepted = []
    for c in cands:
        base = survivors + accepted
        probe = [copy.deepcopy(x) for x in base] + [copy.deepcopy(c)]
        got = Counter(_feature_key(x) for x in
                      _select_pooled(probe, gene_catalog, ir_boundaries, genome_seq))
        want = Counter(_feature_key(x) for x in base)
        want[_feature_key(c)] += 1
        if got == want:
            c.notes.append(RESCUE_NOTE + " -- the candidate preferred over this "
                           "one did not survive selection, and the locus was "
                           "otherwise empty")
            accepted.append(c)
    return accepted


RESCUE_NOTE = "pooled: rescued"


def revoke_implausible_rescues(annotations, genome_seq):
    """Drop a rescued call that the final QC finds implausible as a gene.

    Returns (kept, revoked). Only calls `_rescue` added are candidates; every
    other feature passes untouched, whatever its flag, so the output stays a
    superset of the unrescued one.

    Why here, after finalize_qc, and not at rescue time: ungated, 4 of 13 rescues
    on 44 dev genomes were calls the reference does not annotate as CDS (three
    infA remnants) and they cancelled the gain. A gate at rescue time removed
    them, but it also removed three of the six rescues that matched the reference
    -- two ndhA with exact outer boundaries and an exact ndhF -- whose raw
    candidates are frame-broken until splice refinement and ORF re-derivation
    repair them (parked/gated_rescue/). The final QC sees the repaired call:
    there, NEEDS_REVIEW or an in-frame internal stop separated exactly the three
    infA remnants from the rest.
    """
    L = len(genome_seq) or None
    kept, revoked = [], []
    for a in annotations:
        if (a.gene_type == "CDS" and any(n.startswith(RESCUE_NOTE) for n in a.notes)
                and (a.flag == "NEEDS_REVIEW"
                     or _internal_stops(_coords.extract(genome_seq, a, L)) > 0)):
            revoked.append(a)
        else:
            kept.append(a)
    return kept, revoked


def pooled_candidates(engine_a, engine_b, genome_seq, gene_catalog,
                      orf_pass=POOLED_ORF_PASS, fallbacks=None):
    """The B-primary / A-rescue selection rule, BEFORE scoring and _select.

    Split out so the rule can be tested on its own. Downstream, `_select`
    collapses same-name overlapping features, so an end-to-end test cannot
    distinguish "the conflict branch kept both candidates" from "it kept one" --
    every conflict pair overlaps by construction, and the selector reduces them
    to one either way. What the branch actually decides is WHICH candidate is
    available to the selector, not how many features are emitted.
    """
    bins = match_features(engine_a, engine_b)
    out = []

    for fa, fb, _exon_ovl in bins["AB"]:
        oa = validate_orf(fa, genome_seq, gene_catalog)
        ob = validate_orf(fb, genome_seq, gene_catalog)
        keep = fb if (ob >= orf_pass or ob >= oa) else fa
        keep.notes.append("pooled: engine %s model kept (orf A=%.2f B=%.2f)"
                          % ("B" if keep is fb else "A", oa, ob))
        out.append(keep)
        if fallbacks is not None:
            fallbacks[id(keep)] = fa if keep is fb else fb

    for fa in bins["A_only"]:
        fa.notes.append("pooled: reference-only rescue")
        out.append(fa)
    for fb in bins["B_only"]:
        fb.notes.append("pooled: model-only")
        out.append(fb)

    # A conflict is two predictions for the SAME PHYSICAL LOCUS -- match_features
    # only pairs calls that overlap, so by construction these are not two copies.
    # B-primary therefore applies here exactly as it does to an AB pair: Engine B
    # wins unless its model fails the ORF check and Engine A's is better.
    #
    # There is deliberately NO gene-name special case. An earlier version kept
    # both candidates when the gene appeared in an "IR genes" list, but that list
    # is a POLICY, not biology: IR content varies between species with repeat
    # expansion and contraction, so a name cannot establish copy number. Real IR
    # copies sit at two different positions, never overlap, and so never reach
    # this bin at all -- they arrive as A_only and B_only and are kept
    # independently on their geometry.
    for fa, fb, reason in bins["conflict"]:
        oa = validate_orf(fa, genome_seq, gene_catalog)
        ob = validate_orf(fb, genome_seq, gene_catalog)
        keep = fb if (ob >= orf_pass or ob >= oa) else fa
        keep.notes.append("pooled: same-locus conflict, engine %s kept (%s)"
                          % ("B" if keep is fb else "A", reason))
        out.append(keep)
        if fallbacks is not None:
            fallbacks[id(keep)] = fa if keep is fb else fb

    return out


def reconcile(
    engine_a: List[Feature],
    engine_b: List[Feature],
    genome_seq: str,
    gene_catalog: dict,
    ir_boundaries: dict = None,   # enables the IR-aware spurious-copy guard in _select
    weights: Tuple = None,   # kept for signature compatibility; unused
    genome_len: int = None,  # required to measure a feature that crosses the origin
) -> List[Feature]:
    """Merge Engine A and B features with provenance and a calibrated confidence."""
    if genome_len is None:
        genome_len = len(genome_seq) or None
    bins = match_features(engine_a, engine_b, genome_len=genome_len)
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
            la = _coords.spliced_length(fa, genome_len)
            lb = _coords.spliced_length(fb, genome_len)
            exp = gene_catalog.get(fa.gene_name, {}).get("expected_len", 0)
            # Whichever engine's length is closer to the expected one donates.
            # This was an absolute floor -- `la / exp < 0.6 and lb > la` -- which
            # only rescued a catastrophically short Engine A call. atpA in
            # NC_039155.1 failed at 0.673: A had 1026 bp of a 1524 bp gene, B had
            # 1521, and A donated anyway. The merged feature then lost a _select
            # paralog cluster to a spurious atpB (1359 bp against atpB's expected
            # 1479 is closer than 1026 against 1524), so the output carried the
            # wrong gene at NEEDS_REVIEW where Engine B alone gave atpA at HIGH.
            # See benchmark_v3/cds/ATPA_TRACE.md.
            if exp and abs(lb - exp) < abs(la - exp):
                donor = fb
        feat_frameshifts = getattr(fa, "frameshifts", 0)
        feat = Feature(
            gene_name=fa.gene_name, gene_type=fa.gene_type,
            product=donor.product or fa.product,
            start=donor.start, end=donor.end, strand=donor.strand,
            exons=donor.exons, protein=donor.protein, engine="AB",
            s_overlap=overlap, s_ref=fa.s_ref, s_model=fb.s_model,
        )
        feat.frameshifts = feat_frameshifts
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
    # FROZEN eight-gene set -- see core.ir_genes. It is a strict subset of the
    # genes the IR actually duplicates, so a conflict on rrn16/trnN-GUU and six
    # others halves the pair here. Not widened: legacy output must stay
    # reproducible against the frozen H4-30 run. The pooled path uses the
    # complete set.
    ir_genes = LEGACY_IR_CONFLICT_GENES
    for fa, fb, reason in bins["conflict"]:
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
                f.notes.append("engines disagreed: %s" % reason)
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
            winner.notes.append("conflict resolved by ORF validity: %s" % reason)
            # Two engines that put the gene in the same place but share no coding
            # base are not confirming each other, and the winner was picked on ORF
            # validity alone. That is a reviewer's call, not a confident annotation.
            if "share no base" in reason:
                winner.flag = "NEEDS_REVIEW"
            results.append(winner)

    return _select(results, gene_catalog, ir_boundaries, genome_len, genome_seq)
