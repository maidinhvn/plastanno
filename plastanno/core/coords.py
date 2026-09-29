"""The coordinate contract: one definition of what a feature's positions mean.

A plastome is circular, so a feature may cross position 1. Such a feature has
`start > end`. Three different questions are asked of a feature's coordinates and
they do NOT have the same answer, which is why a single `arcs()` helper was not
enough:

    occupied_arcs     ascending, non-overlapping intervals the feature's exons
                      actually cover. For measuring: length, overlap of the exon
                      model, containment. Introns are NOT included.
                      Order is meaningless here.

    locus_arcs        the same footprint with the introns filled in — where the
                      locus sits, as one piece per side of the origin. For asking
                      "are these two calls the same locus?", where an engine that
                      found the intron and one that did not must still match.
                      Only ever compare locus_arcs of the SAME gene: filling the
                      trnK-UUU intron makes it cover matK, and filling a
                      trans-spliced rps12 makes it cover 69 kb of unrelated loci.

    transcript_parts  the exons in the order they are transcribed, each with its
                      own strand. For extracting sequence and for writing a
                      location. Order is the whole point, and for a trans-spliced
                      gene the strands differ between parts.

Distances near the origin must use `circular_distance`: 2 bp before position 1 and
1 bp after it are 3 bases apart, not 151,759.
"""


def is_wrapped(feat) -> bool:
    """True when the feature crosses the origin."""
    return feat.start > feat.end


def occupied_arcs(feat, genome_len=None) -> list:
    """Ascending genomic intervals the feature's exons occupy. For measurement.

    A wrapped feature with no exons recorded can only be split at the origin when
    the genome length is known; without it the input is unusable and [] is
    returned, which callers must read as "no overlap" rather than guess at.
    """
    if feat.exons:
        ivs = sorted((s, e) for s, e in feat.exons if e > s)
    elif is_wrapped(feat):
        if not genome_len:
            return []
        ivs = [(0, feat.end), (feat.start, genome_len)]
        ivs = sorted(iv for iv in ivs if iv[1] > iv[0])
    else:
        ivs = [(feat.start, feat.end)]
    # Merge anything that touches or overlaps. Without this a feature whose exons
    # overlap is measured twice over the shared bases: its length exceeds the span
    # it covers and its Jaccard with itself comes out above 1.
    merged = []
    for s, e in ivs:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def locus_arcs(feat, genome_len=None) -> list:
    """The feature's footprint with introns filled in — one arc per side of the
    origin. For locus identity only, and only between calls of the same gene.

    Two engines calling the same gene disagree constantly about the exon model:
    one finds the intron, the other calls a single span across it. Comparing their
    exon-level arcs then gives a low overlap and the reconciler treats one locus as
    two. Comparing footprints answers the question actually being asked — is this
    the same place on the genome.

    A trans-spliced gene is exempt: its two halves are genuinely far apart, and
    filling the gap between them would claim the whole interval in between.
    """
    arcs = occupied_arcs(feat, genome_len)
    if not arcs or getattr(feat, "is_trans_spliced", False):
        return arcs
    if is_wrapped(feat):
        if not genome_len:
            return arcs
        # The footprint runs forward from start, through the origin, to end. An
        # earlier version filled each side of the origin separately and left the
        # gap across it open, which made the footprint depend on where position 1
        # happens to fall: rotating a gene so that one of its introns straddled the
        # origin shrank its footprint and dropped its containment from 1.00 to 0.33.
        return sorted(iv for iv in [(0, feat.end), (feat.start, genome_len)]
                      if iv[1] > iv[0])
    return [(arcs[0][0], arcs[-1][1])]


def locus_jaccard(a, b, genome_len=None) -> float:
    """Jaccard of two footprints. Same-gene locus identity; see locus_arcs."""
    return _jaccard_of(locus_arcs(a, genome_len), locus_arcs(b, genome_len))


def locus_overlap_fraction(a, b, genome_len=None) -> float:
    """Overlap of two footprints over the shorter one. Same-gene use only."""
    aa, bb = locus_arcs(a, genome_len), locus_arcs(b, genome_len)
    shorter = min(_arcs_bp(aa), _arcs_bp(bb))
    return _inter_bp(aa, bb) / shorter if shorter > 0 else 0.0


def transcript_parts(feat, genome_len=None) -> list:
    """[(start, end, strand)] in transcription order. For sequence and locations.

    Per-exon strands are honoured when present (a trans-spliced gene has parts on
    both strands). Otherwise exons run low-to-high on the plus strand and
    high-to-low on the minus strand, and a wrapped feature puts the arc that
    precedes the origin first.
    """
    ivs = [(s, e) for s, e in (feat.exons or [(feat.start, feat.end)]) if e > s]
    if not ivs and is_wrapped(feat) and genome_len:
        ivs = sorted(iv for iv in [(0, feat.end), (feat.start, genome_len)]
                     if iv[1] > iv[0])
    strands = list(getattr(feat, "exon_strands", None) or [])
    if strands and len(strands) == len(feat.exons or []):
        return [(s, e, st) for (s, e), st in zip(feat.exons, strands)]
    if is_wrapped(feat):
        # Everything from feat.start to the end of the sequence is the tail (before
        # the origin); everything from 0 up to feat.end is the head (after it). The
        # split used to ask whether an arc began at exactly 0, which only identifies
        # the head when the head is a single exon: a wrapped gene with two exons
        # after the origin had the second of them sorted in with the tail, and the
        # extracted sequence came out with its middle exon at the end.
        tail = sorted((iv for iv in ivs if iv[0] >= feat.start),
                      reverse=(feat.strand == -1))
        head = sorted((iv for iv in ivs if iv[0] < feat.start),
                      reverse=(feat.strand == -1))
        # On the plus strand transcription runs tail -> origin -> head. On the minus
        # strand it runs the other way, so the head is transcribed first: reading
        # tail-first there prepends the wrong bases to the 5' end.
        ordered = (tail + head) if feat.strand != -1 else (head + tail)
    else:
        ordered = sorted(ivs, reverse=(feat.strand == -1))
    return [(s, e, feat.strand) for s, e in ordered]


def biological_exon_runs(arcs, genome_len=None):
    """Storage arcs -> biological exons.

    An exon that straddles the origin is STORED as two arcs meeting exactly at
    the seam -- the same form GenBank's join() uses -- so a feature with two
    biological exons can occupy three arcs. Counting arcs would report three
    exons, and any consumer comparing exon COUNTS would then see a difference
    that does not exist.

    Two arcs are one exon when one ends at the genome end and the other begins at
    0. Nothing else is merged: genuinely separate exons that happen to abut
    elsewhere stay separate.

    Returns a list of lists of arcs, one inner list per biological exon.
    """
    ivs = sorted((s, e) for s, e in arcs if e > s)
    if not ivs or not genome_len:
        return [[iv] for iv in ivs]
    head = [iv for iv in ivs if iv[0] == 0]
    tail = [iv for iv in ivs if iv[1] == genome_len]
    if head and tail and head[0] is not tail[0]:
        joined = [tail[0], head[0]]
        rest = [iv for iv in ivs if iv is not head[0] and iv is not tail[0]]
        return sorted([[x] for x in rest] + [joined], key=lambda g: g[0][0])
    return [[iv] for iv in ivs]


def biological_exon_count(arcs, genome_len=None) -> int:
    """How many exons a feature has, with an origin-split exon counted once."""
    return len(biological_exon_runs(arcs, genome_len))


def _arcs_bp(arcs) -> int:
    return sum(e - s for s, e in arcs)


def _inter_bp(aa, bb) -> int:
    return sum(max(0, min(e1, e2) - max(s1, s2))
               for s1, e1 in aa for s2, e2 in bb)


def _jaccard_of(aa, bb) -> float:
    inter = _inter_bp(aa, bb)
    union = _arcs_bp(aa) + _arcs_bp(bb) - inter
    return inter / union if union > 0 else 0.0


def spliced_length(feat, genome_len=None) -> int:
    """Bases the feature covers, introns excluded."""
    return _arcs_bp(occupied_arcs(feat, genome_len))


def overlap_bp(a, b, genome_len=None) -> int:
    return _inter_bp(occupied_arcs(a, genome_len), occupied_arcs(b, genome_len))


def jaccard(a, b, genome_len=None) -> float:
    """Jaccard of the two exon models. For locus identity use locus_jaccard."""
    return _jaccard_of(occupied_arcs(a, genome_len), occupied_arcs(b, genome_len))


def overlap_fraction(a, b, genome_len=None) -> float:
    inter = overlap_bp(a, b, genome_len)
    shorter = min(spliced_length(a, genome_len), spliced_length(b, genome_len))
    return inter / shorter if shorter > 0 else 0.0


def shift_feature(feat, delta, genome_len):
    """A copy of `feat` with every coordinate moved `delta` bases around the circle.

    The inverse-liftover used to check that annotation does not depend on where the
    submitter happened to place position 1. Rotating the sequence by delta and
    annotating it should give, after shifting the result back by -delta, the same
    annotation; anything that differs is an origin artefact, not biology.

    An exon that lands across the new origin is split in two, and the pieces are
    kept in transcription order, so the sequence the feature extracts is unchanged.
    """
    import copy
    L = genome_len

    def mv(p, at_end=False):
        q = (p + delta) % L
        return L if (q == 0 and at_end) else q

    out = copy.copy(feat)
    out.notes = list(feat.notes)
    out.start = mv(feat.start)
    out.end = mv(feat.end, at_end=True)

    # A shift can cut one exon in two at the origin, and the inverse shift has to
    # put it back. The join is identifiable, not guessed: the position the previous
    # frame's origin maps to is exactly `delta % L`, so two transcript-adjacent
    # pieces are two halves of one exon only when they meet THERE. Merging any two
    # pieces that happen to abut would fuse two genuinely separate exons whose
    # intron is zero-length — unusual, but not ours to erase.
    seam = delta % L
    parts, strands = [], []
    for s, e, st in transcript_parts(feat, L):
        ns, ne = mv(s), mv(e, at_end=True)
        if ns < ne:
            pieces = [(ns, ne)]
        else:
            # the exon now straddles position 1; on the minus strand the piece
            # after the origin is transcribed first
            pieces = [(0, ne), (ns, L)] if st == -1 else [(ns, L), (0, ne)]
            pieces = [pc for pc in pieces if pc[1] > pc[0]]
        parts.extend(pieces)
        strands.extend([st] * len(pieces))

    joined, jstr = [], []
    for iv, st in zip(parts, strands):
        if (joined and st == jstr[-1] and seam not in (0, L)
                and (iv[0] == joined[-1][1] == seam or iv[1] == joined[-1][0] == seam)):
            joined[-1] = (min(joined[-1][0], iv[0]), max(joined[-1][1], iv[1]))
        else:
            joined.append(iv); jstr.append(st)

    if getattr(feat, "exon_strands", None) or len(set(jstr)) > 1:
        out.exons, out.exon_strands = joined, jstr    # transcript order matters
    else:
        out.exons = sorted(joined)                    # order re-derived from strand
        out.exon_strands = []
    return out


def circular_distance(x, y, genome_len) -> int:
    """Shortest distance between two positions on the circle.

    Comparing boundaries linearly makes a 2 bp disagreement across the origin look
    like the width of the genome, which is how a trivial difference was reported as
    151,688 bp.
    """
    if not genome_len:
        return abs(x - y)
    d = abs(x - y) % genome_len
    return min(d, genome_len - d)


# ── regions on the circle ─────────────────────────────────────────────────────
#
# A region — LSC, SSC, IRb, IRa, or any search window — is a (start, end) pair read
# FORWARD around the circle from start to end. When start > end it crosses the
# origin. Every helper below reads it that way, which the callers that took
# genome_seq[start:end] did not: for a wrapping LSC that slice is empty, and the
# genes inside it simply disappeared from the annotation.

def region_arcs(start, end, genome_len) -> list:
    """The one or two ascending intervals a region occupies."""
    if start == end:
        return [(0, genome_len)] if genome_len else []      # the whole circle
    if start < end:
        return [(start, end)]
    return [iv for iv in [(0, end), (start, genome_len)] if iv[1] > iv[0]]


def region_length(start, end, genome_len) -> int:
    return sum(e - s for s, e in region_arcs(start, end, genome_len))


def region_contains(region, pos, genome_len) -> bool:
    """Is this position inside the region, reading forward around the circle?"""
    s, e = region
    return any(a <= pos < b for a, b in region_arcs(s, e, genome_len))


def region_overlap_bp(r1, r2, genome_len) -> int:
    return _inter_bp(region_arcs(r1[0], r1[1], genome_len),
                     region_arcs(r2[0], r2[1], genome_len))


def feature_in_region(feat, region, genome_len, fraction=0.5) -> bool:
    """Does enough of the feature lie in the region to call it a member?

    Measured as a fraction of the feature, not by testing a single coordinate: a
    gene at a junction has bases on both sides, and asking only about its start
    assigns it by an accident of orientation.
    """
    arcs = occupied_arcs(feat, genome_len)
    total = _arcs_bp(arcs)
    if not total:
        return False
    inside = _inter_bp(arcs, region_arcs(region[0], region[1], genome_len))
    return inside / total >= fraction


class Window:
    """A contiguous arc of the circular sequence, handed to a tool as a linear string.

    A search window near the origin used to be clamped to the end of the sequence
    (`genome_seq[max(0, s):min(L, e)]`), so the buffer that was meant to give the
    aligner room simply vanished and genes at the boundary were truncated by
    250-500 bp. A window wraps instead, and carries the mapping needed to put the
    tool's local coordinates back on the genome.
    """
    __slots__ = ("seq", "start", "genome_len")

    def __init__(self, seq, start, genome_len):
        self.seq, self.start, self.genome_len = seq, start, genome_len

    def __len__(self):
        return len(self.seq)

    @property
    def wraps(self) -> bool:
        return self.start + len(self.seq) > self.genome_len

    def to_genomic(self, local: int) -> int:
        """A 0-based offset in the window -> its position on the genome."""
        return (self.start + local) % self.genome_len

    def interval(self, ls, le):
        """A local half-open interval -> a genomic (start, end), wrapped if it must be."""
        gs = self.to_genomic(ls)
        ge = self.to_genomic(le - 1) + 1
        return gs, (ge % self.genome_len if ge != self.genome_len else self.genome_len)

    def arcs(self, ls, le) -> list:
        """A local half-open interval -> its one or two ascending genomic arcs."""
        gs, ge = self.interval(ls, le)
        return region_arcs(gs, ge, self.genome_len)


def open_window(genome_seq, start, end, genome_len=None, buffer=0):
    """Extract [start, end) plus `buffer` bases each side, reading round the origin.

    Never longer than the genome, and never clamped: a window that would run past
    the end continues from the beginning, which is what the sequence actually does.
    """
    L = genome_len or len(genome_seq)
    span = region_length(start, end, L)
    if span == 0:
        span = L
    total = min(L, span + 2 * buffer)
    ws = (start - buffer) % L
    if total >= L:
        ws, total = 0, L
    seq = genome_seq[ws:ws + total] if ws + total <= L \
        else genome_seq[ws:] + genome_seq[:ws + total - L]
    return Window(seq, ws, L)


def geometry_problems(feat, genome_len=None) -> list:
    """Everything structurally wrong with a feature's exon model, as plain strings.

    The measuring helpers here merge overlapping exons so that a length or a Jaccard
    stays meaningful. `extract` does not — it translates the model as given. A
    feature whose exons overlap therefore reports one length to QC and hands a
    longer sequence to the translator, and the two never meet. Merging quietly and
    carrying on hides that; this reports it so the feature can be flagged.
    """
    bad = []
    exons = list(feat.exons or [])
    strands = list(getattr(feat, "exon_strands", None) or [])

    if strands and len(strands) != len(exons):
        bad.append("%d exon strand(s) recorded for %d exon(s)"
                   % (len(strands), len(exons)))

    for s0, e0 in exons:
        if e0 <= s0:
            bad.append("exon %d-%d has no length" % (s0, e0))
        if s0 < 0 or (genome_len and e0 > genome_len):
            bad.append("exon %d-%d lies outside the sequence" % (s0, e0))

    seen = set()
    for iv in exons:
        if iv in seen:
            bad.append("exon %d-%d is listed twice" % iv)
        seen.add(iv)

    ordered = sorted(iv for iv in exons if iv[1] > iv[0])
    for (s1, e1), (s2, e2) in zip(ordered, ordered[1:]):
        if s2 < e1:
            bad.append("exons %d-%d and %d-%d overlap by %d bp: the spliced length "
                       "and the translated sequence cannot both be right"
                       % (s1, e1, s2, e2, e1 - s2))

    if is_wrapped(feat) and not exons and not genome_len:
        bad.append("crosses the origin but records no exons, and the genome length "
                   "is unknown, so it cannot be placed")

    # A non-wrapped feature's exons must lie inside its own span. When they do not,
    # start/end and exons disagree about where the gene is.
    if exons and not is_wrapped(feat):
        lo = min(s0 for s0, _ in exons)
        hi = max(e0 for _, e0 in exons)
        if lo < feat.start or hi > feat.end:
            bad.append("exons span %d-%d but the feature is recorded as %d-%d"
                       % (lo, hi, feat.start, feat.end))
    return bad


def circular_midpoint(feat, genome_len) -> int:
    """The middle of the feature, measured along the circle.

    (start + end) // 2 puts a feature that crosses the origin on the far side of
    the genome from where it actually is, which is how an origin-crossing tRNA was
    tested for inverted-repeat membership against the wrong region.
    """
    if not genome_len or not is_wrapped(feat):
        return (feat.start + feat.end) // 2
    span = (feat.end - feat.start) % genome_len
    return (feat.start + span // 2) % genome_len


def boundary_distance(a, b, genome_len) -> int:
    """Largest circular disagreement between two features' outer boundaries."""
    return max(circular_distance(a.start, b.start, genome_len),
               circular_distance(a.end, b.end, genome_len))


def extract(genome_seq, feat, genome_len=None) -> str:
    """Coding-orientation sequence, origin and per-exon strands handled."""
    from Bio.Seq import Seq
    out = []
    for s, e, st in transcript_parts(feat, genome_len):
        piece = genome_seq[s:e]
        out.append(str(Seq(piece).reverse_complement()) if st == -1 else piece)
    return "".join(out)


def tbl_intervals(feat, genome_len=None) -> list:
    """(start, end) pairs for an NCBI feature table: 1-based, minus strand written
    descending, transcription order preserved."""
    return [((e, s + 1) if st == -1 else (s + 1, e))
            for s, e, st in transcript_parts(feat, genome_len)]
