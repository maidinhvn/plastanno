"""Boundary refinement for intron-bearing tRNA, by glocal exon alignment.

Why this exists
---------------
Intron-bearing tRNA are the tool's worst class by a wide margin: 18.7% of them
have exactly right coordinates against 86.6% for intron-free ones, and that one
class accounts for the whole tRNA deficit against PGA. Detection is not the
problem -- ARAGORN has a call at essentially every such locus -- the boundaries
are.

Two routes were measured and closed before this one:

  benchmark_v3/intron_twopass/   tRNAscan-SE does not find these introns at all.
                                 What it calls an intron is 15-40 bp; a plastid
                                 group II intron is 700-950 bp.
  benchmark_v3/groupII_probe/    No Rfam group II model spans one. They are
                                 domain models of CLEN 77-225, and 0 of 70 test
                                 introns had both ends placed by a CM hit.

This route does not model the intron at all, and that is the point. The exons
are the conserved object; the intron is simply the gap between them. Only seven
plastid tRNA carry these introns, and `database/exon_db/` already holds ~50,000
well-formed reference exon pairs for exactly those seven.

Why not BLAST, which is already wired up
----------------------------------------
BLAST is a *local* aligner. An HSP ends where similarity tails off, not where
the exon ends -- which is precisely the error being measured. Here each
reference exon is instead required to align end to end (global on the query,
free ends on the genome), so the boundary is set by reference geometry.

Why the donors vote
-------------------
A single donor imposes its own exon lengths on the target, and raw alignment
score actively *prefers* the longer donor because it collects more matches. On
Arabidopsis: best single donor by raw score 0/8 loci exact, by length-normalised
score 2/8, by vote 4/8. Donors therefore vote per coordinate over the
better-scoring half of the placements.

Measured, criteria locked first, in `benchmark_v3/intron_glocal/`: 18.7% ->
61.8% exact over 395 loci in 60 development genomes, with outer-exact rising
71.4% -> 91.9% so junctions are not being fixed by moving the locus. Removing
every donor from the test genome's own genus gives 62.5%, so this is
generalisation and not a near-self hit.
"""
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Align

# The complete set. Plastid group II introns in tRNA occur in these genes and
# no others, which is what makes this a closed catalogue rather than a de novo
# prediction problem.
INTRON_TRNA_GENES = frozenset((
    "trnA-UGC", "trnI-GAU", "trnL-UAA", "trnK-UUU",
    "trnV-UAC", "trnG-UCC", "trnG-GCC",
))

# Reference files sometimes carry the bare name. Anticodon-less forms map to the
# only intron-bearing member of their family.
_BARE = {"trnA": "trnA-UGC", "trnI": "trnI-GAU", "trnL": "trnL-UAA",
         "trnK": "trnK-UUU", "trnV": "trnV-UAC"}

PAD           = 200     # window grown around the existing call, each side
N_DONORS      = 40      # donors sampled per locus, spread across the DB
MIN_INTRON    = 250     # bp; the shortest plastid tRNA group II intron seen
MAX_INTRON    = 3200    # bp; trnK-UUU's is the longest, ~2.5 kb
MIN_SPLICED   = 60      # a mature tRNA is 70-80 nt
MAX_SPLICED   = 100
VOTE_FRACTION = 0.5     # share of best-scoring placements that get a vote

_RC = str.maketrans("ACGTacgtN", "TGCAtgcaN")

_ALN = Align.PairwiseAligner()
_ALN.mode = "global"
_ALN.match_score = 2.0
_ALN.mismatch_score = -1.0
_ALN.open_gap_score = -4.0
_ALN.extend_gap_score = -1.0
# Glocal: the query (a reference exon) must be consumed end to end, while the
# target (the genome window) may overhang at both ends for free.
#
# In Biopython a target overhang is an END GAP IN THE QUERY, so it is
# query_end_gap_score that must be free. With the two swapped, a placement that
# is exactly right scores -22 instead of +20 and comes back as three scattered
# blocks -- which is how this was first written, and it reported 0 of 8.
#
# Biopython renamed both attributes (query_end_gap_score -> end_deletion_score,
# target_end_gap_score -> end_insertion_score) and warns on the old spelling.
# Set the new names where they exist so a correct configuration does not print
# a deprecation warning on every run.
try:
    _ALN.end_deletion_score = 0.0
    _ALN.end_insertion_score = -1000.0
except AttributeError:                            # Biopython < 1.84
    _ALN.query_end_gap_score = 0.0
    _ALN.target_end_gap_score = -1000.0

_DONOR_CACHE = {}


def _rc(s):
    return s.translate(_RC)[::-1]


def canonical_gene(name):
    """Map a tRNA name onto its intron-bearing form, or None if it has none.

    `trnG` is deliberately absent from `_BARE`: both trnG-UCC and trnG-GCC
    exist and only one of them carries the intron, so a bare `trnG` is not
    resolvable and is left alone.
    """
    if not name:
        return None
    n = str(name).strip()
    if n in INTRON_TRNA_GENES:
        return n
    if "-" not in n:
        return _BARE.get(n)
    return None


def load_donors(exon_db):
    """Well-formed reference exon pairs from the exon DB, grouped by gene.

    `exon_db` is the BLAST db prefix used elsewhere in Engine B; the FASTA sits
    beside it. Parsed by hand rather than through Biopython because this is
    137k short records read once per process, and cached for the same reason.

    Entries are filtered to plausible geometry. The DB also contains a small
    number of malformed over-long templates -- see the note in
    `blast_intron_trna`, where BLAST's default `-max_target_seqs` makes the tool
    depend on them. Nothing here depends on them, and they are dropped.
    """
    key = str(exon_db)
    if key in _DONOR_CACHE:
        return _DONOR_CACHE[key]

    fa = Path(str(exon_db))
    if fa.suffix != ".fasta":
        fa = fa.with_suffix(".fasta")
    pairs = defaultdict(dict)
    if fa.exists():
        sid, buf = None, []

        def _flush():
            if sid is None:
                return
            try:
                head, ex = sid.rsplit("_", 1)      # NC_084135_trnK-UUU , e1
                acc, gene = head.rsplit("_", 1)    # NC_084135 , trnK-UUU
            except ValueError:
                return
            if ex in ("e1", "e2"):
                pairs[(acc, gene)][ex] = "".join(buf).upper()

        with open(fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    _flush()
                    sid, buf = line[1:].split()[0], []
                else:
                    buf.append(line.strip())
        _flush()

    out = defaultdict(list)
    for (acc, gene), d in pairs.items():
        if gene not in INTRON_TRNA_GENES:
            continue
        a, b = d.get("e1"), d.get("e2")
        if not a or not b:
            continue
        if not (MIN_SPLICED <= len(a) + len(b) <= MAX_SPLICED):
            continue
        if not (15 <= len(a) <= 60 and 15 <= len(b) <= 60):
            continue
        out[gene].append((acc, a, b))

    out = dict(out)
    _DONOR_CACHE[key] = out
    return out


def _place(query, window):
    """Best glocal placement of `query` in `window`, as (score, start, end).

    `len(alignments)` is deliberately not called: for a short query against a
    long target the alignment count can be astronomically large, and asking for
    it costs more than the alignment itself.
    """
    if not query or len(window) < len(query):
        return None
    try:
        a = _ALN.align(window, query)[0]
    except (IndexError, ValueError, OverflowError):
        return None
    blocks = a.aligned[0]
    if len(blocks) == 0:
        return None
    return (float(a.score), int(blocks[0][0]), int(blocks[-1][1]))


def _pair_on_strand(window, e1, e2):
    """Place an exon pair on one strand. Returns (score, (s1,t1), (s2,t2))."""
    p1 = _place(e1, window)
    if p1 is None:
        return None
    p2 = _place(e2, window)
    if p2 is None:
        return None
    _, s1, t1 = p1
    _, s2, t2 = p2
    if s2 <= t1:                                  # exon 2 must follow exon 1
        return None
    if not (MIN_INTRON <= s2 - t1 <= MAX_INTRON):
        return None
    if not (MIN_SPLICED <= (t1 - s1) + (t2 - s2) <= MAX_SPLICED):
        return None
    return (p1[0] + p2[0], (s1, t1), (s2, t2))


def _placements(window, donors):
    """Every donor's placement, on both strands, in forward-window offsets.

    Scores are per aligned base so that donor length cancels; a longer donor
    otherwise wins on raw score and then imposes its own exon lengths.
    """
    rcw = _rc(window)
    L = len(window)
    out = []
    for _acc, a, b in donors:
        denom = float(len(a) + len(b)) or 1.0
        for strand in (1, -1):
            r = _pair_on_strand(window if strand == 1 else rcw, a, b)
            if r is None:
                continue
            sc, x1, x2 = r
            if strand == -1:                      # back to forward offsets
                x1, x2 = (L - x2[1], L - x2[0]), (L - x1[1], L - x1[0])
            out.append((sc / denom, x1, x2))
    return out


def _vote(placements):
    """Per-coordinate mode over the better-scoring half of the placements."""
    if not placements:
        return None
    ranked = sorted(placements, key=lambda p: -p[0])
    keep = ranked[:max(3, int(len(ranked) * VOTE_FRACTION))]
    coords = []
    for j in range(4):
        c = Counter((p[1] if j < 2 else p[2])[j % 2] for p in keep)
        coords.append(c.most_common(1)[0][0])
    (s1, t1, s2, t2) = coords
    if not (s1 < t1 < s2 < t2):
        return None
    return [(s1, t1), (s2, t2)]


def _sample(pool, n):
    """`n` donors spread across the pool rather than the first `n` of one clade."""
    if len(pool) <= n:
        return pool
    step = max(1, len(pool) // n)
    return pool[::step][:n]


def refine_intron_trna(trnas, genome_seq, exon_db, glen, diagnostics=None):
    """Replace the coordinates of intron-bearing tRNA in place.

    Writes only `.start`, `.end`, `.exons` and `.notes` on features already in
    the caller's list. It never appends, removes or reorders, so it cannot
    change the inventory -- the same property `apply_trnascan_geometry` relies
    on, and for the same reason: every step that decides which loci exist has
    already run by the time this is called.

    Returns the number of features whose coordinates changed.
    """
    diag = diagnostics if diagnostics is not None else []
    donors_by_gene = load_donors(exon_db)
    if not donors_by_gene:
        diag.append("no usable exon donors; intron refinement skipped")
        return 0

    g = str(genome_seq).upper()
    moved = 0
    seen = skipped_wrapped = no_donor = no_vote = 0

    for f in trnas:
        arcs = sorted(f.exons) if getattr(f, "exons", None) else []
        if len(arcs) != 2:
            continue
        gene = canonical_gene(getattr(f, "gene_name", None))
        if gene is None:
            continue
        seen += 1

        # An origin-crossing locus would need the window assembled from two
        # arcs and the result mapped back across the join. Rare, and getting it
        # subtly wrong would move a feature to the far side of the genome, so
        # it is left alone and counted.
        if f.start > f.end or arcs[0][0] > arcs[-1][1]:
            skipped_wrapped += 1
            continue

        pool = donors_by_gene.get(gene) or []
        if not pool:
            no_donor += 1
            continue

        lo = max(0, arcs[0][0] - PAD)
        hi = min(glen, arcs[-1][1] + PAD)
        if hi - lo < MIN_SPLICED:
            continue

        v = _vote(_placements(g[lo:hi], _sample(pool, N_DONORS)))
        if v is None:
            no_vote += 1
            continue

        new = [(s + lo, e + lo) for s, e in v]
        if new == arcs:
            continue

        old = (f.start, f.end, list(arcs))
        f.exons = new
        f.start = new[0][0]
        f.end = new[-1][1]
        f.notes.append("geometry_source=glocal-exon-vote")
        f.notes.append("geometry_was=%d-%d exons=%s"
                       % (old[0], old[1],
                          ",".join("%d-%d" % (s, e) for s, e in old[2])))
        moved += 1

    if seen:
        parts = ["%d/%d intron-bearing tRNA re-placed" % (moved, seen)]
        if skipped_wrapped:
            parts.append("%d crossed the origin and were left alone" % skipped_wrapped)
        if no_donor:
            parts.append("%d had no donor" % no_donor)
        if no_vote:
            parts.append("%d had no plausible exon pair" % no_vote)
        diag.append("; ".join(parts))
    return moved
