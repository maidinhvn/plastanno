"""Hybrid tRNA geometry: ARAGORN detection, tRNAscan-SE boundaries.

Why this exists. A causal trace over eight loci in NC_037507.1 followed raw
ARAGORN -> aragorn_coords -> selected candidate -> final GenBank and found the
converted coordinate equal to the final coordinate at every step: Plastanno does
not move a tRNA boundary. The 1 nt difference against RefSeq is already present
in ARAGORN's own output, and where ARAGORN agrees with tRNAscan-SE the pipeline
agrees with RefSeq exactly. So this is a detector boundary-convention difference,
and the fix is to take the boundary from the detector whose convention matches
the references — not to trim a nucleotide, which would break the 326 of 865
intron-free tRNAs that already match.

Three things this module must not do, each of which was a real failure mode:

* It must not silently degrade. `run_trnascan` used to return [] when
  tRNAscan-SE was missing or crashed, which in hybrid mode would quietly hand
  back legacy coordinates under a hybrid label. Here a missing or failing
  tRNAscan-SE raises.
* It must not lose an origin-crossing tRNA. tRNAscan-SE reads a linear FASTA and
  cannot see across the origin; on a real plastome it reported trnH-GUG as only
  the fragment after position 1. The genome is therefore scanned twice, rotated
  by half its length, and the second pass supplies what the first cut in half.
* It must not let a name change a coordinate. Geometry and identity are decided
  separately and recorded separately.
"""
import os
import re
import shutil
import subprocess
import tempfile

from ..core import coords as _coords

TRNASCAN_ARGS = ("-O",)          # locked; see benchmark_v3/gate_b/CONFIG_LOCK.md
TRNASCAN_TIMEOUT = 1200

# reciprocal overlap below which two calls are not the same physical locus
MIN_RECIPROCAL = 0.5


class TrnascanUnavailable(RuntimeError):
    """tRNAscan-SE is absent or failed. Hybrid mode is fail-closed on this."""


def _arcs_of(feature):
    ex = getattr(feature, "exons", None)
    if ex:
        return sorted((int(s), int(e)) for s, e in ex)
    return [(int(feature.start), int(feature.end))]


def _span_set(arcs, glen):
    out = set()
    for s, e in arcs:
        if s <= e:
            out |= set(range(s, e))
        else:                                   # crosses the origin
            out |= set(range(s, glen)) | set(range(0, e))
    return out


def reciprocal_overlap(a_arcs, b_arcs, glen):
    """min(overlap/len_a, overlap/len_b) on the circle."""
    A, B = _span_set(a_arcs, glen), _span_set(b_arcs, glen)
    if not A or not B:
        return 0.0
    inter = len(A & B)
    return min(inter / len(A), inter / len(B))


# ── running the detector twice, on a circle ──────────────────────────────────

_LINE = re.compile(
    r"^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+([0-9.]+)")


def _parse(path, glen, shift):
    """tRNAscan-SE tabular output -> [dict]. `shift` un-rotates the coordinates.

    A position q in a frame rotated left by `shift` sat at (q + shift) % glen in
    the original, so un-rotating is a modular add, and a locus that the rotation
    moved across the seam comes back as two arcs.
    """
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        m = _LINE.match(line)
        if not m:
            continue
        _seq, _n, b, e, aa, codon, ib, ie, score = m.groups()
        b, e, ib, ie = int(b), int(e), int(ib), int(ie)
        strand = -1 if b > e else 1
        lo, hi = (b, e) if b <= e else (e, b)
        seg = [(lo - 1, hi)]
        has_intron = bool(ib or ie)
        if has_intron:
            ilo, ihi = (ib, ie) if ib <= ie else (ie, ib)
            seg = [(lo - 1, ilo - 1), (ihi, hi)]
        arcs, touches_cut = [], (lo - 1 == 0 or hi == glen)
        for s, t in seg:
            if t <= s:
                continue
            s2 = (s + shift) % glen
            t2 = (t + shift - 1) % glen + 1
            if s2 < t2:
                arcs.append((s2, t2))
            else:
                arcs.extend([(s2, glen), (0, t2)])
        if not arcs:
            continue
        out.append({
            "aa": aa, "anticodon": codon.upper().replace("T", "U"),
            "strand": strand, "arcs": sorted(arcs), "score": float(score),
            "has_intron": has_intron, "touches_cut": touches_cut,
            "from_rotation": shift,
            "call": "%d-%d%s" % (b, e, " intron %d-%d" % (ib, ie) if has_intron else ""),
        })
    return out


def _run_once(seq, tag):
    fa = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    out = fa.name + ".%s.out" % tag
    try:
        fa.write(">query\n%s\n" % seq)
        fa.close()
        r = subprocess.run(
            ["tRNAscan-SE"] + list(TRNASCAN_ARGS) + ["-q", "-o", out, fa.name],
            capture_output=True, text=True, timeout=TRNASCAN_TIMEOUT)
        if r.returncode != 0:
            raise TrnascanUnavailable(
                "tRNAscan-SE exited %d on the %s frame: %s"
                % (r.returncode, tag, (r.stderr or "")[:300]))
        if not os.path.exists(out):
            raise TrnascanUnavailable(
                "tRNAscan-SE produced no output file on the %s frame" % tag)
        return out, fa.name
    except subprocess.TimeoutExpired:
        raise TrnascanUnavailable(
            "tRNAscan-SE timed out after %d s on the %s frame" % (TRNASCAN_TIMEOUT, tag))


def require_trnascan():
    """Raise unless tRNAscan-SE can be run.

    Called once at the start of a run as well as at the point of use. Since
    hybrid became the default the binary is a hard dependency, and discovering
    that six steps and several minutes into an annotation — after Exonerate, the
    HMM search and reconciliation have all completed — is not an acceptable way
    to learn it. The message names the fix rather than only the fault.
    """
    if not shutil.which("tRNAscan-SE"):
        raise TrnascanUnavailable(
            "tRNAscan-SE is not on PATH.\n"
            "  It is required by --trna-mode hybrid, which is the DEFAULT since "
            "tRNA boundary\n"
            "  accuracy was measured to depend on it (acceptor stem 5.63 -> 6.78 of 7 "
            "on the\n"
            "  intron-free class, validated on 99 held-out genomes).\n"
            "  Install it:   conda install -c bioconda trnascan-se\n"
            "  Or opt out:   --trna-mode legacy   (keeps ARAGORN boundaries)")


def run_trnascan_circular(genome_seq):
    """tRNAscan-SE over a circular genome: original frame plus a half-rotation.

    Returns [dict] in original coordinates. Fail-closed: raises
    TrnascanUnavailable rather than returning a short list, because a silently
    empty result in hybrid mode is indistinguishable from legacy output.
    """
    require_trnascan()
    glen = len(genome_seq)
    shift = glen // 2
    rotated = genome_seq[shift:] + genome_seq[:shift]

    hits, tmp = [], []
    try:
        for seq, sh, tag in ((genome_seq, 0, "orig"), (rotated, shift, "rot")):
            out, fapath = _run_once(seq, tag)
            tmp += [out, fapath]
            hits += _parse(out, glen, sh)
    finally:
        for p in tmp:
            try:
                os.unlink(p)
            except OSError:
                pass

    # Deduplicate by circular position, strand and reciprocal overlap. A locus
    # seen in both frames is kept once, preferring the call that sits wholly
    # inside its own frame: the other one was cut by the sequence end.
    hits.sort(key=lambda h: (h["touches_cut"], -h["score"]))
    kept = []
    for h in hits:
        twin = next((k for k in kept
                     if k["strand"] == h["strand"]
                     and reciprocal_overlap(k["arcs"], h["arcs"], glen) > MIN_RECIPROCAL),
                    None)
        if twin is None:
            kept.append(h)
    return kept


# ── geometry replacement over the legacy inventory ───────────────────────────
#
# The first version rebuilt the merge from scratch and lost things: BLAST-only
# loci never reached the output at all, and ARAGORN calls were dropped on an
# exon-DB overlap without checking strand, which deletes one of two IR copies
# facing opposite ways. Rebuilding an inventory is not what this change is for.
#
# So: legacy decides WHICH loci exist, and this decides only WHERE the
# intron-free ones start and end. Inventory and geometry move independently, and
# a boundary experiment cannot silently become an inventory experiment.

GEOM_TRNASCAN = "tRNAscan-SE"
GEOM_LEGACY = "legacy"

MIN_RECIPROCAL_REPLACE = 0.8      # stricter than the merge default: replacing a
                                  # coordinate needs more confidence than noting
                                  # that two calls touch


def _fmt(arcs):
    return ",".join("%d-%d" % (s, e) for s, e in arcs)


def _identity_compatible(feature, hit):
    """Same amino-acid family, or one side unresolved.

    Anticodon is deliberately NOT required to match: tRNAscan-SE cannot separate
    trnI-CAU from trnM-CAU, and the audit found it calling Met where both
    Plastanno and RefSeq say Ile. Requiring the anticodon would refuse exactly
    the loci where the naming work already settled the answer. Family is checked
    so a geometry is never taken from a different tRNA that happens to overlap.
    """
    gene = (getattr(feature, "gene_name", "") or "")
    core = gene.split("-")[0]
    if core.startswith("trn"):
        core = core[3:]
    fam1 = _ONE_TO_THREE.get(core)
    fam2 = (hit.get("aa") or "").capitalize()
    if not fam1 or not fam2:
        return True                      # unresolved on either side: allow
    if fam1 == fam2:
        return True
    # The CAU anticodon is the one place where a family disagreement is expected
    # rather than disqualifying. trnI-CAU, trnM-CAU and trnfM-CAU share it, and
    # tRNAscan-SE does not separate them -- the Gate B audit found it calling Met
    # at a locus both Plastanno and RefSeq call Ile. Refusing these would skip
    # exactly the loci the Stage 2 naming work already settled, so the boundary
    # is allowed through while the name stays with TRNAIdentity.
    anticodon = (hit.get("anticodon") or "").upper().replace("T", "U")
    CAU_FAMILIES = {"Ile", "Met"}
    if anticodon == "CAU" and fam1 in CAU_FAMILIES and fam2 in CAU_FAMILIES:
        feature._cau_exception = True
        return True
    return False


_ONE_TO_THREE = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "fM": "Met",
}


def _is_multi_exon(feature, glen):
    """Biological exon count, not storage-arc count.

    An intron-free tRNA crossing the origin is stored as two arcs meeting at the
    seam. Counting arcs calls it spliced and freezes its geometry, which is the
    opposite of what should happen: it is exactly the kind of locus whose
    boundary this change exists to correct.
    """
    return _coords.biological_exon_count(_arcs_of(feature), glen) > 1


def _hit_is_multi_exon(hit, glen):
    return _coords.biological_exon_count(hit["arcs"], glen) > 1


def _evidence_digest(feature, hit, glen):
    """A stable fingerprint of the evidence, independent of where the origin is.

    Only rotation-invariant quantities go in: the amino-acid family, the
    anticodon, the spliced length of each side and the strand. No genome
    coordinate appears, because renumbering the circle must not reorder a tie.
    """
    import hashlib
    parts = [
        (getattr(feature, "gene_name", "") or ""),
        str(getattr(feature, "strand", 0)),
        str(_coords.biological_exon_count(_arcs_of(feature), glen)),
        str(sum(e - s for s, e in _arcs_of(feature))),
        (hit.get("aa") or ""), (hit.get("anticodon") or ""),
        str(sum(e - s for s, e in hit["arcs"])),
        "%.3f" % float(hit.get("score") or 0.0),
    ]
    return hashlib.sha1("|".join(parts).encode()).hexdigest()


def _rank_key(feature, hit, ro, score, glen):
    """Deterministic, rotation-invariant ordering: overlap, score, then digest."""
    return (-ro, -float(score or 0.0), _evidence_digest(feature, hit, glen))


def _tiebreak_fraction(feature, hit, glen):
    """The digest folded into a value far below any real weight difference.

    It only decides exact ties; it can never outrank overlap or score.
    """
    d = _evidence_digest(feature, hit, glen)
    return int(d[:8], 16) / float(0xFFFFFFFF) * 1e-6


def global_match(features, hits, glen, diagnostics=None):
    """One-to-one assignment, maximum cardinality first. Never greedy.

    A greedy pass makes the result depend on input order whenever two candidates
    compete for the same partner. Ranking is: number of pairs, then identity
    agreement, then overlap, then tRNAscan score, then a positional tie-break
    that does not depend on list order.
    """
    diag = diagnostics if diagnostics is not None else []
    cand = []
    for i, f in enumerate(features):
        fa = _arcs_of(f)
        for j, h in enumerate(hits):
            if h["strand"] != f.strand:
                continue
            if not _identity_compatible(f, h):
                continue
            ro = reciprocal_overlap(fa, h["arcs"], glen)
            if ro < MIN_RECIPROCAL_REPLACE:
                continue
            cand.append((i, j, ro, h.get("score", 0.0)))
    if not cand:
        return []

    try:
        from scipy.optimize import linear_sum_assignment
        import numpy as np
    except ImportError:
        linear_sum_assignment = None

    if linear_sum_assignment is None:
        # Deterministic fallback: sort by the same key and take greedily. Still
        # order-independent because the key is total and breaks ties on position.
        cand.sort(key=lambda c: _rank_key(features[c[0]], hits[c[1]], c[2], c[3], glen))
        seen_f, seen_h, out = set(), set(), []
        for i, j, ro, sc in cand:
            if i in seen_f or j in seen_h:
                continue
            seen_f.add(i); seen_h.add(j); out.append((i, j, ro))
        diag.append("scipy unavailable; used the deterministic greedy fallback")
        return out

    import numpy as np
    n, m = len(features), len(hits)
    BIG = 1e6
    cost = np.zeros((n, m))
    allowed = {}
    for i, j, ro, sc in cand:
        allowed[(i, j)] = (ro, sc)
        # cardinality dominates: every legal pair is worth more than any weight
        # The final term must not be a genome coordinate: rotating the genome
        # renumbers every position, so a coordinate tie-break silently makes the
        # assignment rotation-dependent. It is an evidence digest instead --
        # identity, spliced length and score, all invariant under rotation.
        cost[i][j] = -(BIG + ro * 1000.0 + sc + _tiebreak_fraction(features[i], hits[j], glen))
    ri, cj = linear_sum_assignment(cost)
    out = [(int(i), int(j), allowed[(int(i), int(j))][0])
           for i, j in zip(ri, cj) if (int(i), int(j)) in allowed]

    # A digest breaks exact ties deterministically, but determinism is not the
    # same as correctness: if two pairings are indistinguishable on strand,
    # identity, overlap and score, the digest has picked one arbitrarily and
    # saying so is the honest outcome. Report it rather than let an arbitrary
    # choice pass as a resolved one.
    chosen = {(i, j) for i, j, _ in out}
    ambiguous_pairs = set()
    for i, j, ro in out:
        rivals = [(i2, j2) for (i2, j2), (ro2, sc2) in allowed.items()
                  if (i2, j2) not in chosen
                  and (i2 == i or j2 == j)
                  and abs(ro2 - allowed[(i, j)][0]) < 1e-12
                  and abs(sc2 - allowed[(i, j)][1]) < 1e-12]
        if rivals:
            diag.append("AMBIGUOUS_MATCH: feature %d / hit %d tied with %s on "
                        "overlap and score; legacy coordinates kept" % (i, j, rivals))
            try:
                features[i]._ambiguous_match = True
            except AttributeError:
                pass
            ambiguous_pairs.add((i, j))

    # An arbitrary choice must not become a coordinate. Where nothing
    # biological separates two candidate partners, the replacement is declined
    # and the legacy geometry stands -- the conservative direction, since legacy
    # is what the frozen matrix already contains.
    out = [t for t in out if (t[0], t[1]) not in ambiguous_pairs]
    out.sort(key=lambda t: t[0])
    return out


def replace_geometry_only(feature, hit, glen):
    """Move a feature's boundaries onto the hit's, preserving circular form.

    `start = arcs[0][0]; end = arcs[-1][1]` is wrong for an origin-crossing
    locus: the arcs are sorted arithmetically, so the first begins at 0 and the
    last ends at genome_len, and that assignment turns a 73 nt tRNA into the
    whole genome. A wrapped feature keeps start > end, which is how
    coords.is_wrapped recognises it.
    """
    arcs = sorted(hit["arcs"])
    wrapped = len(arcs) >= 2 and arcs[0][0] == 0 and arcs[-1][1] == glen
    if wrapped:
        feature.start = arcs[-1][0]          # the tail, before the origin
        feature.end = arcs[0][1]             # the head, after it
    else:
        feature.start = arcs[0][0]
        feature.end = arcs[-1][1]
    if getattr(feature, "exons", None):
        feature.exons = [tuple(x) for x in arcs]
    return feature


def apply_trnascan_geometry(base, trnascan, glen, diagnostics=None,
                            intron_pass=False):
    """Take intron-free boundaries from tRNAscan-SE; change nothing else.

    `base` is the legacy tRNA inventory. Every locus in it comes back — a locus
    is never dropped here, and no locus is added: tRNAscan-only calls are
    recorded as diagnostics and left out, so this round measures boundary
    accuracy without also moving the inventory.
    """
    diag = diagnostics if diagnostics is not None else []
    # Multi-exon calls are excluded on BOTH sides, and the reason is stronger
    # than "tRNAscan-SE does not model plastid group II introns". Admitting them
    # was measured (benchmark_v3/intron_geometry/): the intron-bearing class
    # barely moved, 19.6% to 18.2% exact, while the intron-FREE class collapsed
    # from 86.6% to 37.2% and its acceptor stem from 6.77 to 5.70 of 7.
    # global_match is a one-to-one assignment over every locus, so adding spliced
    # candidates re-partitions the whole thing: a spliced call takes an
    # intron-free feature and gives it a two-arc geometry. The exclusion protects
    # the path where 87% of the accuracy lives.
    eligible, frozen = [], []
    for f in base:
        (frozen if _is_multi_exon(f, glen) else eligible).append(f)

    usable  = [h for h in trnascan if not _hit_is_multi_exon(h, glen)]
    spliced = [h for h in trnascan if _hit_is_multi_exon(h, glen)]
    if not intron_pass:
        for h in spliced:
            diag.append("trnascan spliced call ignored (group II introns are not "
                        "modelled by it): %s" % _fmt(h["arcs"]))

    def _one_pass(feats, hits, label):
        """Match one pool against one pool and annotate the outcome.

        Factored out so the intron-bearing loci can have their OWN assignment
        over their OWN candidates. Merging them into a single assignment was
        measured and rejected (benchmark_v3/intron_geometry/): global_match is
        one-to-one over everything it is given, so a spliced candidate took an
        intron-free feature and intron-free exactness fell 86.6% -> 37.2%.
        Two disjoint pools cannot do that to each other.
        """
        pairs = global_match(feats, hits, glen, diagnostics=diag)
        matched_f = {i for i, _j, _ro in pairs}
        matched_h = {j for _i, j, _ro in pairs}

        for i, j, _ro in pairs:
            f, h = feats[i], hits[j]
            before = _fmt(_arcs_of(f))
            replace_geometry_only(f, h, glen)
            f.notes.append("geometry_source=%s" % GEOM_TRNASCAN)
            f.notes.append("identity_source=%s" % _identity_source_of(f))
            if getattr(f, "_cau_exception", False):
                f.notes.append("identity_match=CAU_AMBIGUOUS")
            if getattr(f, "_ambiguous_match", False):
                f.notes.append("geometry_match=AMBIGUOUS_MATCH")
            f.notes.append("aragorn_call=%s" % before)
            f.notes.append("trnascan_call=%s" % _fmt(h["arcs"]))

        for k, f in enumerate(feats):
            if k in matched_f:
                continue
            f.notes.append("geometry_source=%s" % GEOM_LEGACY)
            if getattr(f, "_ambiguous_match", False):
                f.notes.append("geometry_match=AMBIGUOUS_MATCH_DECLINED")
            f.notes.append("identity_source=%s" % _identity_source_of(f))
            f.notes.append("aragorn_call=%s" % _fmt(_arcs_of(f)))
            f.notes.append("trnascan_call=none")

        for j, h in enumerate(hits):
            if j not in matched_h:
                diag.append("trnascan-only, not added this round (%s): %s-%s %s "
                            "strand %+d" % (label, h.get("aa"), h.get("anticodon"),
                                            _fmt(h["arcs"]), h["strand"]))

    _one_pass(eligible, usable, "intron-free")

    if intron_pass:
        _one_pass(frozen, spliced, "intron-bearing")
    else:
        for f in frozen:
            f.notes.append("geometry_source=%s" % GEOM_LEGACY)
            f.notes.append("identity_source=%s" % _identity_source_of(f))
            f.notes.append("aragorn_call=%s" % _fmt(_arcs_of(f)))
            f.notes.append("trnascan_call=none")

    return frozen + eligible


# ── where a name came from ───────────────────────────────────────────────────

IDENT_ARAGORN = "ARAGORN"
IDENT_BLAST = "BLAST"
IDENT_EXON_DB = "EXON_DB"
IDENT_UNKNOWN = "unknown"


def tag_identity_sources(aragorn, blast, intron):
    """Stamp each candidate with where its NAME came from, before the legacy
    merge collapses the three lists into one.

    Reading it back afterwards is not possible: the merge records the winning
    caller only inside its alternatives list, not on the surviving feature. So
    the tag is attached first, and the hybrid step reports what the tag says
    rather than asserting a single source for everything.
    """
    for f in intron:
        f._identity_src = IDENT_EXON_DB
    for f in aragorn:
        f._identity_src = IDENT_ARAGORN
    for f in blast:
        f._identity_src = IDENT_BLAST


def _identity_source_of(feature):
    return getattr(feature, "_identity_src", IDENT_UNKNOWN)
