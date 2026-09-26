"""
IR / LSC / SSC boundary detection via self-BLASTN (clean rewrite).

Plastomes are usually quadripartite: LSC - IRb - SSC - IRa, where IRb and IRa
are large inverted (reverse-complement) repeats. The IR pair is detected as the
longest MINUS-strand self-BLASTN HSP that is at least `min_ir_len` long; the four
regions are then derived from the two repeat intervals.

Returns a dict of 0-based half-open (start, end) tuples for LSC, IRb, SSC, IRa,
or None when no IR pair >= min_ir_len is found (e.g. IR-lacking plastomes).

Convention (standard NCBI orientation, sequence rotated to start in the LSC):
    LSC = [0, IRb_start)
    IRb = [IRb_start, IRb_end)
    SSC = [IRb_end, IRa_start)
    IRa = [IRa_start, IRa_end)        with IRa_end == genome_len

Nothing here may assume that orientation, though: where position 1 falls and which
strand was deposited are the submitter's choices, and both vary. Two consequences
are handled explicitly.

A repeat that crosses the origin is reported by BLAST as two shorter HSPs, one at
each end of the sequence. Taking the longest single HSP then finds only the larger
half: on a test genome rotated into IRb, the repeat came back 16,020 bp instead of
26,492, ten kilobases of it were reassigned to the LSC, and the genes there
(rpl2, rpl23, ycf2) were searched in the wrong compartment. Collinear HSPs are
therefore joined across the origin before the longest is chosen.

Which single-copy region is which is decided by SIZE, not by position. The old code
called the lower-coordinate one the LSC, which is only true in the deposited
orientation; on a reverse-complemented genome it swapped the two, labelling a
17.9 kb region "LSC" and an 84.0 kb region "SSC" — and every membership test that
follows, rps12's 5' exon among them, then asked about the wrong compartment.

The default threshold is 10 kb (per the v1 lesson; 5 kb produced spurious calls).
"""
import os
import subprocess
import tempfile
from typing import Dict, Optional, Tuple


def _self_blastn(genome_seq: str, min_ir_len: int):
    """Return minus-strand HSPs as (qs, qe, ss, se, length, pident)."""
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    try:
        tmp.write(">g\n%s\n" % genome_seq)
        tmp.close()
        proc = subprocess.run(
            ["blastn", "-query", tmp.name, "-subject", tmp.name,
             "-dust", "no", "-evalue", "1e-20", "-word_size", "11",
             "-outfmt", "6 qstart qend sstart send length pident sstrand"],
            capture_output=True, text=True, timeout=300,
        )
    finally:
        try:
            os.unlink(tmp.name)
        except OSError:
            pass

    hsps = []
    for line in proc.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 7:
            continue
        qs, qe, ss, se, length = int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4])
        # Short pieces are kept here and only filtered after the join below: the two
        # halves of a repeat split by the origin can each be under the threshold.
        if p[6] != "minus":
            continue
        hsps.append((qs, qe, ss, se, length, float(p[5])))
    return hsps


# Slack allowed when deciding that two HSPs are the same repeat continuing across
# the origin. BLAST rarely takes an alignment to the last base, and the two halves
# are reported independently, so a few tens of bases of gap or overlap at the join
# is normal. Far smaller than any real gap between distinct repeats.
_JOIN_SLACK = 200


# The full IUPAC ambiguity alphabet. Complementing only ACGTN leaves R, Y, K, M,
# B, D, H, V and S/W untouched, so a repeat containing a single ambiguity code
# would give a different key from its own reverse complement — defeating the one
# thing this key exists to guarantee. Deposited plastomes do contain these codes.
_RCMAP = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                       "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def _repeat_key(seq, qs, qe):
    """A key for one repeat that a rotation or a strand flip cannot change.

    The repeat's own bases, taken in whichever orientation sorts first, so the same
    physical repeat gives the same key however the sequence was deposited.
    """
    import hashlib
    n = len(seq)
    a, b = (qs - 1) % n, ((qe - 1) % n) + 1
    sub = seq[a:b] if a < b else seq[a:] + seq[:b]
    # U is normalised to T first: complementing U gives A, whose complement is T,
    # so without this the key of a sequence containing U would not survive a round
    # trip and the reverse-complement symmetry this function exists for would fail.
    sub = sub.upper().replace("U", "T")
    rc = sub.translate(_RCMAP)[::-1]
    return hashlib.sha256(min(sub, rc).encode()).hexdigest()


def _join_across_origin(hsps, n):
    """Join HSP halves that are one repeat continuing across position 1.

    A minus-strand HSP pairs an ascending query interval with a descending subject
    interval. When the repeat crosses the origin, BLAST reports a head piece whose
    query starts at 1 and a tail piece whose query ends at n, and their subject
    intervals meet at the corresponding place. Joined, the pair is one wrapped
    repeat, described here by a query interval with start > end.
    """
    out = list(hsps)
    for a in hsps:
        for b in hsps:
            if a is b:
                continue
            aqs, aqe, ass, ase, alen, apid = a
            bqs, bqe, bss, bse, blen, bpid = b
            # a ends at the sequence end, b begins at the start
            if aqe < n - _JOIN_SLACK or bqs > 1 + _JOIN_SLACK:
                continue
            # their partners must meet too: a's subject continues into b's
            lo_a, hi_a = min(ass, ase), max(ass, ase)
            lo_b, hi_b = min(bss, bse), max(bss, bse)
            if abs(lo_a - hi_b) > _JOIN_SLACK and abs(lo_b - hi_a) > _JOIN_SLACK:
                continue
            out.append((aqs, bqe + n, min(lo_a, lo_b), max(hi_a, hi_b),
                        alen + blen, (apid * alen + bpid * blen) / (alen + blen)))
    return out


def detect_ir_boundaries(
    genome_seq: str,
    min_ir_len: int = 10000,
) -> Optional[Dict[str, Tuple[int, int]]]:
    """Detect LSC/IRb/SSC/IRa boundaries. Returns None if no IR pair is found."""
    n = len(genome_seq)
    raw = _self_blastn(genome_seq, min_ir_len)
    if not raw:
        return None
    hsps = [h for h in _join_across_origin(raw, n) if h[4] >= min_ir_len]
    if not hsps:
        return None

    # Longest minus-strand HSP = the IR pair. The two HSPs describing one pair have
    # IDENTICAL length AND identity by definition, so this tie always happens — but
    # the two are the same pair seen from either end (query/subject swapped), and
    # _regions below tries both assignments anyway, so it does not matter which is
    # taken. What must not decide it is a query coordinate: that changes under
    # rotation and reverse complement. Ties beyond that are broken by the repeat's
    # own sequence, normalised so that a repeat and its reverse complement give the
    # same key.
    def _rank(h):
        return (h[4], h[5], _repeat_key(genome_seq, h[0], h[1]))
    qs, qe, ss, se, length, pident = max(hsps, key=_rank)

    # two repeat intervals (blast 1-based inclusive -> 0-based half-open). A joined
    # HSP has qe > n; taking it mod n gives the wrapped (start > end) form.
    iv1 = ((qs - 1) % n, ((qe - 1) % n) + 1)
    iv2 = ((min(ss, se) - 1) % n, ((max(ss, se) - 1) % n) + 1)

    def _len(iv):
        return (iv[1] - iv[0]) if iv[1] > iv[0] else (n - iv[0] + iv[1])

    def _mid(iv):
        return (iv[0] + _len(iv) // 2) % n

    # IRb is the repeat whose downstream neighbour is the SMALLER single-copy
    # region; naming them by coordinate order is only right in the deposited
    # orientation. Try both assignments and keep the one where the region
    # following IRb is the smaller of the two.
    def _regions(first, second):
        return {"IRb": first, "SSC": (first[1] % n, second[0]),
                "IRa": second, "LSC": (second[1] % n, first[0])}

    cand_a = _regions(iv1, iv2)
    cand_b = _regions(iv2, iv1)
    best = min((cand_a, cand_b), key=lambda c: _len(c["SSC"]))

    # A region of length 0 means the two repeats abut; report the whole remainder
    # rather than an empty compartment.
    for key in ("SSC", "LSC"):
        if _len(best[key]) == 0:
            best[key] = (best[key][0], best[key][0])

    # An LSC that happens to begin at 0 is written [0, IRb_start) rather than
    # (n, IRb_start), so downstream code sees the familiar non-wrapping form.
    if best["LSC"][0] % n == 0:
        best["LSC"] = (0, best["LSC"][1])
    return {"LSC": best["LSC"], "IRb": best["IRb"],
            "SSC": best["SSC"], "IRa": best["IRa"]}


def _region_len(region: Tuple[int, int], n: int) -> int:
    s, e = region
    return (e - s) if e >= s else (n - s + e)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("usage: python ir_detector.py genome.fasta")
        sys.exit(1)
    seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
    n = len(seq)
    res = detect_ir_boundaries(seq)
    if res is None:
        print("No IR pair >= 10 kb found (possibly IR-lacking). Genome = %d bp" % n)
    else:
        print("Genome length: %d bp" % n)
        for k in ("LSC", "IRb", "SSC", "IRa"):
            s, e = res[k]
            print("  %-4s %9d - %-9d  (%d bp)" % (k, s, e, _region_len((s, e), n)))
