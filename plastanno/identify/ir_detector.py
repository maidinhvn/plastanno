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
        if p[6] != "minus" or length < min_ir_len:
            continue
        hsps.append((qs, qe, ss, se, length, float(p[5])))
    return hsps


def detect_ir_boundaries(
    genome_seq: str,
    min_ir_len: int = 10000,
) -> Optional[Dict[str, Tuple[int, int]]]:
    """Detect LSC/IRb/SSC/IRa boundaries. Returns None if no IR pair is found."""
    n = len(genome_seq)
    hsps = _self_blastn(genome_seq, min_ir_len)
    if not hsps:
        return None

    # longest minus-strand HSP = the IR pair
    qs, qe, ss, se, length, pident = max(hsps, key=lambda h: h[4])

    # two repeat intervals (blast 1-based inclusive -> 0-based half-open)
    iv1 = (min(qs, qe) - 1, max(qs, qe))
    iv2 = (min(ss, se) - 1, max(ss, se))
    low, high = sorted([iv1, iv2], key=lambda x: x[0])

    irb = (low[0], low[1])
    ssc = (low[1], high[0])
    ira = (high[0], high[1])

    # LSC wraps from end-of-IRa to start-of-IRb; for a sequence rotated to start
    # in the LSC (IRa ends at the genome end), that simplifies to [0, IRb_start).
    if (n - ira[1]) <= 100:
        lsc = (0, irb[0])
    else:
        lsc = (ira[1], irb[0])  # start > end signals a wrap (non-rotated genome)

    return {"LSC": lsc, "IRb": irb, "SSC": ssc, "IRa": ira}


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
