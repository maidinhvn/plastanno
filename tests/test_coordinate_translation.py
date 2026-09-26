#!/usr/bin/env python3
"""GATE A — coordinate-translation fidelity.

Does Plastanno carry coordinates correctly along
    ARAGORN (1-based inclusive) -> internal (0-based half-open) -> GenBank / GFF3 ?

This gate is TECHNICAL ONLY. It cannot say whether ARAGORN picked the
biologically right boundary; that is Gate B. It detects off-by-one, strand
inversion, lost or duplicated exons, wrong transcript order, and mishandled
origin crossing.

Eight fixtures: +/- strand x (no intron | intron) x (no origin crossing |
origin crossing), with coordinates worked out by hand below and NOT produced by
calling the code under test.

Run directly; exit 0 = all pass. No skips.
"""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.identify.engine_b import exons_from_aragorn
from plastanno.core.feature import Feature
from plastanno.output import writers as W

FAIL, RUN = [], [0]
L = 1000                      # genome length for every fixture


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-62s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


def span(arcs):
    """Nucleotides supported, summed over arcs."""
    return sum(e - s for s, e in arcs)


# ── the eight fixtures, hand-computed ────────────────────────────────────────
# ARAGORN prints 1-based inclusive [start,end]; a wrapped feature as [high,low].
#   name, s1, e1, strand, intron(off,len) or None, expected 0-based half-open arcs,
#   expected nucleotide support
FIXTURES = [
    ("plus  / no intron / no origin", 101, 172, +1, None,
     [(100, 172)], 72),
    ("minus / no intron / no origin", 101, 172, -1, None,
     [(100, 172)], 72),
    # 101..272 inclusive = 172 nt; intron 100 nt at offset 30 -> mature 72
    #   + strand: 5' exon 101..130 (30 nt), 3' exon 231..272 (42 nt)
    ("plus  / intron    / no origin", 101, 272, +1, (30, 100),
     [(100, 130), (230, 272)], 72),
    #   - strand: 5' exon at the HIGH end 243..272 (30 nt), 3' exon 101..142 (42)
    ("minus / intron    / no origin", 101, 272, -1, (30, 100),
     [(100, 142), (242, 272)], 72),
    # wrapped: ARAGORN writes [971,40]; 971..1000 is 30 nt, 1..40 is 40 nt
    ("plus  / no intron / origin crossing", 971, 40, +1, None,
     [(970, 1000), (0, 40)], 70),
    ("minus / no intron / origin crossing", 971, 40, -1, None,
     [(970, 1000), (0, 40)], 70),
    # wrapped AND intron-bearing, with an EXON split across the origin.
    # Unrolled 951..1120 = 170 nt; intron 60 -> mature 110.
    #  + strand, offset 80: 5' exon 951..1030 (80 nt) STRADDLES the origin and
    #    becomes two storage arcs; 3' exon 1091..1120 (30 nt) -> 90..120.
    ("plus  / intron    / origin crossing", 951, 120, +1, (80, 60),
     [(0, 30), (90, 120), (950, 1000)], 110),
    #  - strand, offset 30: 5' exon is at the high end 1091..1120 -> 90..120;
    #    3' exon 951..1030 (80 nt) straddles the origin.
    ("minus / intron    / origin crossing", 951, 120, -1, (30, 60),
     [(0, 30), (90, 120), (950, 1000)], 110),
]

print("--- 1. ARAGORN -> internal 0-based half-open ---")
internal = {}
for name, s1, e1, st, intr, want_arcs, want_nt in FIXTURES:
    # call the PRODUCTION decision, not a copy of its branching
    arcs, exons, wrapped = exons_from_aragorn(
        s1, e1, st, (intr[0] if intr else None), (intr[1] if intr else None), L)
    got = exons if exons else arcs
    check("%-38s arcs" % name, [tuple(x) for x in got], [tuple(x) for x in want_arcs])
    check("%-38s nucleotides" % name, span(got), want_nt)
    internal[name] = (got, st, wrapped, exons is not None)

print("--- 2. no nucleotide is lost or duplicated at a boundary ---")
for name, s1, e1, st, intr, want_arcs, want_nt in FIXTURES:
    got, _st, wrapped, _has = internal[name]
    flat = []
    for s, e in got:
        flat.extend(range(s, e))
    check("%-38s positions are distinct" % name, len(flat), len(set(flat)))
    check("%-38s count matches the span" % name, len(flat), want_nt)

print("--- 3. an intron survives the origin ---")
# The previous version of this gate ASSERTED that a wrapped tRNA loses its
# intron. That locked a defect in: the same tRNA was a two-exon model in one
# rotation and a one-exon model in another, and a correct fix would have failed
# the test. The origin is a presentation choice, not biology.
for name in ("plus  / intron    / origin crossing",
             "minus / intron    / origin crossing"):
    got, _st, wrapped, has_exons = internal[name]
    check("%-38s is flagged wrapped" % name, wrapped, True)
    check("%-38s KEEPS its intron" % name, has_exons, True)
    check("%-38s 3 storage arcs for 2 biological exons" % name, len(got), 3)
    # the split exon's two arcs meet exactly at the seam, which is how the
    # biological exon count stays recoverable from the storage form
    touches = [(s_, e_) for s_, e_ in got if s_ == 0 or e_ == L]
    check("%-38s the split exon meets at the seam" % name, len(touches), 2)

print("--- 4. internal -> GenBank/GFF3 and back, independently ---")
from Bio import SeqIO
import random
rng = random.Random(5)
GEN = "".join(rng.choice("ACGT") for _ in range(L))


def roundtrip(name):
    got, st, wrapped, has_exons = internal[name]
    if not got:
        raise ValueError("no coordinates produced")
    f = Feature(gene_name="trnG-UCC", gene_type="tRNA", product="tRNA-Gly",
                start=min(s for s, _ in got), end=max(e for _, e in got),
                strand=st)
    f.exons = [tuple(x) for x in got]
    f.has_intron = len(got) > 1
    d = tempfile.mkdtemp()
    W.write_all(annotations=[f], genome_seq=GEN, accession="T", genome_len=L,
                ir_boundaries={"LSC": (0, L)}, relatives=[], out_dir=d,
                prefix="T", no_plot=True)
    rec = next(SeqIO.parse(open(os.path.join(d, "T.gb")), "genbank"))
    gb = None
    for x in rec.features:
        if x.type == "tRNA":
            gb = ([(int(p.start), int(p.end)) for p in x.location.parts],
                  1 if x.location.strand != -1 else -1)
    gff = []
    gstrand = None
    for line in open(os.path.join(d, "T.gff3")):
        p = line.split("\t")
        if len(p) > 6 and p[2] == "tRNA":
            gff.append((int(p[3]) - 1, int(p[4])))
            gstrand = 1 if p[6] == "+" else -1
    # NOTE: the ORDER is returned unsorted as well. An earlier version of this
    # gate compared sorted() arcs only, and so could not see a reversed
    # transcript order -- one of the properties it exists to protect.
    return gb, (sorted(gff), gstrand), gff


for name, _s1, _e1, st, _i, want_arcs, want_nt in FIXTURES:
    try:
        (gb_arcs, gb_st), (gff_arcs, gff_st), gff_order = roundtrip(name)
    except Exception as e:
        # a degenerate or missing exon model is a FAILURE of this fixture, not a
        # reason for the gate to abort and leave the rest unexamined
        check("%-38s round-trips at all" % name, "%s" % type(e).__name__, "no error")
        continue
    check("%-38s GenBank arcs" % name, sorted(gb_arcs), sorted(tuple(x) for x in want_arcs))
    check("%-38s GenBank strand" % name, gb_st, st)
    check("%-38s GFF3 arcs" % name, gff_arcs, sorted(tuple(x) for x in want_arcs))
    check("%-38s GFF3 strand" % name, gff_st, st)
    check("%-38s GenBank nt == GFF3 nt == expected" % name,
          (span(gb_arcs), span(gff_arcs)), (want_nt, want_nt))
    check("%-38s exon count preserved" % name, len(gb_arcs), len(want_arcs))
    # TRANSCRIPT ORDER. On the plus strand transcription runs along ascending
    # coordinates, on the minus strand descending. A feature crossing the origin
    # is exempt: its order is set by the seam, not by coordinate sort.
    # "origin" not in name would also exclude "no origin" -- the substring
    # matches both. Match the crossing case explicitly.
    if len(gb_arcs) > 1 and "origin crossing" not in name:
        want_order = sorted(gb_arcs, reverse=(st == -1))
        check("%-38s GenBank transcript order" % name, gb_arcs, want_order)
        check("%-38s GFF3 transcript order" % name, gff_order,
              sorted(gff_order, reverse=(st == -1)))

print("--- 5. ROTATION INVARIANCE: the same tRNA, wrapped and unwrapped ---")
# The decisive check. One intron tRNA is expressed twice: at a position where it
# does not cross the origin, and after rotating the genome so that it does.
# Mapped back, the two must agree on spliced sequence, exon model and strand.
from plastanno.core import coords as _coords

def build_feature(arcs, exons, wrapped, strand):
    """Build the Feature exactly as run_aragorn does.

    A wrapped feature takes start/end from the ARCS, not from min/max of the
    exons. When the INTRON spans the origin, no exon touches the seam, so
    min/max silently unwraps the feature and transcript_parts then orders the
    exons backwards -- which is what an earlier version of this test did before
    blaming the pipeline for the resulting sequence mismatch.
    """
    if wrapped:
        f_start, f_end = arcs[0][0], arcs[-1][1]
    elif exons:
        f_start, f_end = exons[0][0], exons[-1][1]
    else:
        f_start, f_end = arcs[0][0], arcs[-1][1]
    f = Feature(gene_name="t", gene_type="tRNA", strand=strand,
                start=f_start, end=f_end)
    f.exons = [tuple(x) for x in (exons or arcs)]
    f.has_intron = len(f.exons) > 1
    return f


for st in (+1, -1):
    lab = "plus" if st == 1 else "minus"
    # unrotated: 1-based 101..272, intron i(30,100)
    arcs_u, ex_u, wrapped_u = exons_from_aragorn(101, 272, st, 30, 100, L)
    check("%-5s unrotated is not wrapped" % lab, wrapped_u, False)
    # rotate by 200: genomic p -> (p - 200) mod L. 0-based 100..272 -> 900..72
    DELTA = 200
    rot_genome = GEN[DELTA:] + GEN[:DELTA]
    arcs_w, ex_w, wrapped_w = exons_from_aragorn(901, 72, st, 30, 100, L)
    check("%-5s rotated IS wrapped" % lab, wrapped_w, True)
    check("%-5s rotated still has an exon model" % lab, ex_w is not None, True)
    if ex_w is None or ex_u is None:
        # without an exon model the remaining comparisons are meaningless; record
        # them as failures rather than crashing and leaving the strand untested
        for what in ("mapped back == unrotated exon model", "same nucleotide support",
                     "spliced sequence identical", "the wrapped feature is recognised"):
            check("%-5s %s" % (lab, what), "no exon model", "a comparable model")
        continue
    # map the rotated exons back onto the original coordinate frame
    back = []
    for a, b in ex_w:
        a2, b2 = (a + DELTA) % L, (b + DELTA) % L
        if b2 == 0:
            b2 = L
        if a2 < b2:
            back.append((a2, b2))
        else:
            back.append((a2, L)); back.append((0, b2))
    # merge arcs that meet at the seam, which is where the split exon rejoins
    back.sort()
    merged = []
    for a, b in back:
        if merged and a == merged[-1][1]:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append([a, b])
    merged = [tuple(x) for x in merged]
    check("%-5s mapped back == unrotated exon model" % lab, merged, [tuple(x) for x in ex_u])
    check("%-5s same nucleotide support" % lab,
          sum(b - a for a, b in merged), sum(b - a for a, b in ex_u))
    # and the spliced sequence itself, read from each genome in its own frame
    fu = build_feature(arcs_u, ex_u, wrapped_u, st)
    fw = build_feature(arcs_w, ex_w, wrapped_w, st)
    check("%-5s the wrapped feature is recognised as wrapped" % lab,
          _coords.is_wrapped(fw), True)
    # INDEPENDENT ORACLE. The expected mature sequence is cut BY HAND from the
    # unrotated genome using the known exon coordinates, with a hand-written
    # reverse complement. Comparing the production extractor against itself in
    # two representations would pass even if it were symmetrically wrong.
    def hand_rc(x):
        comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
        return "".join(comp[c] for c in reversed(x))

    # exons on the unrotated genome, written out literally rather than computed
    if st == 1:
        oracle = GEN[100:130] + GEN[230:272]          # 5' exon then 3' exon
    else:
        oracle = hand_rc(GEN[242:272]) + hand_rc(GEN[100:142])
    check("%-5s the hand-cut oracle is the mature length" % lab, len(oracle), 72)

    s_u = _coords.extract(GEN, fu, L)
    s_w = _coords.extract(rot_genome, fw, L)
    check("%-5s unrotated output == hand-cut oracle" % lab, s_u, oracle)
    check("%-5s rotated output   == hand-cut oracle" % lab, s_w, oracle)
    check("%-5s the two representations agree" % lab, s_u, s_w)

print("--- 6. three storage arcs are still TWO biological exons ---")
from plastanno.core.coords import biological_exon_runs, biological_exon_count
check("an exon split at the seam counts once",
      biological_exon_count([(0, 30), (90, 120), (950, 1000)], L), 2)
check("  and the two halves are grouped together",
      biological_exon_runs([(0, 30), (90, 120), (950, 1000)], L),
      [[(90, 120)], [(950, 1000), (0, 30)]])
check("two ordinary exons count as two",
      biological_exon_count([(100, 130), (230, 272)], L), 2)
check("exons abutting AWAY from the seam are NOT merged",
      biological_exon_count([(100, 130), (130, 200)], L), 2)
check("a single arc is one exon", biological_exon_count([(100, 172)], L), 1)
check("a wrapped single exon is one exon",
      biological_exon_count([(0, 40), (970, 1000)], L), 1)
check("without a genome length nothing is merged",
      biological_exon_count([(0, 30), (950, 1000)], None), 2)
# the two wrapped-intron fixtures must read as two exons, not three
for name in ("plus  / intron    / origin crossing",
             "minus / intron    / origin crossing"):
    got, _st, _w, _h = internal[name]
    check("%-38s reads as 2 biological exons" % name,
          biological_exon_count(got, L), 2)

print("--- 7. naming may change; geometry may not ---")
base = internal["plus  / intron    / no origin"][0]
geoms = []
for gene, prod in (("trnG-UCC", "tRNA-Gly"), ("trnG-GCC", "tRNA-Gly"),
                   ("trnG", "tRNA-Gly"), ("trnX", "tRNA-OTHER")):
    f = Feature(gene_name=gene, gene_type="tRNA", product=prod,
                start=base[0][0], end=base[-1][1], strand=1)
    f.exons = [tuple(x) for x in base]
    d = tempfile.mkdtemp()
    W.write_all(annotations=[f], genome_seq=GEN, accession="T", genome_len=L,
                ir_boundaries={"LSC": (0, L)}, relatives=[], out_dir=d, prefix="T",
                no_plot=True)
    rec = next(SeqIO.parse(open(os.path.join(d, "T.gb")), "genbank"))
    for x in rec.features:
        if x.type == "tRNA":
            geoms.append((sorted((int(p.start), int(p.end)) for p in x.location.parts),
                          1 if x.location.strand != -1 else -1))
check("four different identities, one geometry", len(set(map(str, geoms))), 1)
check("and it is the expected geometry", geoms[0][0],
      sorted(tuple(x) for x in base))

print()
print("checks run: %d" % RUN[0])
if FAIL:
    print("%d FAILED:" % len(FAIL))
    for x in FAIL:
        print("   - %s" % x)
    sys.exit(1)
print("all %d checks passed" % RUN[0])
sys.exit(0)
