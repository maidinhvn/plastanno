#!/usr/bin/env python3
"""Glocal exon refinement for intron-bearing tRNA.

The measured claim lives in `benchmark_v3/intron_glocal/`. What is pinned here
is the behaviour that measurement depends on, and the things that would make it
silently wrong:

  * the glocal orientation. Biopython's end-gap scores are easy to swap, and
    swapped they do not fail -- they place a correct exon at score -22 in three
    scattered blocks instead of +20 in one. The first version of this module
    reported 0 of 8 loci for exactly that reason.
  * the donor vote. A single donor imposes its own exon lengths, and raw score
    prefers the longer donor, so a long outlier must not be able to drag the
    junction on its own.
  * that refinement touches coordinates and nothing else. It runs after
    reconciliation precisely because moving a feature earlier turned a boundary
    change into an inventory change once already.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plastanno.core.feature import Feature                       # noqa: E402
from plastanno.identify import intron_refine as IR               # noqa: E402

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %s %s" % ("ok  " if ok else "FAIL", label))
    if not ok:
        FAIL.append("%s: got %r want %r" % (label, got, want))


# ── name mapping ──────────────────────────────────────────────────────────────
print("gene names")
check("a full name passes through", IR.canonical_gene("trnA-UGC"), "trnA-UGC")
check("a bare name resolves", IR.canonical_gene("trnI"), "trnI-GAU")
# trnG-UCC carries the intron and trnG-GCC mostly does not, so a bare 'trnG'
# cannot be resolved and must be left alone rather than guessed.
check("bare trnG is NOT guessed", IR.canonical_gene("trnG"), None)
check("a non-intron tRNA is rejected", IR.canonical_gene("trnF-GAA"), None)
check("empty input is rejected", IR.canonical_gene(""), None)
check("None is rejected", IR.canonical_gene(None), None)


# ── glocal orientation ────────────────────────────────────────────────────────
print("\nglocal placement")
_t = "AAAAAAAAAA" + "GGTTCGATTC" + "A" * 20
p = IR._place("GGTTCGATTC", _t)
check("an exact query lands on its own coordinates", (p[1], p[2]), (10, 20))
check("  ... at the full match score", p[0], 20.0)
# The failure mode that matters is not an exception, it is a plausible-looking
# wrong answer. A swapped configuration scores this same input negative.
check("  ... with a positive score, i.e. end gaps are free on the GENOME",
      p[0] > 0, True)
check("a query longer than the window is refused",
      IR._place("A" * 50, "ACGT"), None)
check("an empty query is refused", IR._place("", _t), None)


# ── the vote ──────────────────────────────────────────────────────────────────
print("\ndonor vote")
# Seven placements agree; one outlier sits five bases to the right on the
# junction and carries the best score. Raw "take the best" would follow it.
agree = [(1.0, (100, 140), (900, 940))] * 7
outlier = [(9.9, (100, 145), (905, 940))]
v = IR._vote(agree + outlier)
check("the majority junction wins over a better-scoring outlier",
      v, [(100, 140), (900, 940)])
check("an empty placement list yields nothing", IR._vote([]), None)
# Coordinates that do not form an ordered exon pair are not a boundary at all.
check("a non-monotonic vote is rejected",
      IR._vote([(1.0, (900, 940), (100, 140))] * 5), None)


# ── refinement leaves everything but coordinates alone ────────────────────────
print("\nrefine_intron_trna: what it may and may not touch")

E1 = "GGGCTATTAGCTCAGTGGTAGAGCGCGCCCCT"
E2 = "TTCACGGGCGAGGTCTCTGGTTCAAGTCCAGGATGGCCCA"
INTRON = "GTGCG" + "TTACGATCGATTGCA" * 50 + "TTCAC"          # 760 bp
GENOME = "TTGCA" * 120 + E1 + INTRON + E2 + "ACGTT" * 120
TRUE_S = 600
TRUE_E1 = (TRUE_S, TRUE_S + len(E1))
TRUE_E2 = (TRUE_S + len(E1) + len(INTRON),
           TRUE_S + len(E1) + len(INTRON) + len(E2))

DONORS = {"trnI-GAU": [("NC_000001", E1, E2)] * 5}
IR._DONOR_CACHE["__test__"] = DONORS


def mk(name, exons, strand=1):
    f = Feature(gene_name=name, gene_type="tRNA")
    f.exons = list(exons)
    f.start, f.end, f.strand = exons[0][0], exons[-1][1], strand
    return f


# A call that is three bases off at both ends of the junction.
off = mk("trnI-GAU", [(TRUE_E1[0], TRUE_E1[1] + 3), (TRUE_E2[0] - 3, TRUE_E2[1])])
single = mk("trnI-GAU", [(50, 122)])
wrong_gene = mk("trnF-GAA", [(TRUE_E1[0], TRUE_E1[1] + 3),
                             (TRUE_E2[0] - 3, TRUE_E2[1])])
wrapped = mk("trnA-UGC", [(len(GENOME) - 40, len(GENOME)), (0, 35)])
wrapped.start, wrapped.end = len(GENOME) - 40, 35

feats = [off, single, wrong_gene, wrapped]
before = [(f.gene_name, f.gene_type, f.strand) for f in feats]
diag = []
n = IR.refine_intron_trna(feats, GENOME, "__test__", len(GENOME), diagnostics=diag)

check("the off-by-three call is re-placed", n, 1)
check("  ... onto the true exon pair", off.exons, [TRUE_E1, TRUE_E2])
check("  ... and start/end follow the exons",
      (off.start, off.end), (TRUE_E1[0], TRUE_E2[1]))
check("  ... recording where it came from",
      any(x == "geometry_source=glocal-exon-vote" for x in off.notes), True)
check("  ... and what it replaced",
      any(x.startswith("geometry_was=") for x in off.notes), True)

check("a single-exon tRNA is untouched", single.exons, [(50, 122)])
check("a tRNA with no group II intron is untouched",
      wrong_gene.exons, [(TRUE_E1[0], TRUE_E1[1] + 3), (TRUE_E2[0] - 3, TRUE_E2[1])])
check("an origin-crossing locus is left alone, not corrupted",
      (wrapped.start, wrapped.end), (len(GENOME) - 40, 35))
check("  ... and is reported rather than dropped silently",
      any("origin" in d for d in diag), True)

check("the list is not appended to or removed from", len(feats), 4)
check("no identity, type or strand changed",
      [(f.gene_name, f.gene_type, f.strand) for f in feats], before)

# Idempotence: a second pass over an already-correct call must be a no-op, or
# repeated runs would drift.
n2 = IR.refine_intron_trna([off], GENOME, "__test__", len(GENOME))
check("re-running on a correct call changes nothing", n2, 0)
check("  ... and leaves the coordinates alone", off.exons, [TRUE_E1, TRUE_E2])

# A missing donor set must degrade to doing nothing, not to an exception.
IR._DONOR_CACHE["__empty__"] = {}
d2 = []
check("an empty donor DB is a no-op",
      IR.refine_intron_trna([mk("trnI-GAU", [(600, 640), (1400, 1440)])],
                            GENOME, "__empty__", len(GENOME), diagnostics=d2), 0)
check("  ... and says so", any("donor" in x for x in d2), True)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
