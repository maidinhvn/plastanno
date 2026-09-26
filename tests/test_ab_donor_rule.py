#!/usr/bin/env python3
"""Which engine donates coordinates when both find the same CDS.

`reconcile()` merges an A/B pair into one feature and has to take the coordinates
from one of them. The rule used to be an absolute floor: hand over to Engine B
only when A's length fell below 0.6 of `expected_len`. atpA in NC_039155.1 failed
at 0.673 -- Engine A had 1026 bp of a 1524 bp gene, Engine B had 1521, and A
donated anyway. The merged feature then lost a `_select` paralog cluster to a
spurious atpB and the gene was reported as absent. See
benchmark_v3/cds/ATPA_TRACE.md and benchmark_v3/donor_rule/RESULT.md.

The rule is now a comparison: whichever length is closer to `expected_len`.

These tests pin the comparison, the case that motivated it, the old floor case
that must keep working, and the boundaries of the branch -- it is CDS-only and it
needs an `expected_len` to consult.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plastanno.core.feature import Feature        # noqa: E402
from plastanno.core.reconcile import reconcile    # noqa: E402

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-66s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


GLEN = 20_000
# A genome of real bases. Content does not steer the donor -- only the lengths do
# -- but validate_orf reads it, so it must exist and be long enough.
GENOME = ("ATG" + "GCT" * 3000 + "TAA") * 2
GENOME = (GENOME + "A" * GLEN)[:GLEN]


def cat(gene, exp, gtype="CDS"):
    return {gene: {"type": gtype, "region": "LSC", "n_exons": 1,
                   "synonym_of": None, "product": gene, "expected_len": exp}}


def pair(gene, a_span, b_span, gtype="CDS"):
    """One Engine A and one Engine B call of the same gene, overlapping."""
    fa = Feature(gene_name=gene, gene_type=gtype, start=a_span[0], end=a_span[1],
                 strand=1, exons=[a_span], engine="A", s_ref=0.95)
    fb = Feature(gene_name=gene, gene_type=gtype, start=b_span[0], end=b_span[1],
                 strand=1, exons=[b_span], engine="B", s_model=0.95)
    return fa, fb


def donor_span(gene, a_span, b_span, exp, gtype="CDS"):
    """Run the merge and report which call's coordinates came out."""
    fa, fb = pair(gene, a_span, b_span, gtype)
    out = reconcile(engine_a=[fa], engine_b=[fb], genome_seq=GENOME,
                    gene_catalog=cat(gene, exp, gtype),
                    ir_boundaries={"LSC": (0, GLEN)})
    if not out:
        return "DROPPED"
    f = out[0]
    got = (f.start, f.end)
    return "A" if got == tuple(a_span) else "B" if got == tuple(b_span) else got


print("--- 1. the case that motivated the change ---")
# atpA: A has 1026 bp, B has 1521, expected 1524. |1521-1524|=3 beats |1026-1524|=498.
check("B donates when it is far closer to expected",
      donor_span("atpA", (2000, 3026), (1500, 3021), 1524), "B")

print("--- 2. the floor case the old rule handled must still work ---")
# A catastrophically short: 300 of 1524. The old rule caught this and so must this one.
check("B donates when A is a fragment",
      donor_span("atpA", (2000, 2300), (1500, 3021), 1524), "B")

print("--- 3. Engine A still donates when it is the closer one ---")
check("A donates when A is closer to expected",
      donor_span("atpA", (1500, 3021), (1500, 2800), 1524), "A")
check("A donates when both are equal distance (no gratuitous change)",
      donor_span("rbcL", (1000, 2400), (1000, 2600), 1500), "A")

print("--- 4. the branch is CDS-only ---")
# A tRNA pair with the same length relationship must NOT consult expected_len.
check("a non-CDS pair keeps Engine A's coordinates",
      donor_span("trnK-UUU", (1000, 1040), (1000, 1073), 72, gtype="tRNA"), "A")

print("--- 5. no expected_len means no comparison to make ---")
fa, fb = pair("mystery", (2000, 3026), (1500, 3021))
out = reconcile(engine_a=[fa], engine_b=[fb], genome_seq=GENOME,
                gene_catalog={"mystery": {"type": "CDS", "region": "LSC",
                                          "n_exons": 1, "synonym_of": None,
                                          "product": "x"}},
                ir_boundaries={"LSC": (0, GLEN)})
check("Engine A donates when the catalog has no expected_len",
      (out[0].start, out[0].end) if out else "DROPPED", (2000, 3026))

print("--- 6. the merged feature is marked as confirmed by both engines ---")
fa, fb = pair("atpA", (2000, 3026), (1500, 3021))
out = reconcile(engine_a=[fa], engine_b=[fb], genome_seq=GENOME,
                gene_catalog=cat("atpA", 1524), ir_boundaries={"LSC": (0, GLEN)})
check("engine is AB", out[0].engine if out else None, "AB")
check("s_ref comes from A", round(out[0].s_ref, 2) if out else None, 0.95)
check("s_model comes from B", round(out[0].s_model, 2) if out else None, 0.95)


print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
