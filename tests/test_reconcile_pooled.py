#!/usr/bin/env python3
"""Unit tests for the pooled (B-primary / A-rescue) reconciliation path.

Run DIRECTLY and check $?. Exit 0 = all pass.

These pin the five behaviours the architecture depends on, plus the guarantee
that the legacy path is untouched.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.core.feature import Feature
from plastanno.core.reconcile import (match_features, reconcile, reconcile_pooled,
                                      pooled_candidates, validate_orf)
from plastanno.core.ir_genes import (IR_DUPLICATED_GENES, LEGACY_IR_CONFLICT_GENES,
                                     BUILD_CATALOG_IR_GENES, LEGACY_IR_CONFLICT_GAP)

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-64s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


# ── a genome with a clean ORF at 1000 and a stopless one at 2000 ─────────────
def build_genome(n=6000):
    import random
    rng = random.Random(11)
    # avoid accidental in-frame stops in the filler
    g = list("".join(rng.choice("ACG") for _ in range(n)))
    def put(pos, s):
        g[pos:pos + len(s)] = list(s)
    clean = "ATG" + "".join("GCT" for _ in range(98)) + "TAA"      # 300 bp, clean
    put(1000, clean)
    broken = "ATG" + "".join("GCT" for _ in range(98)) + "GCT"     # 300 bp, no stop
    put(2000, broken)
    return "".join(g)


GENOME = build_genome()
CAT = {"ndhB": {"expected_len": 300, "n_exons": 1},
       "psbA": {"expected_len": 300, "n_exons": 1}}


def F(name, s, e, strand=1, engine="A", **kw):
    f = Feature(gene_name=name, gene_type="CDS", start=s, end=e, strand=strand,
                engine=engine, **kw)
    f.exons = [(s, e)]
    return f


print("--- the ORF fixtures are what the tests assume ---")
check("the ORF at 1000 is clean (validate_orf == 1.0)",
      validate_orf(F("ndhB", 1000, 1300), GENOME, CAT), 1.0)
check("the ORF at 2000 is not (no terminal stop)",
      validate_orf(F("ndhB", 2000, 2300), GENOME, CAT) < 1.0, True)

print("--- 1. the conflict bin holds 3-tuples, not pairs ---")
# same strand, partial overlap below the 0.5 containment threshold
A = [F("ndhB", 1000, 2000)]
B = [F("ndhB", 1700, 2700, engine="B", s_model=0.9)]
bins = match_features(A, B)
check("one conflict entry", len(bins["conflict"]), 1)
check("it is a 3-tuple (fa, fb, reason)", len(bins["conflict"][0]), 3)
check("the third element is a human-readable reason",
      isinstance(bins["conflict"][0][2], str), True)

print("--- 2. a conflict is ONE locus: B wins, whatever the gene is called ---")
# match_features only pairs calls that overlap, so a conflict is by construction
# two predictions of the same physical locus -- never two copies. There is no
# gene-name special case: an "IR genes" list is a POLICY, not biology, and IR
# content varies between species with repeat expansion and contraction.
CATR = dict(CAT, rrn16={"expected_len": 300, "n_exons": 1})
for gene, cat in (("psbA", CAT), ("rrn16", CATR), ("ndhB", CAT)):
    cand = pooled_candidates([F(gene, 1000, 2000)],
                             [F(gene, 1700, 2700, engine="B", s_model=0.9)],
                             GENOME, cat)
    check("%-6s conflict yields ONE candidate, not two" % gene, len(cand), 1)
check("an IR-listed gene is treated no differently from any other",
      "rrn16" in IR_DUPLICATED_GENES, True)

print("--- 3. a same-locus conflict is resolved B-primary, A only on ORF failure ---")
# A holds the clean ORF (1000-1300), B a span that is not one -> A substitutes
A2 = [F("psbA", 1000, 1300)]
B2 = [F("psbA", 1200, 2200, engine="B", s_model=0.9)]
bins2 = match_features(A2, B2)
check("these land in the conflict bin", len(bins2["conflict"]), 1)
cand2 = pooled_candidates(list(A2), list(B2), GENOME, CAT)
check("one candidate is handed over", len(cand2), 1)
check("and it is the one with the valid ORF",
      (cand2[0].start, cand2[0].end), (1000, 1300))
check("the decision is recorded as a same-locus conflict",
      any("same-locus conflict" in n for n in cand2[0].notes), True)
# and the reverse: when B's model is the clean one, B wins
A2b = [F("psbA", 1200, 2200)]
B2b = [F("psbA", 1000, 1300, engine="B", s_model=0.9)]
cand2b = pooled_candidates(A2b, B2b, GENOME, CAT)
check("when Engine B holds the valid ORF, B wins",
      (cand2b[0].start, cand2b[0].end), (1000, 1300))

print("--- 4. B is primary where both engines agree on the locus ---")
A3 = [F("psbA", 1000, 1300, s_ref=0.95)]
B3 = [F("psbA", 1000, 1300, engine="B", s_model=0.90)]
check("this is an AB pair", len(match_features(A3, B3)["AB"]), 1)
out3 = reconcile_pooled(list(A3), list(B3), GENOME, CAT)
check("one feature results", len(out3), 1)
check("Engine B's model is kept when it passes the ORF check",
      any("engine B model kept" in n for n in out3[0].notes), True)
check("no cross-engine agreement term is credited", out3[0].s_overlap, 0.0)
# a feature arriving with a non-zero agreement term must have it cleared --
# with a fresh Feature the field is already 0.0, so the clearing is invisible
A3b = [F("psbA", 1000, 1300, s_ref=0.95, s_overlap=0.77)]
B3b = [F("psbA", 1000, 1300, engine="B", s_model=0.90, s_overlap=0.77)]
out3b = reconcile_pooled(A3b, B3b, GENOME, CAT)
check("an incoming agreement term is CLEARED, not carried through",
      out3b[0].s_overlap, 0.0)

print("--- 5. an Engine-A-only locus is rescued, not dropped ---")
A4 = [F("psbA", 1000, 1300, s_ref=0.95)]
out4 = reconcile_pooled(list(A4), [], GENOME, CAT)
check("the A-only locus survives", len(out4), 1)
check("and is marked as a rescue",
      any("reference-only rescue" in n for n in out4[0].notes), True)

print("--- 6. two real IR copies are kept on GEOMETRY, and survive a shift ---")
# Each copy has its own A and B call, disagreeing slightly on the boundary, so
# every copy exercises the merge independently. Copy 1 at 1000 (+), copy 2 at
# 4000 (-). They never overlap, so they are two loci, not a conflict.
def two_copies(shift=0):
    g = GENOME[shift:] + GENOME[:shift] if shift else GENOME
    m = lambda x: (x - shift) % len(GENOME)
    A = [F("ndhB", m(1000), m(1300), strand=1, s_ref=0.95),
         F("ndhB", m(4000), m(4300), strand=-1, s_ref=0.95)]
    B = [F("ndhB", m(1002), m(1300), engine="B", strand=1, s_model=0.9),
         F("ndhB", m(4000), m(4298), engine="B", strand=-1, s_model=0.9)]
    return g, A, B

g0, A5, B5 = two_copies()
b5 = match_features(A5, B5)
check("each copy pairs with its own partner, none conflict",
      (len(b5["AB"]), len(b5["conflict"])), (2, 0))
out5 = reconcile_pooled(A5, B5, g0, CAT)
check("exactly two ndhB loci survive",
      sum(1 for f in out5 if f.gene_name == "ndhB"), 2)
check("at two distinct positions",
      len({(f.start, f.end) for f in out5 if f.gene_name == "ndhB"}), 2)
check("on opposite strands",
      sorted(f.strand for f in out5 if f.gene_name == "ndhB"), [-1, 1])

# the same construction shifted along the molecule must give the same answer.
# (A full rotation matrix, including origin-crossing, is the Stage 2 gate; this
# pins that the merge itself carries no absolute-position dependence.)
for shift in (500, 2500):
    gs, As, Bs = two_copies(shift)
    outs = reconcile_pooled(As, Bs, gs, CAT)
    check("shift %-5d still yields exactly two ndhB loci" % shift,
          sum(1 for f in outs if f.gene_name == "ndhB"), 2)

print("--- 7. the legacy path is untouched ---")
leg = reconcile(engine_a=[F("psbA", 1000, 1300, s_ref=0.95)],
                engine_b=[F("psbA", 1000, 1300, engine="B", s_model=0.90)],
                genome_seq=GENOME, gene_catalog=CAT)
check("legacy still produces a feature", len(leg), 1)
check("legacy still credits cross-engine agreement", leg[0].s_overlap > 0, True)
check("legacy still marks it confirmed by both engines",
      any("confirmed by both" in n for n in leg[0].notes), True)
check("pooled and legacy differ on exactly that term",
      out3[0].s_overlap != leg[0].s_overlap, True)

print("--- 8. one owner for the IR gene sets ---")
check("the legacy subset is strictly smaller",
      LEGACY_IR_CONFLICT_GENES < IR_DUPLICATED_GENES, True)
check("the documented gap is exactly the 8 known genes",
      sorted(LEGACY_IR_CONFLICT_GAP),
      ["rrn16", "rrn23", "rrn4.5", "rrn5", "trnL-CAA", "trnN-GUU", "trnR-ACG",
       "trnV-GAC"])
check("the build-catalog variant differs only by orf70",
      BUILD_CATALOG_IR_GENES - IR_DUPLICATED_GENES, {"orf70"})
check("the pooled path uses the COMPLETE set",
      all(g in IR_DUPLICATED_GENES for g in ("rrn16", "trnN-GUU")), True)

print()
print("checks run: %d" % RUN[0])
if FAIL:
    print("%d FAILED:" % len(FAIL))
    for f in FAIL:
        print("   - %s" % f)
    sys.exit(1)
print("all %d checks passed" % RUN[0])
sys.exit(0)
