#!/usr/bin/env python3
"""The IR-lacking plastome, which nothing in this suite used to exercise.

About one plastome in ten has no inverted repeat — conifers, many bryophytes,
some lycophytes have genuinely lost one copy. Measured here: 4 of the 30
development genomes and 9 of the 99 holdout genomes.

Every defect found in the 2026-09-23 audit lived on that path:

  * rps12 exon1 was chosen as "the LSC fragment nearest the LSC start", but with
    no IR the pipeline sets LSC = the whole genome, so the tiebreak degenerated to
    lowest coordinate and picked the 3' half of the gene;
  * the rps12 reconstruction loop iterated over IRb and IRa and therefore never
    executed at all, reporting a failure it had not attempted;
  * the whole-genome anchor search returns exon1 as well as exon2, and on
    NC_039155.1 exon1 scored higher, so the anchor landed on it.

None of the three could have survived a single IR-lacking genome in the tests.
This file is that genome, at the level the defects actually lived: the selection
logic, exercised through the real functions with Exonerate stubbed out, so it runs
in milliseconds and always runs.

Real IR-lacking accessions in the development set, for whoever builds a sample:
NC_037507.1 (Marchantiaceae), NC_039155.1 (Araucariaceae), NC_039564.1,
NC_079849.1 (Frullaniaceae).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plastanno.core.feature import Feature              # noqa: E402
from plastanno.identify import ir_detector as IRD       # noqa: E402
import plastanno.annotate.special_cases as SC           # noqa: E402
import plastanno.identify.engine_a as EA                # noqa: E402

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-68s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


GLEN = 120_000
NO_IR = {"LSC": (0, GLEN)}          # what pipeline.run substitutes when detection fails
WITH_IR = {"LSC": (0, 60_000), "IRb": (60_000, 80_000),
           "SSC": (80_000, 90_000), "IRa": (90_000, 110_000)}


print("--- 1. the contract: detection returning None means a single LSC ---")
# pipeline.run does: ir_boundaries = detect(...) or {"LSC": (0, genome_len)}
fallback = None or {"LSC": (0, GLEN)}
check("the fallback names exactly one region", sorted(fallback), ["LSC"])
check("and it spans the whole genome", fallback["LSC"], (0, GLEN))
check("so IRb and IRa are absent, not empty",
      ("IRb" in fallback, "IRa" in fallback), (False, False))

print("--- 2. Engine A falls back to the whole genome, it does not skip ---")
# An IR gene with no IR present must still be searched, over everything.
regs = EA.get_search_regions("ycf2", NO_IR, GLEN)
check("an IR gene still gets a region to search", len(regs) >= 1, True)
check("and that region is the whole genome", regs, [(0, GLEN)])
with_ir = EA.get_search_regions("ycf2", WITH_IR, GLEN)
check("with a real IR it searches both copies", len(with_ir), 2)

print("--- 3/4/5. handle_rps12 on an IR-lacking genome, through the real path ---")
# Two rps12 fragments, both "in the LSC" because with no IR the LSC is everything.
# The 3' block sits at the LOWER coordinate, exactly as in NC_037507.1 where it was
# at 92..851 against a real exon1 at 66686..66800.
#
# _reconstruct_rps12_copy and the Exonerate anchor search are stubbed so this runs
# in milliseconds; what is under test is the SELECTION, which is where all three
# defects lived. Stubbing the reconstruction is also what makes defect 2 visible:
# if the loop never runs, the stub is never called.
GEN = list("ACGT" * (GLEN // 4))
def _put(seq, at):
    for i, ch in enumerate(seq):
        GEN[at + i] = ch
_put("ATG" + "GCT" * 37, 66_686)          # a Met-starting 114 nt exon1
_put("ACG" + "CCT" * 37, 92)              # a 3' block that does NOT start with Met
GENOME = "".join(GEN)

frag_lo = Feature(gene_name="rps12", gene_type="CDS", start=92, end=92 + 114,
                  strand=1, exons=[(92, 92 + 114)], engine="A", s_ref=0.9)
frag_hi = Feature(gene_name="rps12", gene_type="CDS", start=66_686, end=66_800,
                  strand=1, exons=[(66_686, 66_800)], engine="A", s_ref=0.9)

calls = []
real_recon, real_region = SC._reconstruct_rps12_copy, None
import plastanno.identify.engine_a as _EA_mod

def fake_recon(genome, exon1, ir_strand, exon2_start_g, refs):
    calls.append({"exon1": exon1, "anchor": exon2_start_g})
    return None                            # selection is what is under test

def fake_region(genome_seq, protein_seq, gene_name, prot_acc,
                region_start, region_end, genome_len=None):
    # Return BOTH blocks, as a whole-genome search does. exon1 scores higher,
    # which is what made the anchor land on it in NC_039155.1.
    hi = Feature(gene_name="rps12", gene_type="CDS", start=66_686, end=66_800,
                 strand=1, exons=[(66_686, 66_800)], engine="A", s_ref=0.895)
    lo = Feature(gene_name="rps12", gene_type="CDS", start=92, end=92 + 301,
                 strand=1, exons=[(92, 92 + 301)], engine="A", s_ref=0.880)
    return [hi, lo]

SC._reconstruct_rps12_copy = fake_recon
_EA_mod.run_exonerate_region = fake_region
try:
    SC.handle_rps12([frag_lo, frag_hi], GENOME, NO_IR, protein_db="database/protein_db")
finally:
    SC._reconstruct_rps12_copy = real_recon

check("the reconstruction IS attempted with no IR (defect 2)", len(calls) > 0, True)
if calls:
    check("exon1 is the Met-starting fragment, not the lowest coordinate (defect 1)",
          calls[0]["exon1"][0], 66_686)
    e1s, e1e = calls[0]["exon1"][0], calls[0]["exon1"][1]
    check("the anchor is NOT inside exon1 (defect 3)",
          not (e1s <= calls[0]["anchor"] < e1e), True)
    check("the anchor is the 3' block", calls[0]["anchor"], 92)

print("--- 5b. an IR question is ANSWERED with no IR, not abandoned ---")
# The distinction that made one of three identical-looking guards a bug. This one
# is correct: "mirror this span into the other repeat" has a right answer when
# there is no other repeat, and it returns it.
check("mirroring into the other repeat is None when there is no repeat",
      SC._mirror_into_other_ir((1000, 2000), NO_IR, GLEN), None)
check("and it still works when a real IR is present",
      SC._mirror_into_other_ir((61_000, 62_000), WITH_IR, GLEN) is not None, True)

print("--- 6. the detector reports absence rather than inventing a repeat ---")
# A sequence with no large inverted repeat must return None, not a spurious pair.
check("detect_ir_boundaries returns None on a non-repetitive sequence",
      IRD.detect_ir_boundaries("ACGT" * 500), None)


print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
