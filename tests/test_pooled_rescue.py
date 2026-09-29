#!/usr/bin/env python3
"""The pooled rule must not lose a locus that Engine A found completely.

`pooled_candidates` keeps one candidate of each AB or conflict pair. `_select`
applies its length filter AFTER the per-cluster pick, so when the kept candidate
is a fragment it is deleted and the locus vanishes although a complete call
existed. Measured on 30 development genomes drawn from splits/dev_set.txt: 6 of
2222 discarded Engine A calls left their gene absent from the output, 5 of them
ndhA. The ndhA case, e.g. NC_053736.1: Engine A 1107 bp spliced over two exons,
Engine B a 561 bp one-exon fragment; B was kept, its validate_orf score being no
lower than A's, and B then fell under the 0.6 x 1092 bp length filter.

The fix leaves `_select`'s input untouched and only afterwards rescues a
discarded candidate whose preferred partner did not survive, accepting it only
if `_select` itself — run on copies — keeps it without changing any survivor.

The fixture reproduces the ndhA mechanism exactly: A is full length but lacks a
start codon (0.75), B is a 55% fragment with a clean start (0.75), they tie, B
is preferred, and B falls under the 0.6x filter.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.core.feature import Feature                    # noqa: E402
from plastanno.core import reconcile as RC                    # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-66s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


def build_genome(n=8000):
    import random
    rng = random.Random(7)
    g = list("".join(rng.choice("ACG") for _ in range(n)))  # no T: no stray stops

    def put(pos, s):
        g[pos:pos + len(s)] = list(s)
    # gX at 1000: 600 bp, NO start codon, one in-frame ATG at +270, stop at end
    gx = "GCT" + "GCT" * 89 + "ATG" + "GCT" * 108 + "TAA"
    assert len(gx) == 600
    put(1000, gx)
    # gY at 3000: a clean 300 bp ORF, the control
    gy = "ATG" + "GCT" * 98 + "TAA"
    put(3000, gy)
    # gZ at 5000: 600 bp with no start, and a 150 bp fragment — both too short
    # to survive when A is also made short below
    gz = "GCT" + "GCT" * 189 + "ATG" + "GCT" * 8 + "TAA"
    assert len(gz) == 600
    put(5000, gz)
    # gW at 6500: 600 bp, clean start and stop, but an in-frame TAA at +300 -- a
    # pseudogene remnant; the B fragment after the stop is clean and too short
    gw = "ATG" + "GCT" * 99 + "TAA" + "GCT" + "ATG" + "GCT" * 96 + "TAA"
    assert len(gw) == 600, len(gw)
    put(6500, gw)
    return "".join(g)


GENOME = build_genome()
CAT = {"gX": {"expected_len": 600, "n_exons": 1},
       "gY": {"expected_len": 300, "n_exons": 1},
       "gZ": {"expected_len": 600, "n_exons": 1},
       "gW": {"expected_len": 600, "n_exons": 1}}


def F(name, s, e, engine, **kw):
    f = Feature(gene_name=name, gene_type="CDS", start=s, end=e, strand=1,
                engine=engine, **kw)
    f.exons = [(s, e)]
    return f


def keys(feats):
    return sorted(RC._feature_key(f) for f in feats)


# ── the fixture is what the test assumes ─────────────────────────────────────
print("--- fixture ---")
A_x = F("gX", 1000, 1600, "A", s_ref=0.93)
B_x = F("gX", 1270, 1600, "B", s_model=0.90)
check("A is full length but lacks a start: validate_orf 0.75",
      RC.validate_orf(A_x, GENOME, CAT), 0.75)
check("B is a 330 bp fragment with a clean start: validate_orf 0.75",
      RC.validate_orf(B_x, GENOME, CAT), 0.75)
check("they are paired (AB or conflict), not left single",
      len(RC.match_features([A_x], [B_x])["AB"])
      + len(RC.match_features([A_x], [B_x])["conflict"]), 1)
check("B is under the 0.6x filter (330/600 = 0.55)", 330 / 600 < 0.6, True)

# ── 1. the ndhA mechanism: the locus is recovered ────────────────────────────
print("--- 1. B preferred, B filtered out, A complete -> A is rescued ---")
out1 = RC.reconcile_pooled([F("gX", 1000, 1600, "A", s_ref=0.93)],
                           [F("gX", 1270, 1600, "B", s_model=0.90)],
                           GENOME, CAT)
gx1 = [f for f in out1 if f.gene_name == "gX"]
check("exactly one gX feature comes out", len(gx1), 1)
check("  ... and it is Engine A's full-length call",
      (gx1[0].start, gx1[0].end) if gx1 else None, (1000, 1600))
check("  ... carrying the rescue note",
      any("rescued" in n for n in gx1[0].notes) if gx1 else False, True)

# ── 2. the ordinary case must not change ─────────────────────────────────────
print("--- 2. preferred candidate survives -> nothing is rescued ---")
out2 = RC.reconcile_pooled([F("gY", 3000, 3300, "A", s_ref=0.9)],
                           [F("gY", 3000, 3300, "B", s_model=0.9)],
                           GENOME, CAT)
gy2 = [f for f in out2 if f.gene_name == "gY"]
check("exactly one gY feature, as before", len(gy2), 1)
check("  ... with no rescue note", any("rescued" in n for n in gy2[0].notes), False)

# ── 3. a fallback that would itself be filtered is not rescued ───────────────
print("--- 3. fallback below the length filter -> not rescued ---")
out3 = RC.reconcile_pooled([F("gZ", 5000, 5200, "A", s_ref=0.9)],
                           [F("gZ", 5000, 5150, "B", s_model=0.9)],
                           GENOME, CAT)
check("nothing is emitted for gZ", [f for f in out3 if f.gene_name == "gZ"], [])

# ── 4. _rescue refuses a fallback that overlaps a same-gene survivor ─────────
print("--- 4. fallback overlapping a surviving same-gene call -> refused ---")
# The SURVIVOR must be the better call here, so that inside the probe it wins the
# cluster and _select writes "selected over ..." onto it. That is the write a
# probe on real objects would leak into the output. A first version of this case
# let the fallback win, so the note landed on the fallback and a probe without
# deepcopy passed unnoticed.
kept = F("gX", 1270, 1600, "B", s_model=0.9)
fb = F("gX", 1000, 1480, "A", s_ref=0.9)        # 480/600: passes, but farther
survivor = F("gX", 1000, 1600, "A", s_ref=0.9)  # 600/600: wins the probe cluster
for f in (kept, fb, survivor):
    RC._score_pooled(f, GENOME, CAT)
notes_before = list(survivor.notes)
res4 = RC._rescue([kept], [survivor], {id(kept): fb}, GENOME, CAT, None)
check("no rescue: the locus is already covered", res4, [])
check("the survivor's notes are untouched by the probe",
      survivor.notes, notes_before)

# ── 5. a survivor of the preferred candidate means no rescue at all ──────────
print("--- 5. the preferred candidate is alive -> its fallback is never tried ---")
# "Never tried", not merely "not rescued": the probe would reject such a fallback
# anyway, since a pair always overlaps, so an output check cannot tell the two
# apart. What the liveness check buys is not probing ~74 fallbacks per genome.
# Count the selector calls instead.
kept5 = F("gY", 3000, 3300, "B", s_model=0.9)
fb5 = F("gY", 3000, 3300, "A", s_ref=0.9)
RC._score_pooled(kept5, GENOME, CAT)
_calls = [0]
_real_sel = RC._select_pooled
def _counting(*a, **k):
    _calls[0] += 1
    return _real_sel(*a, **k)
RC._select_pooled = _counting
try:
    res5 = RC._rescue([kept5], [kept5], {id(kept5): fb5}, GENOME, CAT, None)
finally:
    RC._select_pooled = _real_sel
check("nothing rescued", res5, [])
check("  ... and the selector was never invoked for it", _calls[0], 0)

# ── 6. pooled_candidates is unchanged when no fallback dict is passed ────────
print("--- 6. pooled_candidates: identical output with and without the dict ---")
a6 = [F("gX", 1000, 1600, "A", s_ref=0.93), F("gY", 3000, 3300, "A", s_ref=0.9)]
b6 = [F("gX", 1270, 1600, "B", s_model=0.9), F("gY", 3000, 3300, "B", s_model=0.9)]
plain = RC.pooled_candidates(a6, b6, GENOME, CAT)
a6b = [F("gX", 1000, 1600, "A", s_ref=0.93), F("gY", 3000, 3300, "A", s_ref=0.9)]
b6b = [F("gX", 1270, 1600, "B", s_model=0.9), F("gY", 3000, 3300, "B", s_model=0.9)]
fbd = {}
withd = RC.pooled_candidates(a6b, b6b, GENOME, CAT, fallbacks=fbd)
check("same features returned", keys(plain), keys(withd))
check("one fallback recorded per pair", len(fbd), 2)

# ── 7. the rescue changes nothing else in a mixed genome ─────────────────────
print("--- 7. mixed: gX rescued, gY exactly as it would be alone ---")
out7 = RC.reconcile_pooled(
    [F("gX", 1000, 1600, "A", s_ref=0.93), F("gY", 3000, 3300, "A", s_ref=0.9)],
    [F("gX", 1270, 1600, "B", s_model=0.9), F("gY", 3000, 3300, "B", s_model=0.9)],
    GENOME, CAT)
check("gY is emitted exactly as in case 2",
      keys([f for f in out7 if f.gene_name == "gY"]), keys(gy2))
check("gX is the rescued full-length call",
      keys([f for f in out7 if f.gene_name == "gX"]), keys(gx1))

# ── 8. the gate: after the final QC, a rescued call must be a plausible gene ──
print("--- 8. gate: revoke an implausible RESCUED call, and nothing else ---")
aw, bw = F("gW", 6500, 7100, "A", s_ref=0.93), F("gW", 6806, 7100, "B", s_model=0.9)
check("fixture: A carries one in-frame internal stop",
      RC._internal_stops(RC._coords.extract(GENOME, aw, len(GENOME))), 1)
out8 = RC.reconcile_pooled([aw], [bw], GENOME, CAT)
gw8 = [f for f in out8 if f.gene_name == "gW"]
check("reconcile_pooled still rescues it (the gate comes later)", len(gw8), 1)
kept, revoked = RC.revoke_implausible_rescues(out8, GENOME)
check("the final gate revokes the internal-stop rescue",
      ([f.gene_name for f in revoked], [f for f in kept if f.gene_name == "gW"]), (["gW"], []))
low = RC.reconcile_pooled([F("gX", 1000, 1600, "A", s_ref=0.01)],
                          [F("gX", 1270, 1600, "B", s_model=0.9)], GENOME, CAT)
low_gx = [f for f in low if f.gene_name == "gX"]
check("a low-confidence fallback is rescued as NEEDS_REVIEW",
      [f.flag for f in low_gx], ["NEEDS_REVIEW"])
k2, r2 = RC.revoke_implausible_rescues(low, GENOME)
check("  ... and revoked", [f.gene_name for f in r2], ["gX"])
k1, r1 = RC.revoke_implausible_rescues(out1, GENOME)
check("a clean rescue is kept (case 1)", (len(r1), [(f.start, f.end) for f in k1 if f.gene_name == "gX"]),
      (0, [(1000, 1600)]))
plain = F("gW", 6500, 7100, "A", s_ref=0.9)
plain.flag = "NEEDS_REVIEW"
k3, r3 = RC.revoke_implausible_rescues([plain], GENOME)
check("a feature the rescue did NOT add is never revoked, whatever its flag", (len(k3), r3), (1, []))
check("validate_orf still counts internal stops the same way",
      RC.validate_orf(F("gW", 6500, 7100, "A"), GENOME, CAT), 0.75)
import inspect                                                  # noqa: E402
import plastanno.pipeline as PL                                 # noqa: E402
src = inspect.getsource(PL.run)
i_qc, i_rv, i_ap, i_w = (src.find("finalize_qc(annotations"), src.find("revoke_implausible_rescues(annotations"),
                         src.find("assign_cds_products("), src.find("write_all("))
check("pipeline.run revokes after finalize_qc and before products and writing",
      0 <= i_qc < i_rv < i_ap < i_w, True)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
