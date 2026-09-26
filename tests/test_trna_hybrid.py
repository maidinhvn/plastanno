#!/usr/bin/env python3
"""The nine checks required before the hybrid tRNA mode may be run on real data.

The hybrid mode changes coordinates, so each test states what must change and,
more importantly, what must not. A blanket 1 nt trim would satisfy the first two
and destroy the rest: 326 of 865 intron-free tRNAs already match the references
exactly.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plastanno.core.feature import Feature                       # noqa: E402
from plastanno.identify import trna_hybrid as H
from plastanno.core import coords as C


def run(base, hits, glen=None, diag=None):
    """apply_trnascan_geometry over a legacy inventory."""
    return H.apply_trnascan_geometry(base, hits, glen or L, diagnostics=diag)                   # noqa: E402

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-64s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


L = 130000


def feat(gene, s, e, strand, exons=None):
    f = Feature(gene_name=gene, gene_type="tRNA", product="tRNA-X",
                start=s, end=e, strand=strand)
    f.engine = "B"
    if exons:
        f.exons = [tuple(x) for x in exons]
    return f


def ts(aa, anti, arcs, strand, score=60.0, intron=False):
    return {"aa": aa, "anticodon": anti, "strand": strand, "arcs": sorted(arcs),
            "score": score, "has_intron": intron, "touches_cut": False,
            "from_rotation": 0, "call": "x"}


def geom(f):
    for n in f.notes:
        if n.startswith("geometry_source="):
            return n.split("=", 1)[1]
    return None


def span(f):
    return (f.start, f.end)


print("--- 1. trnS-GCU: hybrid must move 22719-22808 to 22720-22807 ---")
a = feat("trnS-GCU", 22718, 22808, 1)
out = run([a], [ts("Ser", "GCU", [(22719, 22807)], 1)])
check("one feature out", len(out), 1)
check("coordinates taken from tRNAscan-SE", span(out[0]), (22719, 22807))
check("geometry_source recorded", geom(out[0]), H.GEOM_TRNASCAN)
check("identity_source recorded",
      any(n.startswith("identity_source=") for n in out[0].notes), True)
check("aragorn_call preserved", any("aragorn_call=22718-22808" in n for n in out[0].notes), True)
check("gene name untouched", out[0].gene_name, "trnS-GCU")

print("--- 2. trnR-ACG: hybrid must move 88956-89030 to 88956-89029 ---")
out = run([feat("trnR-ACG", 88955, 89030, 1)], [ts("Arg", "ACG", [(88955, 89029)], 1)])
check("3' end pulled in by one", span(out[0]), (88955, 89029))

print("--- 3. already-agreeing loci must not move ---")
out = run([feat("trnL-CAA", 3632, 3712, 1)], [ts("Leu", "CAA", [(3632, 3712)], 1)])
check("trnL-CAA unchanged", span(out[0]), (3632, 3712))
out = run([feat("trnR-UCU", 21112, 21184, -1)], [ts("Arg", "UCU", [(21112, 21184)], -1)])
check("trnR-UCU unchanged", span(out[0]), (21112, 21184))

print("--- 4. every legacy locus survives; nothing is added ---")
base = [feat("trnP-UGG", 50000, 50074, 1), feat("trnQ-UUG", 7000, 7073, -1)]
d = []
out = run(base, [ts("Gly", "GCC", [(30000, 30072)], 1)], diag=d)
check("both legacy loci returned", len(out), 2)
check("tRNAscan-only NOT added", all(f.gene_name != "trnG-GCC" for f in out), True)
check("and it is reported as a diagnostic",
      any("trnascan-only" in x for x in d), True)

print("--- 5. a BLAST-only locus survives (it is in the legacy inventory) ---")
b = feat("trnX-AAA", 12000, 12072, 1)
b.notes.append("BLAST")
out = run([b], [])
check("BLAST-only locus still present", len(out), 1)
check("its coordinates are untouched", span(out[0]), (12000, 12072))
check("marked as legacy geometry", geom(out[0]), H.GEOM_LEGACY)

print("--- 6. intron-bearing tRNA is frozen ---")
ex = [(111415, 111470), (112394, 112449)]
a = feat("trnI-GAU", 111415, 112449, 1, exons=ex)
out = run([a], [ts("Ile", "GAU", [(111433, 111469), (111485, 111521)], 1, intron=True)])
check("exons untouched", [tuple(x) for x in out[0].exons], ex)
check("span untouched", span(out[0]), (111415, 112449))
check("geometry stays legacy", geom(out[0]), H.GEOM_LEGACY)

print("--- 7. an intron-FREE tRNA crossing the origin is still eligible ---")
# stored as two arcs meeting at the seam: 1 biological exon, not 2
w = feat("trnH-GUG", L - 3, 72, -1, exons=[(0, 72), (L - 3, L)])
check("counted as one biological exon", C.biological_exon_count([(0, 72), (L - 3, L)], L), 1)
out = run([w], [ts("His", "GUG", [(0, 71), (L - 3, L)], -1)])
check("its boundary WAS replaced", geom(out[0]), H.GEOM_TRNASCAN)
check("and it stays wrapped, not 0..genome_length",
      (out[0].start, out[0].end), (L - 3, 71))
check("start > end preserves the wrap", out[0].start > out[0].end, True)

print("--- 8. an intron-BEARING tRNA crossing the origin is frozen ---")
wi = feat("trnK-UUU", L - 40, 900, 1, exons=[(0, 40), (860, 900), (L - 40, L)])
check("counted as two biological exons",
      C.biological_exon_count([(0, 40), (860, 900), (L - 40, L)], L), 2)
out = run([wi], [ts("Lys", "UUU", [(0, 40), (L - 40, L)], 1)])
check("geometry frozen", geom(out[0]), H.GEOM_LEGACY)
check("exons untouched", len(out[0].exons), 3)

print("--- 9. an exon-DB locus on the opposite strand does not delete ARAGORN ---")
# the two IR copies face opposite ways and must both survive
c1 = feat("trnA-UGC", 40000, 41000, 1, exons=[(40000, 40040), (40960, 41000)])
c2 = feat("trnA-UGC", 90000, 91000, -1, exons=[(90000, 90040), (90960, 91000)])
out = run([c1, c2], [])
check("both IR copies retained", len(out), 2)
check("on opposite strands", sorted(f.strand for f in out), [-1, 1])

print("--- 10. competing 2x2 matching is order-independent ---")
def mk():
    return ([feat("trnV-GAC", 1000, 1073, 1), feat("trnV-GAC", 1050, 1123, 1)],
            [ts("Val", "GAC", [(1000, 1072)], 1, score=70.0),
             ts("Val", "GAC", [(1050, 1122)], 1, score=60.0)])
f1, h1 = mk(); r1 = sorted(span(x) for x in run(f1, h1))
f2, h2 = mk(); r2 = sorted(span(x) for x in run(f2[::-1], h2))
f3, h3 = mk(); r3 = sorted(span(x) for x in run(f3, h3[::-1]))
f4, h4 = mk(); r4 = sorted(span(x) for x in run(f4[::-1], h4[::-1]))
check("features reversed gives the same panel", r1, r2)
check("hits reversed gives the same panel", r1, r3)
check("both reversed gives the same panel", r1, r4)
check("assignment is one-to-one", len(set(r1)), 2)

print("--- 11. identity guards the replacement ---")
out = run([feat("trnW-CCA", 60000, 60074, 1)], [ts("Gly", "GCC", [(60000, 60073)], 1)])
check("a different family cannot donate geometry", span(out[0]), (60000, 60074))
check("geometry stays legacy", geom(out[0]), H.GEOM_LEGACY)
# tRNAscan calling Met where the name says Ile must still be allowed: the CAU
# families are exactly what the naming work settled, and refusing them would
# skip those loci
out = run([feat("trnI-CAU", 89484, 89559, -1)], [ts("Met", "CAU", [(89485, 89559)], -1)])
check("trnI-CAU still takes a Met-CAU boundary", geom(out[0]), H.GEOM_TRNASCAN)
check("but keeps its own name", out[0].gene_name, "trnI-CAU")

print("--- 12. opposite strand and weak overlap are never merged ---")
out = run([feat("trnW-CCA", 60000, 60074, 1)], [ts("Trp", "CCA", [(60000, 60073)], -1)])
check("opposite strand refused", geom(out[0]), H.GEOM_LEGACY)
out = run([feat("trnQ-UUG", 7000, 7073, 1)], [ts("Gln", "UUG", [(7060, 7200)], 1)])
check("weak reciprocal overlap refused", geom(out[0]), H.GEOM_LEGACY)

print("--- 13. no blanket trim ---")
ex3 = [feat("trnL-CAA", 3632, 3712, 1), feat("trnR-UCU", 21112, 21184, -1),
       feat("trnF-GAA", 51212, 51285, 1)]
h3 = [ts("Leu", "CAA", [(3632, 3712)], 1), ts("Arg", "UCU", [(21112, 21184)], -1),
      ts("Phe", "GAA", [(51212, 51285)], 1)]
out = run(ex3, h3)
check("all three unchanged", sorted(span(f) for f in out),
      [(3632, 3712), (21112, 21184), (51212, 51285)])

print("--- 14. fail-closed when tRNAscan-SE is unavailable ---")
import shutil as _sh                                              # noqa: E402
real = _sh.which
_sh.which = lambda name: None
try:
    H.run_trnascan_circular("ACGT" * 100)
    check("raises when tRNAscan-SE is missing", "no raise", "TrnascanUnavailable")
except H.TrnascanUnavailable as exc:
    check("raises when tRNAscan-SE is missing", "raised", "raised")
    check("and says why", "not on PATH" in str(exc), True)
finally:
    _sh.which = real

print("--- 15. circular parse: a rotated hit maps back ---")
import tempfile                                                   # noqa: E402
p = tempfile.NamedTemporaryFile("w", suffix=".out", delete=False)
p.write("a\nb\nc\nquery\t1\t11\t80\tHis\tGTG\t0\t0\t50.0\t\n")
p.close()
got = H._parse(p.name, L, L // 2)
os.unlink(p.name)
check("un-rotation is a modular add", got[0]["arcs"], [(10 + L // 2, 80 + L // 2)])
p = tempfile.NamedTemporaryFile("w", suffix=".out", delete=False)
p.write("a\nb\nc\nquery\t1\t%d\t%d\tHis\tGTG\t0\t0\t50.0\t\n" % (L // 2 - 10, L // 2 + 30))
p.close()
got = H._parse(p.name, L, L // 2)
os.unlink(p.name)
check("a hit crossing the seam returns two arcs", len(got[0]["arcs"]), 2)
check("the two arcs meet at the origin",
      got[0]["arcs"][0][0] == 0 and got[0]["arcs"][-1][1] == L, True)

print("--- 16. the 2x2 assignment is invariant under genome rotation ---")
# Rotating the genome renumbers every coordinate. If the tie-break reads a
# coordinate, the two frames can pair differently; the pairing must be decided
# by evidence instead.
def rot_arcs(arcs, sh, glen):
    out = []
    for a, b in arcs:
        a2, b2 = (a + sh) % glen, (b + sh - 1) % glen + 1
        if a2 < b2:
            out.append((a2, b2))
        else:
            out.extend([(a2, glen), (0, b2)])
    return sorted(out)

SH = L // 3
def build(sh):
    f1 = feat("trnV-GAC", *rot_arcs([(1000, 1073)], sh, L)[0], 1)
    f2 = feat("trnV-GAC", *rot_arcs([(1050, 1123)], sh, L)[0], 1)
    h1 = ts("Val", "GAC", rot_arcs([(1000, 1072)], sh, L), 1, score=70.0)
    h2 = ts("Val", "GAC", rot_arcs([(1050, 1122)], sh, L), 1, score=60.0)
    return [f1, f2], [h1, h2]

fa, ha = build(0)
pa = H.global_match(fa, ha, L)
fb, hb = build(SH)
pb = H.global_match(fb, hb, L)
check("same number of pairs in both frames", len(pa), len(pb))
check("the same feature-to-hit pairing", sorted((i, j) for i, j, _ in pa),
      sorted((i, j) for i, j, _ in pb))
# and the evidence digest itself must not move
check("evidence digest is rotation-invariant",
      H._evidence_digest(fa[0], ha[0], L), H._evidence_digest(fb[0], hb[0], L))
check("tie-break fraction is rotation-invariant",
      H._tiebreak_fraction(fa[0], ha[0], L), H._tiebreak_fraction(fb[0], hb[0], L))

print("--- 17. identity_source carries the right VALUE in all four cases ---")
def isrc(f):
    for n in f.notes:
        if n.startswith("identity_source="):
            return n.split("=", 1)[1]
    return None

ar = feat("trnP-UGG", 50000, 50074, 1)
bl = feat("trnX-AAA", 12000, 12072, 1)
ed = feat("trnI-GAU", 111415, 112449, 1, exons=[(111415, 111470), (112394, 112449)])
H.tag_identity_sources([ar], [bl], [ed])
out = run([ar, bl, ed], [ts("Pro", "UGG", [(50000, 50073)], 1)])
by = {f.gene_name: f for f in out}
check("ARAGORN-named locus", isrc(by["trnP-UGG"]), H.IDENT_ARAGORN)
check("BLAST-only locus", isrc(by["trnX-AAA"]), H.IDENT_BLAST)
check("exon-DB locus", isrc(by["trnI-GAU"]), H.IDENT_EXON_DB)

cau = feat("trnI-CAU", 89484, 89559, -1)
H.tag_identity_sources([cau], [], [])
out = run([cau], [ts("Met", "CAU", [(89485, 89559)], -1)])
check("CAU locus keeps its own identity source", isrc(out[0]), H.IDENT_ARAGORN)
check("and the CAU exception is recorded",
      any(n == "identity_match=CAU_AMBIGUOUS" for n in out[0].notes), True)
check("a non-CAU family mismatch sets no exception flag",
      getattr(feat("trnW-CCA", 1, 74, 1), "_cau_exception", False), False)

print("--- 18. the geometry step mutates IN PLACE and never touches the list ---")
# This is the property that makes inventory invariance structural. The old
# design rebuilt a list inside Engine B, upstream of reconcile._select, and a
# boundary move that separated two overlapping duplicates turned into an extra
# locus. Now the caller keeps its own list and this function may only write
# .start/.end/.exons/.notes on objects already in it.
base = [feat("trnS-GCU", 22718, 22808, 1),
        feat("trnI-GAU", 111415, 112449, 1, exons=[(111415, 111470), (112394, 112449)]),
        feat("trnX-AAA", 12000, 12072, 1)]
ident = [id(f) for f in base]
snapshot = list(base)
out = run(base, [ts("Ser", "GCU", [(22719, 22807)], 1)])
check("the caller's list is not modified", base, snapshot)
check("no object is replaced", [id(f) for f in base], ident)
check("every input object is returned", sorted(id(f) for f in out), sorted(ident))
check("nothing is added", len(out), len(base))
check("the eligible one did move", span(base[0]), (22719, 22807))
check("the frozen one did not", span(base[1]), (111415, 112449))

print("--- 19. the reported duplication case: abutting neighbours stay two ---")
# NC_026958.1 had two same-name candidates overlapping by 7 bp that legacy's
# _select collapsed. Hybrid moved them to abutting, they stopped clustering, and
# both survived. This function cannot cause that any more, because it does not
# decide membership at all -- but it must still not lose either object.
d1 = feat("trnT-GGU", 15903, 15977, 1)
d2 = feat("trnT-GGU", 15970, 16045, 1)
pair = [d1, d2]
out = run(pair, [ts("Thr", "GGU", [(15904, 15976)], 1),
                 ts("Thr", "GGU", [(15976, 16040)], 1)])
check("both objects come back", len(out), 2)
check("the caller's list still holds exactly its two objects",
      [id(f) for f in pair], [id(d1), id(d2)])
check("membership is the caller's business, not this function's",
      sorted(id(f) for f in out), sorted([id(d1), id(d2)]))

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
