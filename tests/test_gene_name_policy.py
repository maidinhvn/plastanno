#!/usr/bin/env python3
"""Output gene names follow the majority GenBank usage: psbN, ycf3, ycf4, clpP.

Run DIRECTLY and check $?. Exit 0 = all pass.

Decision of 2026-09-29, on parked/arbiter_oracle/gene_name_timeline.txt: the newer
names pbf1, pafI, pafII and clpP1 appear from 2014-15 and reach 23-28% of records
submitted in 2024; the majority still write psbN, ycf3, ycf4 and clpP. The
synonym table used to map psbN -> pbf1, the minority name.

The rename happens in step 6 only. Doing it before reconciliation was tried and
failed its criterion: it changed which call wins, not only its name
(parked/gene_name_policy/RESULT.md). Because the databases name some genes twice
(pbf1/psbN, pafI/ycf3, pafII/ycf4), a renamed call can land on a locus that
already has a call of the target name; that renamed call is dropped.
"""
import inspect
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.core.feature import Feature                     # noqa: E402
from plastanno.core import reconcile as RC                     # noqa: E402
from plastanno.annotate import special_cases as SC             # noqa: E402
from plastanno import paths                                    # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-66s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


def F(name, s, e, engine, **kw):
    f = Feature(gene_name=name, gene_type="CDS", start=s, end=e, strand=1,
                engine=engine, **kw)
    f.exons = [(s, e)]
    return f


print("--- the table ---")
for old, new in (("pbf1", "psbN"), ("pafI", "ycf3"), ("pafII", "ycf4"),
                 ("clpP1", "clpP"), ("clpP2", "clpP")):
    check("%s -> %s" % (old, new), SC.SYNONYMS.get(old), new)
check("no majority name is renamed away",
      [g for g in ("psbN", "ycf3", "ycf4", "clpP") if g in SC.SYNONYMS], [])
check("no target is itself renamed (one pass is final)",
      [v for v in SC.SYNONYMS.values() if v in SC.SYNONYMS], [])

print("--- apply_synonyms ---")
fs = [F("pbf1", 100, 232, "A"), F("psbA", 500, 800, "B"), F("pafII", 900, 1455, "B")]
n = SC.apply_synonyms(fs)
check("renames exactly the synonyms", ([f.gene_name for f in fs], n),
      (["psbN", "psbA", "ycf4"], 2))
check("keeps the engine's name in the notes",
      any("pbf1" in x for x in fs[0].notes), True)
check("leaves other features' notes alone", fs[1].notes, [])

print("--- the catalog the tool reads ---")
cat = json.load(open(os.path.join(str(paths.config_dir()), "gene_catalog.json")))
check("psbN carries pbf1's expected_len", cat["psbN"].get("expected_len"),
      cat["pbf1"].get("expected_len"))
check("psbN carries pbf1's exon count", cat["psbN"].get("n_exons"),
      cat["pbf1"].get("n_exons"))
check("psbN's product is the protein name", cat["psbN"].get("product"),
      "photosystem II protein N")
for g in ("psbN", "ycf3", "ycf4", "clpP"):
    check("%s has an expected_len for the length filter" % g,
          bool(cat[g].get("expected_len")), True)

print("--- step 6 renames, and drops a renamed duplicate ---")
def N(feats):
    out, ch = SC.normalize_names(feats, {"LSC": (0, 10000)})
    return sorted((f.gene_name, f.start, f.end) for f in out), ch
got, ch = N([F("pbf1", 100, 232, "A")])
check("a lone pbf1 becomes psbN and is kept", got, [("psbN", 100, 232)])
got, ch = N([F("ycf3", 1000, 1507, "A"), F("pafI", 1100, 1300, "B")])
check("ycf3 + overlapping pafI -> the original ycf3 only", got, [("ycf3", 1000, 1507)])
check("  ... and the drop is recorded", any("Dropped" in c for c in ch), True)
got, ch = N([F("psbN", 100, 232, "B"), F("pbf1", 106, 232, "A")])
check("psbN + overlapping pbf1 -> the original psbN only", got, [("psbN", 100, 232)])
got, ch = N([F("ycf3", 1000, 1507, "A"), F("pafI", 6000, 6507, "B")])
check("non-overlapping copies are both kept (an IR pair is not a duplicate)",
      got, [("ycf3", 1000, 1507), ("ycf3", 6000, 6507)])

print("--- attempt 1 must not come back ---")
# Renaming the engines' calls before reconciliation moved coordinates in 8 of
# 30 dev genomes and lost an exact psbN (parked/gene_name_policy/RESULT.md).
import plastanno.pipeline as PL                                # noqa: E402
src = inspect.getsource(PL.run)
check("pipeline.run does not rename before step 5", "apply_synonyms" in src, False)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
