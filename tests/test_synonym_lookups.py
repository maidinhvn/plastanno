#!/usr/bin/env python3
"""Every gene renamed in step 6 still finds its data under the new name.

Run DIRECTLY and check $?. Exit 0 = all pass.

Step 6 renames genes to their majority names (SYNONYMS) before anything else. Every later step
looks its data up by the NEW name: reference proteins (ORF completion, the start-codon pass),
the gene catalog (expected length, exon count, protein name), the special start codons, the
exon panel and its templates. psbN's reference proteins lived only under pbf1, and for two
releases every psbN silently skipped ORF completion. This test walks every SYNONYMS pair and
checks that each lookup made after the rename finds at least what the old name would have
found, so a new synonym, or a rebuilt database, cannot reopen that gap unnoticed.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.annotate.special_cases import SYNONYMS, _ref_protein_file   # noqa: E402
from plastanno.core.reconcile import SPECIAL_START_CODONS                 # noqa: E402
from plastanno.paths import config_dir, db_root                            # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-72s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


cat = json.load(open(config_dir() / "gene_catalog.json"))
meta = json.load(open(config_dir() / "boundary_db" / "meta.json"))
tpl = json.load(open(config_dir() / "exon_templates.json"))
pdb = db_root() / "protein_db"
have_db = pdb.is_dir() and any(pdb.glob("*.fasta"))

for old, new in sorted(SYNONYMS.items()):
    print("--- %s -> %s ---" % (old, new))
    co, cn = cat.get(old) or {}, cat.get(new) or {}
    check("the catalog has an entry for %s" % new, bool(cn), True)
    lost = sorted(k for k, v in co.items() if k != "synonym_of" and v is not None and cn.get(k) is None)
    check("... with every field the %s entry has" % old, lost, [])
    check("special start codons: nothing under %s that %s lacks" % (old, new),
          sorted(set(SPECIAL_START_CODONS.get(old, ())) - set(SPECIAL_START_CODONS.get(new, ()))), [])
    check("exon panel: %s is not keyed only under %s" % (new, old), old in meta and new not in meta, False)
    check("exon templates: %s is not keyed only under %s" % (new, old), old in tpl and new not in tpl, False)
    if have_db:
        has_old = (pdb / ("%s.fasta" % old)).exists()
        found = _ref_protein_file(str(pdb), new)
        check("reference proteins: a %s lookup finds a file whenever %s has one" % (new, old),
              found is not None or not has_old, True)

if not have_db:
    print("  (no protein database at %s: the protein lookups were not checked)" % pdb)
print("\n%d checks, %d failed" % (RUN[0], len(FAIL)))
sys.exit(1 if FAIL else 0)
