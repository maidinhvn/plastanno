#!/usr/bin/env python3
"""The GFF3 gene line for a feature that crosses the origin.

Found by the cross-format verifier on NC_045219.1 during the hybrid holdout, and
reproduced deterministically:

    trnH-GUG        .gb 1 copy [((0,62),(154576,154579))], .gff3 0
    trnH-GUG_154576 .gb 0 copies,                          .gff3 1

Identical coordinates, different names. `write_gff3` built the gene line's spans
from `ann.exons`, so a wrapped feature whose exons list is empty — every
intron-free tRNA that `--trna-mode hybrid` re-seats across the origin — got NO
gene line. With no gene line there is no `ID -> Name` entry, so `verify.py` could
not resolve `Parent=` and fell back to reading the feature's own ID as its gene
name.

The same expression was wrong a second way: a wrapped multi-exon gene got one
gene line per exon, with the introns cut out. A GFF3 gene line is the outer span,
introns included.

These tests pin both, and pin that unwrapped features did not move.
"""
import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plastanno.core.feature import Feature          # noqa: E402
from plastanno.output.writers import write_gff3     # noqa: E402

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-68s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


GLEN = 154_600


def gff(annotations, glen=GLEN):
    """Write the annotations and return the parsed non-comment rows."""
    fd, path = tempfile.mkstemp(suffix=".gff3")
    os.close(fd)
    try:
        write_gff3(annotations, glen, "TEST.1", path)
        rows = []
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) == 9:
                    rows.append(c)
        return rows
    finally:
        os.unlink(path)


def attrs(col9):
    return dict(kv.split("=", 1) for kv in col9.split(";") if "=" in kv)


def at0(rows, key):
    """attrs of the first row, or a marker -- a missing row is a FAIL to report,
    not an IndexError that hides every later section."""
    return attrs(rows[0][8]).get(key, "<no row>") if rows else "<no row>"


def gene_rows(rows):
    return [r for r in rows if r[2] == "gene"]


def spans(rows):
    """1-based inclusive (start, end) pairs, as written."""
    return sorted((int(r[3]), int(r[4])) for r in rows)


# The exact feature that failed: intron-free, wrapped, exons never populated.
WRAPPED_BARE = Feature(
    gene_name="trnH-GUG", gene_type="tRNA", product="tRNA-His",
    start=154_576, end=62, strand=-1, exons=[],
    engine="B", confidence=0.90, flag="HIGH",
)

# A wrapped gene that DOES carry exons, with an intron between them.
WRAPPED_SPLICED = Feature(
    gene_name="trnK-UUU", gene_type="tRNA", product="tRNA-Lys",
    start=154_000, end=500, strand=1,
    exons=[(154_000, 154_037), (100, 500)],
    engine="B", confidence=0.85, flag="HIGH", has_intron=True,
)

ORDINARY = Feature(
    gene_name="psbA", gene_type="CDS", product="photosystem II protein D1",
    start=1_000, end=2_062, strand=1, exons=[],
    engine="A", confidence=0.95, flag="HIGH",
)


print("--- 1. a wrapped feature with NO exons still gets its gene line ---")
rows = gff([WRAPPED_BARE])
g = gene_rows(rows)
check("a gene line exists at all (this is the bug)", len(g) > 0, True)
check("it is split into exactly two arcs", len(g), 2)
check("the arcs are the outer span cut at the origin",
      spans(g), [(1, 62), (154_577, 154_600)])
check("both arcs share one ID", len({attrs(r[8])["ID"] for r in g}), 1)
check("the ID carries the disambiguating suffix",
      at0(g, "ID"), "trnH-GUG_154576")
check("the Name is the bare gene name",
      at0(g, "Name"), "trnH-GUG")

print("--- 2. Parent= resolves to a Name=, which is what verify.py walks ---")
# Reproduce verify.py's two passes rather than trusting that they agree.
name_of = {}
for r in rows:
    if r[2] == "gene":
        a = attrs(r[8])
        name_of[a["ID"]] = a.get("Name", a["ID"])
resolved = set()
for r in rows:
    if r[2] in ("CDS", "tRNA", "rRNA"):
        a = attrs(r[8])
        ident = a.get("Parent") or a.get("ID", "?")
        resolved.add(name_of.get(ident, a.get("Name", ident)))
check("the tRNA lines resolve to the gene name, not to the ID",
      resolved, {"trnH-GUG"})
check("  ... and never to the suffixed form",
      any("_154576" in n for n in resolved), False)

print("--- 3. a wrapped SPLICED gene: two gene arcs, introns included ---")
rows = gff([WRAPPED_SPLICED])
g = gene_rows(rows)
check("two gene lines, not one per exon", len(g), 2)
check("the arcs span the whole gene, intron included",
      spans(g), [(1, 500), (154_001, 154_600)])
# The intron must still be absent from the tRNA lines: only the GENE line fills it.
t = [r for r in rows if r[2] == "tRNA"]
check("the tRNA lines keep the exon structure", spans(t),
      [(101, 500), (154_001, 154_037)])

print("--- 4. unwrapped features are untouched ---")
rows = gff([ORDINARY])
g = gene_rows(rows)
check("one gene line", len(g), 1)
check("spanning start..end", spans(g), [(1_001, 2_062)])
check("ID keeps the familiar <gene>_<start> form",
      at0(g, "ID"), "psbA_1000")

print("--- 5. a feature ending exactly at the origin emits no empty arc ---")
seam = Feature(gene_name="trnX-XXX", gene_type="tRNA", start=154_500, end=0,
               strand=1, exons=[], engine="B", confidence=0.8, flag="HIGH")
g = gene_rows(gff([seam]))
check("one arc only, no zero-length second arc", len(g), 1)
check("it runs to the end of the genome", spans(g), [(154_501, 154_600)])

print("--- 6. repeated names still get unique IDs (the rps12 regression) ---")
a = Feature(gene_name="rps12", gene_type="CDS", start=100, end=200, strand=1,
            exons=[], engine="A", confidence=0.9, flag="HIGH")
b = Feature(gene_name="rps12", gene_type="CDS", start=100, end=200, strand=1,
            exons=[], engine="A", confidence=0.9, flag="HIGH")
ids = [attrs(r[8])["ID"] for r in gene_rows(gff([a, b]))]
check("two distinct IDs", len(set(ids)), 2)
check("the second is suffixed", sorted(ids), ["rps12_100", "rps12_100.2"])


print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
