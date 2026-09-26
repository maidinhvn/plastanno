#!/usr/bin/env python3
"""Engine A tries the query's closest relatives first, and changes nothing else.

Step 3 BLASTs the genome against genus_reps and hands the ranked list to
run_engine_a, which used to discard it: `relatives` appeared in a comment and in
a signature and never in a body. Every genome was annotated against whichever
five references came first in the FASTA.

The fix reorders. These tests pin that it is ONLY a reorder — the same proteins,
no additions, no removals — and that a genome with nothing to reorder keeps the
order it had. That second property is what caught nothing in the 30-genome run
(19 of 19 such genomes produced byte-identical CDS output) and is the cheapest
way to notice if the ordering ever starts doing something else.

The reordering is exercised through the real code path by monkeypatching the
Exonerate call to record which proteins it was handed, rather than by reproducing
the sort here — a test that reimplements the logic it checks proves nothing.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Bio.Seq import Seq                              # noqa: E402
from Bio.SeqRecord import SeqRecord                  # noqa: E402
import plastanno.identify.engine_a as EA             # noqa: E402

RUN = [0]
FAIL = []


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-66s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


GENE = "atpA"
# Five DB proteins in file order. NC_000005 is the one a relative search finds.
DB = [SeqRecord(Seq("M" * 100), id="NC_000001.1_%s" % GENE, description=""),
      SeqRecord(Seq("M" * 100), id="NC_000002.1_%s" % GENE, description=""),
      SeqRecord(Seq("M" * 100), id="NC_000003.1_%s" % GENE, description=""),
      SeqRecord(Seq("M" * 100), id="NC_000004.1_%s" % GENE, description=""),
      SeqRecord(Seq("M" * 100), id="NC_000005.1_%s" % GENE, description="")]

CATALOG = {GENE: {"type": "CDS", "region": "LSC", "n_exons": 1,
                  "synonym_of": None, "product": GENE, "expected_len": 300}}


def run_with(relatives, db=DB, cap_seen=None):
    """Call the real run_exonerate_gene; record the order Exonerate was handed."""
    seen = []

    def fake_region(genome_seq, protein_seq, gene_name, prot_acc,
                    region_start, region_end, genome_len=None):
        seen.append(prot_acc)
        return []

    def fake_parse(path, fmt):
        return list(db)

    real_region, real_parse = EA.run_exonerate_region, None
    EA.run_exonerate_region = fake_region
    import Bio.SeqIO as SIO
    real_parse = SIO.parse
    SIO.parse = fake_parse
    # the function resolves the fasta path and checks .exists(); point it at a
    # directory that has the file so the branch is taken
    try:
        EA.run_exonerate_gene(
            genome_seq="A" * 5000, gene_name=GENE,
            protein_db="database/protein_db", ir_boundaries={"LSC": (0, 5000)},
            gene_catalog=CATALOG, threads=1, relatives=relatives)
    finally:
        EA.run_exonerate_region = real_region
        SIO.parse = real_parse
    return seen


print("--- 1. with no relatives, file order is preserved ---")
seen = run_with(None)
check("the first five in file order are tried",
      seen, ["NC_00000%d.1_%s" % (i, GENE) for i in range(1, 6)])

print("--- 2. a relative is promoted to the front ---")
seen = run_with([{"accession": "NC_000005.1"}])
check("the relative is tried first", seen[0], "NC_000005.1_%s" % GENE)
check("and the rest keep file order",
      seen[1:], ["NC_00000%d.1_%s" % (i, GENE) for i in range(1, 5)])

print("--- 3. several relatives are promoted in rank order ---")
seen = run_with([{"accession": "NC_000004.1"}, {"accession": "NC_000002.1"}])
check("rank 1 first, rank 2 second",
      seen[:2], ["NC_000004.1_%s" % GENE, "NC_000002.1_%s" % GENE])
check("the remainder keeps file order",
      seen[2:], ["NC_000001.1_%s" % GENE, "NC_000003.1_%s" % GENE,
                 "NC_000005.1_%s" % GENE])

print("--- 4. a relative that is not in the database changes nothing ---")
seen = run_with([{"accession": "NC_999999.1"}])
check("order is exactly the unrelativised order",
      seen, ["NC_00000%d.1_%s" % (i, GENE) for i in range(1, 6)])

print("--- 5. it is a reorder: nothing is added, dropped or duplicated ---")
for rel in (None, [{"accession": "NC_000005.1"}],
            [{"accession": "NC_000003.1"}, {"accession": "NC_999999.1"}]):
    seen = run_with(rel)
    check("five proteins tried, all distinct, all from the DB (%s)"
          % ("none" if rel is None else rel[0]["accession"]),
          (len(seen), len(set(seen)), set(seen) <= {r.id for r in DB}),
          (5, 5, True))

print("--- 6. the cap is still five even with more relatives than that ---")
big = [SeqRecord(Seq("M" * 100), id="NC_0000%02d.1_%s" % (i, GENE), description="")
       for i in range(1, 11)]
seen = run_with([{"accession": "NC_0000%02d.1" % i} for i in (9, 8, 7)], db=big)
check("still five", len(seen), 5)
check("and they start with the three relatives in rank order",
      seen[:3], ["NC_000009.1_%s" % GENE, "NC_000008.1_%s" % GENE,
                 "NC_000007.1_%s" % GENE])


print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
