"""The single definition of which genes the inverted repeat duplicates.

Before this module the set existed twice, with different contents:

    identify/engine_a.py  IR_GENES   16 genes, used to decide which genes to
                                     SEARCH in both IR copies
    core/reconcile.py     ir_genes    8 genes, a local, used to decide which
                                     genes may KEEP both copies out of the
                                     conflict bin

The second is a strict subset of the first: it omits rrn16, rrn23, rrn4.5, rrn5,
trnR-ACG, trnN-GUU, trnL-CAA and trnV-GAC. The consequence is a latent defect in
the legacy path -- when the two engines disagree about the placement of one of
those eight, the conflict handler keeps a single copy and the IR pair is halved.

That defect is NOT fixed in the legacy path here. Legacy output is frozen: the
H4-30 audit and every regression baseline were produced with the eight-gene set,
and widening it would silently change annotations the whole benchmark record
depends on. So the historical set is preserved verbatim as
LEGACY_IR_CONFLICT_GENES, named and documented rather than hidden in a function
body, and the pooled path uses the complete set.

CLAUDE.md already warned that the two sets "must stay consistent". They were not.
"""

#: Genes carried in both copies of the inverted repeat. The canonical set; used
#: for IR-aware search and by the pooled reconciliation path.
IR_DUPLICATED_GENES = frozenset({
    "rpl2", "rpl23", "ndhB", "rps7", "ycf2",
    "trnI-CAU", "trnI-GAU", "trnA-UGC",
    "trnR-ACG", "trnN-GUU", "trnL-CAA", "trnV-GAC",
    "rrn16", "rrn23", "rrn4.5", "rrn5",
})

#: FROZEN. The eight genes the legacy conflict handler has always used. Do not
#: extend this: legacy output must stay reproducible against the frozen H4-30
#: run. New work belongs in IR_DUPLICATED_GENES.
LEGACY_IR_CONFLICT_GENES = frozenset({
    "rpl2", "rpl23", "ndhB", "rps7", "ycf2",
    "trnI-CAU", "trnI-GAU", "trnA-UGC",
})

#: The database build script assigns catalog regions from its own copy, which
#: additionally contains orf70. Content is preserved exactly: changing it would
#: change gene_catalog.json on the next build, and the database is frozen.
BUILD_CATALOG_IR_GENES = IR_DUPLICATED_GENES | {"orf70"}

#: What the legacy set is missing, kept as data so a test can assert the gap is
#: exactly this and has not silently grown or shrunk.
LEGACY_IR_CONFLICT_GAP = IR_DUPLICATED_GENES - LEGACY_IR_CONFLICT_GENES
