#!/usr/bin/env python3
"""A gene renamed in step 6 still finds its reference proteins.

Run DIRECTLY and check $?. Exit 0 = all pass.

Step 6 renames pbf1 to psbN before anything else, but the protein database keeps psbN's
reference proteins in pbf1.fasta. Every later lookup was by the new name, found nothing, and
was skipped silently: from 3.0.1 on no psbN was ORF-completed, and a psbN call four codons short
of its stop stayed 12 bp short (the regression the 3.0.1 benchmark saw). The lookup now falls back
to the file of a name SYNONYMS maps to the gene; the new name wins when both exist.
"""
import copy
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.core.feature import Feature                                   # noqa: E402
from plastanno.core import coords as C                                       # noqa: E402
from plastanno.annotate import special_cases as SC                           # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-74s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


PSBN = "METATLVAIFISGLLVSFTGYALYTAFGQPSQQLRDPFEEHGD"      # 43 aa
CODON = {"M": "ATG", "E": "GAA", "T": "ACT", "A": "GCT", "L": "CTT", "V": "GTT", "I": "ATT",
         "F": "TTT", "S": "TCT", "G": "GGT", "Y": "TAT", "Q": "CAA", "P": "CCT", "R": "CGT",
         "D": "GAT", "H": "CAT"}


def db_with(*names, prot=PSBN):
    d = tempfile.mkdtemp(prefix="renamed_db_")
    for n in names:
        with open(os.path.join(d, n + ".fasta"), "w") as fh:
            for i in range(4):
                fh.write(">%s_%d\n%s\n" % (n, i, prot))
    return d


name = lambda p: os.path.basename(str(p)) if p else None

print("--- which file is read ---")
d = db_with("pbf1")
check("psbN with only pbf1.fasta: pbf1.fasta", name(SC._ref_protein_file(d, "psbN")), "pbf1.fasta")
check("... and its profile loads", SC._load_ref_profile(d, "psbN")[0][1], 4)
check("pbf1 itself: pbf1.fasta", name(SC._ref_protein_file(d, "pbf1")), "pbf1.fasta")
d2 = db_with("pbf1", "psbN")
check("both files: the new name wins", name(SC._ref_protein_file(d2, "psbN")), "psbN.fasta")
d3 = db_with("rbcL")
check("neither file: None", SC._ref_protein_file(d3, "psbN"), None)
check("a gene with no synonym and no file: None", SC._ref_protein_file(d3, "atpA"), None)
check("a gene with its own file: that file", name(SC._ref_protein_file(d3, "rbcL")), "rbcL.fasta")
d4 = db_with("clpP1")
check("clpP falls back to clpP1.fasta (a synonym of clpP)", name(SC._ref_protein_file(d4, "clpP")), "clpP1.fasta")

print("--- the regression: a psbN call four codons short of its stop ---")
L = 3000
G = list("GCC" * (L // 3))                          # no stop and no ATG in any frame
orf = "".join(CODON[a] for a in PSBN) + "TAA"
G[900:900 + len(orf)] = list(orf)
SEQ = "".join(G)
call = Feature(gene_name="pbf1", gene_type="CDS", start=900, end=900 + len(orf) - 12, strand=1,
               engine="A", flag="HIGH")
call.exons = [(call.start, call.end)]
anns, _ = SC.normalize_names([copy.deepcopy(call)], None, genome_seq=SEQ)
check("step 6 renames the call to psbN first", anns[0].gene_name, "psbN")
anns, ch = SC.complete_single_exon_orfs(anns, SEQ, {}, protein_db=d)
f = anns[0]
check("ORF completion now reaches the stop", (f.start, f.end, C.extract(SEQ, f, L)[-3:]),
      (900, 900 + len(orf), "TAA"))
check("... and says so", any("ORF boundary completed" in str(n) for n in f.notes), True)
anns, _ = SC.normalize_names([copy.deepcopy(call)], None, genome_seq=SEQ)
anns, ch = SC.complete_single_exon_orfs(anns, SEQ, {}, protein_db=d3)
check("without any psbN proteins it is left as it was", (anns[0].end, ch), (900 + len(orf) - 12, []))

print("\n%d checks, %d failed" % (RUN[0], len(FAIL)))
sys.exit(1 if FAIL else 0)
