#!/usr/bin/env python3
"""For a panel gene marked keep_ends (rpl2), the exon panel moves the junctions, never the ends.

Run DIRECTLY and check $?. Exit 0 = all pass.

refine_splice transfers exon boundaries from a panel of reference exons. For rpl2 the panel fixes
the junction (Exonerate starts exon 2 one codon late against almost every reference), but its
references start on ATG where many lineages start on an edited ACG. Transferring the start made
those rpl2 look complete when they were not, which is why an earlier rpl2 panel was not accepted.
With "keep_ends" in the panel's meta.json, the CDS keeps the 5' and 3' ends it came in with and
only its junctions follow the panel; the start-codon and terminal-stop passes see to the ends.

The panel search itself (BLAST) is replaced by a fixed hit list, so each case knows its answer.
"""
import copy
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.core.feature import Feature                                   # noqa: E402
from plastanno.core import coords as C                                       # noqa: E402
from plastanno.annotate import refine_splice as RS                           # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-76s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
L = 6000
G = list("GCC" * (L // 3))                        # no stop codon in any frame on either strand
E1 = "ATG" + "GCC" * 3                             # exon 1: [1000, 1012)
INTRON = "GTGCG" + "GCC" * 60 + "TTAC"             # [1012, 1201)
E2 = "GCC" * 10 + "TAA"                            # exon 2: [1201, 1234), ends on a stop
gene = E1 + INTRON + E2
G[1000:1000 + len(gene)] = list(gene)
SEQ = "".join(G)
SEQM = rc(SEQ)                                     # the same gene on the minus strand
D, A, END = 1012, 1201, 1234                       # true donor, acceptor, 3' end


def cds(name, ex, strand=1):
    if strand == -1:
        ex = sorted((L - e, L - s) for s, e in ex)
    f = Feature(gene_name=name, gene_type="CDS", start=ex[0][0], end=ex[-1][1], strand=strand,
                engine="A", flag="HIGH")
    f.exons = list(ex)
    return f


# the panel: exon 1 from 1003 to the true donor, exon 2 from the true acceptor to 3 bp before the
# stop -- right junction, different ends. Coordinates are in the CDS window, which starts PAD bp
# before the CDS's 5' end, read in coding orientation; w0 is that window's start on the plus copy.
HIT = {}
RS._blast_exons = lambda Cseq, gene: {0: (1003 - HIT["w0"], D - HIT["w0"]), 1: (A - HIT["w0"], END - 3 - HIT["w0"])}
RS._meta = lambda: {"tK": {"n_exons": 2, "keep_ends": True}, "tN": {"n_exons": 2}}
RS._tpl = lambda gene: None


def refine(f, seq=SEQ, five=1000):
    """five: the input CDS's 5' end on the plus copy, which fixes where its window starts."""
    HIT["w0"] = five - RS.PAD
    g = copy.deepcopy(f)
    return RS._refine_one(g, seq), g


off = [(1000, D + 3), (A + 3, END)]                # engine junction 3 bp late on both sides
print("--- keep_ends: the junction follows the panel, the ends do not ---")
ch, f = refine(cds("tK", off))
check("plus: exons are the input's ends with the panel's junction",
      (ch, f.exons), (True, [(1000, D), (A, END)]))
check("... so the CDS still starts on its ATG and ends on its stop",
      (C.extract(SEQ, f, L)[:3], C.extract(SEQ, f, L)[-3:]), ("ATG", "TAA"))
ch, f = refine(cds("tK", off, strand=-1), seq=SEQM)
check("minus: the same", (ch, sorted(f.exons)), (True, sorted([(L - END, L - A), (L - D, L - 1000)])))
check("... reading ATG ... TAA on the minus strand",
      (C.extract(SEQM, f, L)[:3], C.extract(SEQM, f, L)[-3:]), ("ATG", "TAA"))
bad_start = [(1003, D + 3), (A + 3, END)]          # the input starts on GCC, not a start codon
ch, f = refine(cds("tK", bad_start), five=1003)
check("an input start that is not a start codon is kept for the start pass to handle",
      (ch, f.exons), (True, [(1003, D), (A, END)]))
ch, f = refine(cds("tK", [(1000, D), (A, END)]))
check("an input whose junction is already the panel's is left as it is", (ch, f.exons), (False, [(1000, D), (A, END)]))

print("--- without keep_ends the panel's ends are transferred, as before ---")
ch, f = refine(cds("tN", off))
check("plus: exons are the panel's", (ch, f.exons), (True, [(1003, D), (A, END - 3)]))

print("\n%d checks, %d failed" % (RUN[0], len(FAIL)))
sys.exit(1 if FAIL else 0)
