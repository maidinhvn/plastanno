#!/usr/bin/env python3
"""A gene that an expanded inverted repeat has taken out of its catalog region is still found,
and what is found there is marked for review.

Run DIRECTLY and check $?. Exit 0 = all pass.

Engine A searches a gene outside the IR gene set only in its catalog region. In genomes whose
inverted repeat has expanded over most of the SSC, SSC genes sit in the two IR copies, and
ndhA, the one SSC gene with an intron, was lost whole (Engine B's profile hits cover one exon
each). Now, when the catalog region gives no hit reaching 0.6x the gene's expected length, the
gene is searched in IRb and IRa too:

  - "no hit at all" is not the trigger: the search window around the SSC can reach a fragment of
    an IR copy near the boundary (Cyperaceae ndhA), and that fragment used to block the IR search;
  - IR hits replace the region hits they overlap;
  - each carries an "[IR search]" note, and finalize.submission_check keeps such a CDS
    NEEDS_REVIEW (in Knoxia the IR search found a ycf1 its reference does not annotate at all).

Exonerate is replaced by a fake that "finds" the gene wherever its copies lie inside the region
searched, and records every region it was asked about.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.core.feature import Feature                                   # noqa: E402
from plastanno.core.finalize import submission_check, SUBMIT, IR_SEARCH  # noqa: E402
from plastanno.identify import engine_a as EA                                # noqa: E402
from plastanno.output.writers import write_report                           # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-76s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


L = 200000
G = list("GCC" * (L // 3 + 1))[:L]           # no stop codon in any frame on either strand
G[124000:124003] = list("ATG"); G[126097:126100] = list("TAA")      # ndhA, IRb copy (+)
G[145097:145100] = list("CAT"); G[143000:143003] = list("TTA")      # ndhA, IRa copy (-)
SEQ = "".join(G)
BUF = 2000            # Exonerate is run on the region plus this much either side (engine_a.BUFFER)
DB = tempfile.mkdtemp(prefix="irfb_db_")
for g in ("ndhA", "ndhC", "ndhF", "rpl2", "psbA", "orfX"):
    with open(os.path.join(DB, g + ".fasta"), "w") as fh:
        fh.write(">X_%s\nMSTARTPROTEIN\n" % g)
CAT = {"ndhA": {"region": "SSC", "expected_len": 2000}, "ndhC": {"region": "SSC", "expected_len": 2000},
       "ndhF": {"region": "SSC", "expected_len": 2000}, "rpl2": {"region": "IRb"},
       "psbA": {"region": "LSC", "expected_len": 1000}, "orfX": {"region": "SSC"}}
# an IR that has taken almost the whole SSC: SSC is 3 kb
IRS = {"LSC": (0, 70000), "IRb": (70000, 133000), "SSC": (133000, 136000), "IRa": (136000, 199000)}
COPIES = {}           # gene -> [(start, end, strand)] where the fake Exonerate "finds" it
ASKED = []            # every region searched: (gene, start, end)


def fake_exonerate(genome_seq, protein_seq, gene_name, prot_acc, region_start, region_end, genome_len):
    """Like Exonerate on the region plus its window: a copy inside the window is found whole, a
    copy the window cuts is found as the fragment inside it."""
    ASKED.append((gene_name, region_start, region_end))
    ws, we = region_start - BUF, region_end + BUF
    out = []
    for s, e, st in COPIES.get(gene_name, []):
        a, b = max(s, ws), min(e, we)
        if b - a >= 100:
            f = Feature(gene_name=gene_name, gene_type="CDS", start=a, end=b, strand=st, engine="A", flag="HIGH")
            f.exons = [(a, b)]
            f.s_ref = 0.9
            out.append(f)
    return out


EA.run_exonerate_region = fake_exonerate


def run(gene, irb=IRS):
    ASKED.clear()
    return EA.run_exonerate_gene(SEQ, gene, DB, irb, CAT, threads=1)


note = lambda f: [n for n in f.notes if str(n).startswith(IR_SEARCH)]
spans = lambda fs: sorted((f.start, f.end, f.strand) for f in fs)
asked = lambda: sorted({a[1:] for a in ASKED})

print("--- a gene the IR has taken ---")
COPIES["ndhA"] = [(124000, 126100, 1), (143000, 145100, -1)]      # one full copy in each IR
feats = run("ndhA")
check("both IR copies are found", spans(feats), [(124000, 126100, 1), (143000, 145100, -1)])
check("... each with the note", [note(f) for f in feats],
      [[IR_SEARCH + "found by searching the inverted repeats, because its catalog region (SSC) "
        "gave no complete hit"]] * 2)
check("... after the catalog region was searched first", [a[1:] for a in ASKED][:1], [IRS["SSC"]])

print("--- a fragment of an IR copy inside the SSC window does not block the search ---")
# the IR copies sit 1.3 kb from the SSC: its search window reaches 700 bp of each (Cyperaceae)
COPIES["ndhC"] = [(129500, 131700, 1), (137300, 139500, -1)]
check("the SSC search alone sees only 700 bp fragments",
      spans(EA.run_exonerate_gene(SEQ, "ndhC", DB, {"SSC": IRS["SSC"]}, CAT, threads=1)),
      [(131000, 131700, 1), (137300, 138000, -1)])
feats = run("ndhC")
check("the IRs are searched although the SSC gave hits", IRS["IRb"] in asked(), True)
check("the full copies are found and replace the fragments they overlap",
      spans(feats), [(129500, 131700, 1), (137300, 139500, -1)])
COPIES["ndhF"] = [(133200, 135400, -1)]                            # 2200 bp: a usable SSC hit
feats = run("ndhF")
check("a usable hit in the catalog region: no IR search, no note",
      (spans(feats), asked(), [note(f) for f in feats]), ([(133200, 135400, -1)], [IRS["SSC"]], [[]]))
COPIES["orfX"] = [(133100, 133300, 1), (124000, 126000, 1)]          # no expected length: v1 rule
feats = run("orfX")
check("a gene with no expected length keeps the v1 rule (any hit blocks the IR search)",
      (spans(feats), asked()), ([(133100, 133300, 1)], [IRS["SSC"]]))

print("--- nothing changes elsewhere ---")
COPIES["rpl2"] = [(80000, 81500, -1), (190000, 191500, 1)]
feats = run("rpl2")
check("an IR gene is searched in both IR copies as before, without the note",
      (len(feats), asked(), [note(f) for f in feats]), (2, [IRS["IRb"], IRS["IRa"]], [[], []]))
COPIES["psbA"] = []
feats = run("psbA")
check("a gene found nowhere: the IRs are searched, nothing is returned",
      (feats, asked()), ([], sorted({IRS["LSC"], IRS["IRb"], IRS["IRa"]})))
feats = run("ndhA", irb={"LSC": (0, L)})
check("an IR-lacking genome searches the whole sequence, as before, with no note",
      (len(feats), asked(), [note(f) for f in feats]), (2, [(0, L)], [[], []]))
feats = run("ndhA", irb={"LSC": (0, 70000), "IRb": (70000, 133000), "SSC": (133000, 136000)})
check("only one IR copy detected: no IR search", (feats, asked()), ([], [IRS["SSC"]]))

print("--- the review policy ---")
cds = run("ndhA")
anns = cds
bad = submission_check(anns, SEQ)
check("both are valid ORFs: nothing NCBI would reject", [p for _, p in bad],
      [["found by searching the inverted repeats"]] * 2)
check("submission_check keeps both IR-found CDS NEEDS_REVIEW", [a.flag for a in cds], ["NEEDS_REVIEW"] * 2)
sub = [str(n) for n in cds[0].notes if str(n).startswith(SUBMIT)]
check("... with a note saying where it was found",
      any("it was found by searching the inverted repeats, because its catalog region (SSC) gave no "
          "complete hit" in n for n in sub), True)
with tempfile.TemporaryDirectory() as d:
    rep = os.path.join(d, "r.report")
    write_report(anns, "TEST", L, {}, [], 1.0, rep)
    txt = open(rep).read()
check("the report lists them apart", "pass, but were found by searching the inverted repeats" in txt, True)

print("\n%d checks, %d failed" % (RUN[0], len(FAIL)))
sys.exit(1 if FAIL else 0)
