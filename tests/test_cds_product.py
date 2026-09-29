#!/usr/bin/env python3
"""A CDS is written with its protein name as /product, never its gene symbol.

Run DIRECTLY and check $?. Exit 0 = all pass.

Until this fix no CDS but rps12 reached the writers with a product, and each of
the three writers substituted the gene symbol: /product="psbA" instead of
"photosystem II protein D1" in 97.7% of CDS, in the .gb, the .gff3 and the .tbl
meant for NCBI submission (parked/product_qualifier/FINDINGS.md). References use
the bare symbol in 15 of 33,222 CDS. The catalog held the protein names all along.

The product is now set once, in finalize.assign_cds_products, before any writer
runs, so the three formats cannot disagree.
"""
import inspect
import json
import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from plastanno.core.feature import Feature                     # noqa: E402
from plastanno.core.finalize import assign_cds_products        # noqa: E402
from plastanno.output.writers import write_genbank, write_gff3, write_tbl  # noqa: E402
from plastanno import paths                                    # noqa: E402

FAIL, RUN = [], [0]


def check(label, got, want):
    RUN[0] += 1
    ok = got == want
    print("  %-66s %s" % (label, "ok" if ok else "FAIL got=%r want=%r" % (got, want)))
    if not ok:
        FAIL.append(label)


CAT = json.load(open(os.path.join(str(paths.config_dir()), "gene_catalog.json")))


def ycf3():
    """A three-exon CDS: .tbl writes its 2nd and 3rd intervals as bare lines."""
    f = Feature(gene_name="ycf3", gene_type="CDS", start=4100, end=4900, strand=1, engine="A")
    f.exons = [(4100, 4224), (4400, 4628), (4747, 4900)]
    return f


def feats():
    def F(name, s, e, gtype="CDS", **kw):
        f = Feature(gene_name=name, gene_type=gtype, start=s, end=e, strand=1,
                    engine="B", **kw)
        f.exons = [(s, e)]
        return f
    return [F("psbA", 100, 1162),
            F("rps12", 2000, 2372, product="ribosomal protein S12"),
            F("orfX", 3000, 3300),
            F("atpA", 3400, 4000, product="set by an earlier step"),
            ycf3(),
            F("trnH-GUG", 4000, 4075, gtype="tRNA", product="tRNA-His"),
            F("rrn16", 5000, 6491, gtype="rRNA", product="16S ribosomal RNA")]


print("--- assign_cds_products ---")
fs = feats()
named, unnamed = assign_cds_products(fs, CAT)
by = {f.gene_name: f for f in fs}
check("psbA gets its protein name", by["psbA"].product, "photosystem II protein D1")
check("a product already set is kept (rps12)", by["rps12"].product, "ribosomal protein S12")
check("a product that differs from the catalog is not overwritten",
      by["atpA"].product, "set by an earlier step")
check("tRNA and rRNA products are untouched",
      (by["trnH-GUG"].product, by["rrn16"].product), ("tRNA-His", "16S ribosomal RNA"))
check("a name the catalog lacks keeps no product", by["orfX"].product, "")
check("  ... and says so in its notes", any("no protein name" in n for n in by["orfX"].notes), True)
check("counts: 2 named, 1 unnamed", (named, unnamed), (2, 1))
pseudo = Feature(gene_name="ycf1", gene_type="CDS", start=1, end=500, strand=1, engine="A")
pseudo.is_pseudogene = True
assign_cds_products([pseudo], CAT)
check("a pseudogene is skipped (it gets no product feature)", pseudo.product, "")

print("--- the three writers agree ---")
G = "".join("ACGT"[(i * 7) % 4] for i in range(7000))
fs = feats()
assign_cds_products(fs, CAT)
d = tempfile.mkdtemp(prefix="prod_")
write_genbank(fs, G, "TEST.1", len(G), {"LSC": (0, len(G))}, [], os.path.join(d, "t.gb"))
write_gff3(fs, len(G), "TEST.1", os.path.join(d, "t.gff3"))
write_tbl(fs, "TEST.1", os.path.join(d, "t.tbl"), genome_len=len(G))
def parse_gb(path):
    """gene -> product for every CDS; continuation lines start in column 22."""
    t = re.sub(r"\n {21}(?=[^/\s])", " ", open(path).read())
    out = {}
    for blk in re.split(r"\n     (?=\S)", t):
        if blk.startswith("CDS "):
            g = re.search(r'/gene="([^"]+)"', blk)
            pr = re.search(r'/product="([^"]+)"', blk)
            out[g.group(1)] = pr.group(1) if pr else None
    return out


def parse_gff3(path):
    out = {}
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if len(c) == 9 and c[2] == "CDS":
            a = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            out[a["ID"][4:].rsplit("_", 1)[0]] = a.get("product")
    return out


def parse_tbl(path):
    """A feature's 2nd and later intervals are bare 'start<TAB>end' lines; a first
    version of this parser read them as new features and dropped every
    multi-exon CDS, which the single-exon fixture could not notice."""
    out, cur = {}, None
    def flush():
        if cur is not None and "gene" in cur:
            out[cur["gene"]] = cur.get("product")
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("\t\t\t"):
            if cur is not None:
                k, _, v = line.strip("\t").partition("\t")
                cur.setdefault(k, v)
        elif re.match(r"^[<>]?\d+\t[<>]?\d+\t\S", line):
            flush()
            cur = {} if line.split("\t")[2] == "CDS" else None
        elif not re.match(r"^[<>]?\d+\t[<>]?\d+$", line):
            flush()
            cur = None
    flush()
    return out


gb_cds = parse_gb(os.path.join(d, "t.gb"))
gff_cds = parse_gff3(os.path.join(d, "t.gff3"))
tbl_cds = parse_tbl(os.path.join(d, "t.tbl"))
for fmt, got in (("gb", gb_cds), ("gff3", gff_cds), ("tbl", tbl_cds)):
    check(".%s: the parser finds all five CDS (not a vacuous pass)" % fmt,
          sorted(got), ["atpA", "orfX", "psbA", "rps12", "ycf3"])
for fmt, got in (("gb", gb_cds), ("gff3", gff_cds), ("tbl", tbl_cds)):
    check(".%s: psbA /product is the protein name" % fmt, got.get("psbA"),
          "photosystem II protein D1")
for fmt, got in (("gb", gb_cds), ("gff3", gff_cds), ("tbl", tbl_cds)):
    check(".%s: the three-exon ycf3 carries its protein name" % fmt, got.get("ycf3"),
          "photosystem I assembly protein Ycf3")
check("no format writes a known CDS with its symbol as product",
      [(f, g) for f, got in (("gb", gb_cds), ("gff3", gff_cds), ("tbl", tbl_cds))
       for g, p in got.items() if p == g and g != "orfX"], [])

print("--- the catalog can name every CDS the tool emits ---")
EMITTED = ("accD atpA atpB atpE atpF atpH atpI ccsA cemA chlB chlL chlN clpP infA matK "
           "ndhA ndhB ndhC ndhD ndhE ndhF ndhG ndhH ndhI ndhJ ndhK petA petB petD petG "
           "petL petN psaA psaB psaC psaI psaJ psbA psbB psbC psbD psbE psbF psbH psbI "
           "psbJ psbK psbL psbM psbN psbT psbZ rbcL rpl14 rpl16 rpl2 rpl20 rpl21 rpl22 "
           "rpl23 rpl32 rpl33 rpl36 rpoA rpoB rpoC1 rpoC2 rps11 rps12 rps14 rps15 rps16 "
           "rps18 rps19 rps2 rps3 rps4 rps7 rps8 ycf1 ycf2 ycf3 ycf4 ycf66").split()
check("every emitted CDS name has a product other than its symbol",
      [g for g in EMITTED if not (CAT.get(g) or {}).get("product")
       or CAT[g]["product"] == g], [])
check("chlN is spelled protochlorophyllide, as chlB and chlL are",
      CAT["chlN"]["product"], "protochlorophyllide reductase ChlN subunit")

print("--- pipeline.run assigns products after QC and before writing ---")
import plastanno.pipeline as PL                                # noqa: E402
src = inspect.getsource(PL.run)
i_qc, i_ap, i_w = (src.find("finalize_qc(annotations"), src.find("assign_cds_products("),
                   src.find("write_all("))
check("assign_cds_products is called", i_ap >= 0, True)
check("  ... after finalize_qc and before write_all", i_qc < i_ap < i_w, True)

print()
print("%d checks, %d failed" % (RUN[0], len(FAIL)))
for f in FAIL:
    print("  FAILED: %s" % f)
sys.exit(1 if FAIL else 0)
