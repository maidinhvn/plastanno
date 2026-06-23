#!/usr/bin/env python3
"""
Build protein database for Exonerate (protein2genome).

Approach: Translate CDS nucleotide sequences from all_genes_full.tsv
- Fast: in-memory translation, no file I/O per genome
- Portable: no rawdata directory needed
- Diverse: 1 sequence per order for maximum coverage
"""
import argparse
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ── Constants ─────────────────────────────────────────────────────────────────
STANDARD_CDS = {
    "accD","atpA","atpB","atpE","atpF","atpH","atpI",
    "ccsA","cemA","clpP","infA","matK",
    "ndhA","ndhB","ndhC","ndhD","ndhE","ndhF",
    "ndhG","ndhH","ndhI","ndhJ","ndhK",
    "petA","petB","petD","petG","petL","petN",
    "psaA","psaB","psaC","psaI","psaJ",
    "psbA","psbB","psbC","psbD","psbE","psbF",
    "psbH","psbI","psbJ","psbK","psbL","psbM",
    "pbf1","psbT","psbZ",
    "rbcL","rpl2","rpl14","rpl16","rpl20","rpl22",
    "rpl23","rpl32","rpl33","rpl36",
    "rpoA","rpoB","rpoC1","rpoC2",
    "rps2","rps3","rps4","rps7","rps8","rps11",
    "rps12","rps14","rps15","rps16","rps18","rps19",
    "ycf1","ycf2","ycf3","ycf4",
}

SYNONYMS = {
    "clpP1":"clpP", "clpP2":"clpP",
    "psbN" :"pbf1",
}

MAX_PER_GENE = 500  # max sequences per gene
MIN_AA_LEN   = 30   # minimum protein length

# ── Utilities ─────────────────────────────────────────────────────────────────

def translate_cds(nuc_seq):
    """
    Translate nucleotide to amino acid.
    Returns AA string or None if invalid.
    """
    seq = str(nuc_seq).upper().replace("U","T")
    # Pad to multiple of 3
    remainder = len(seq) % 3
    if remainder:
        seq = seq + "N" * (3 - remainder)
    try:
        aa = str(Seq(seq).translate(to_stop=True))
        if len(aa) >= MIN_AA_LEN:
            return aa
    except Exception:
        pass
    return None

def build_blast_prot_db(fasta_path):
    r = subprocess.run([
        "makeblastdb",
        "-in",              str(fasta_path),
        "-dbtype",          "prot",
        "-out",             str(fasta_path).replace(".fasta",""),
        "-blastdb_version", "4",
    ], capture_output=True, text=True)
    return r.returncode == 0

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes",  required=True,
                        help="Path to all_genes_full.tsv")
    parser.add_argument("--train",  required=True,
                        help="Path to splits/train.csv")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    args = parser.parse_args()

    OUT = Path(args.outdir)
    OUT.mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("Building protein database (translate from TSV)")
    print("="*60)

    # Load data
    print("\nLoading all_genes_full.tsv...")
    df    = pd.read_csv(args.genes, sep="\t", low_memory=False)
    train = pd.read_csv(args.train)

    # Filter CDS
    cds = df[
        (df["gene_type"] == "CDS") &
        (df["sequence"].notna()) &
        (df["length_bp"] >= 90) &
        (df["gene_name"].apply(
            lambda g: SYNONYMS.get(g,g) in STANDARD_CDS
        ))
    ].copy()

    # Normalize gene names
    cds["gene_name"] = cds["gene_name"].apply(
        lambda g: SYNONYMS.get(g, g)
    )

    print(f"  CDS entries : {len(cds):,}")
    print(f"  Unique genes: {cds['gene_name'].nunique()}")

    # order column already in all_genes_full.tsv
    if "order" not in cds.columns:
        cds["order"] = "unknown"
    cds["order"] = cds["order"].fillna("unknown")

    print("\nTranslating and building databases...")
    built = 0
    total_seqs = 0

    for gene, group in cds.groupby("gene_name"):
        # 1 sequence per order → diversity
        seen_orders = {}
        for _, row in group.iterrows():
            order = row["order"]
            if order in seen_orders: continue
            aa = translate_cds(row["sequence"])
            if aa:
                seen_orders[order] = (row["accession"], aa)
            if len(seen_orders) >= MAX_PER_GENE:
                break

        if len(seen_orders) < 3:
            continue  # skip genes with too few sequences

        # Write FASTA
        records = [
            SeqRecord(
                Seq(aa),
                id          = f"{acc}_{gene}",
                description = gene,
            )
            for order, (acc, aa) in seen_orders.items()
        ]

        out_fa = OUT / f"{gene}.fasta"
        SeqIO.write(records, out_fa, "fasta")

        if build_blast_prot_db(out_fa):
            built      += 1
            total_seqs += len(records)
            print(f"  ✅ {gene:<15}: {len(records):>4} sequences "
                  f"({len(seen_orders)} orders)")
        else:
            print(f"  ❌ {gene}: failed")

    print(f"\n{'='*60}")
    print(f"✅ Done!")
    print(f"   Genes built : {built}")
    print(f"   Total seqs  : {total_seqs:,}")
    print(f"   Output      : {OUT}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
