#!/usr/bin/env python3
"""
Plastanno v2 — build_all.py
Build ALL databases from scratch.

Usage:
    python3 scripts/build/build_all.py \
        --genes   /path/to/all_genes_full.tsv \
        --train   /path/to/splits/train.csv \
        --rawdata /path/to/rawdata/ \
        --outdir  /path/to/database/

Databases built:
    1. trna_db/      Hierarchical tRNA (genus/family/global)
    2. exon_db/      Actual exon sequences for intron tRNAs
    3. rrna_db/      Full-length rRNA per gene
    4. blast_db/     Genus representative genomes
    5. gene_catalog.json
"""

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────

STANDARD_TRNA = {
    "trnA-UGC", "trnC-GCA", "trnD-GUC", "trnE-UUC",
    "trnF-GAA", "trnfM-CAU","trnG-GCC", "trnG-UCC",
    "trnH-GUG", "trnI-CAU", "trnI-GAU", "trnK-UUU",
    "trnL-CAA", "trnL-UAA", "trnL-UAG", "trnM-CAU",
    "trnN-GUU", "trnP-UGG", "trnQ-UUG", "trnR-ACG",
    "trnR-UCU", "trnS-GCU", "trnS-GGA", "trnS-UGA",
    "trnT-GGU", "trnT-UGU", "trnV-GAC", "trnV-UAC",
    "trnW-CCA", "trnY-GUA",
}

INTRON_TRNA = {
    "trnK-UUU", "trnG-UCC", "trnL-UAA", "trnV-UAC",
    "trnI-GAU", "trnA-UGC", "trnG-GCC",
}

RRNA_LENGTH = {
    "rrn16" : (1400, 1600),
    "rrn23" : (2600, 3000),
    "rrn4.5": (90,   120),
    "rrn5"  : (110,  135),
}

# ──────────────────────────────────────────────────────────────────────────────
# Utilities
# ──────────────────────────────────────────────────────────────────────────────

def log(msg):
    print(msg, flush=True)

def make_record(seq, rec_id, desc=""):
    seq = str(seq).upper().replace("U", "T")
    return SeqRecord(Seq(seq), id=rec_id, description=desc)

def build_blast_db(fasta_path):
    """Build BLAST nucleotide DB version 4."""
    out = str(fasta_path).replace(".fasta", "")
    r = subprocess.run([
        "makeblastdb",
        "-in",              str(fasta_path),
        "-dbtype",          "nucl",
        "-out",             out,
        "-blastdb_version", "4",
    ], capture_output=True, text=True)
    if r.returncode != 0:
        log(f"    ❌ makeblastdb failed: {r.stderr.strip()}")
        return False
    return True

def write_fasta_and_build(records, fasta_path):
    """Write records to FASTA and build BLAST DB."""
    if not records:
        log(f"    ⚠️  No records for {fasta_path.name}")
        return False
    SeqIO.write(records, fasta_path, "fasta")
    return build_blast_db(fasta_path)

# ──────────────────────────────────────────────────────────────────────────────
# Step 1: tRNA hierarchical database
# ──────────────────────────────────────────────────────────────────────────────

def build_trna_db(df, out_dir):
    """
    Build hierarchical tRNA database:
        genus/  (≥3 genomes per genus)
        family/ (≥5 genomes per family)
        global  (fallback for all)
    """
    log("\n[1/4] Building hierarchical tRNA database...")
    (out_dir / "genus").mkdir(parents=True, exist_ok=True)
    (out_dir / "family").mkdir(parents=True, exist_ok=True)

    # Filter: standard tRNA, valid sequence, reasonable length
    trna = df[
        (df["gene_type"]  == "tRNA") &
        (df["sequence"].notna()) &
        (df["length_bp"]  >= 60) &
        (df["length_bp"]  <= 3000) &
        (df["gene_name"].isin(STANDARD_TRNA)) &
        (df["family"].notna()) &
        (df["genus"].notna())
    ].copy()
    log(f"  Standard tRNA entries: {len(trna):,}")

    def make_recs(subset, max_per_gene):
        recs = []
        for gene, grp in subset.groupby("gene_name"):
            for _, row in grp.head(max_per_gene).iterrows():
                recs.append(make_record(
                    row["sequence"],
                    f"{row['accession']}_{gene}",
                    gene,
                ))
        return recs

    # Global
    global_recs = make_recs(trna, max_per_gene=20)
    write_fasta_and_build(global_recs, out_dir / "global.fasta")
    log(f"  Global: {len(global_recs)} sequences")

    # Genus (≥3 genomes)
    genus_built = []
    counts = trna.groupby("genus")["accession"].nunique()
    for genus, n in counts.items():
        if n < 3 or not isinstance(genus, str): continue
        recs = make_recs(trna[trna["genus"] == genus], max_per_gene=5)
        if write_fasta_and_build(recs, out_dir / "genus" / f"{genus}.fasta"):
            genus_built.append(genus)
    log(f"  Genus DBs: {len(genus_built)}")

    # Family (≥5 genomes)
    family_built = []
    counts = trna.groupby("family")["accession"].nunique()
    for family, n in counts.items():
        if n < 5 or not isinstance(family, str): continue
        recs = make_recs(trna[trna["family"] == family], max_per_gene=10)
        if write_fasta_and_build(recs, out_dir / "family" / f"{family}.fasta"):
            family_built.append(family)
    log(f"  Family DBs: {len(family_built)}")

    # Index
    with open(out_dir / "index.json", "w") as f:
        json.dump({
            "genus_dbs" : genus_built,
            "family_dbs": family_built,
        }, f, indent=2)
    log(f"  ✅ tRNA DB complete")

# ──────────────────────────────────────────────────────────────────────────────
# Step 2: tRNA exon database (actual exon sequences from GenBank)
# ──────────────────────────────────────────────────────────────────────────────

def build_exon_db(train_csv, rawdata_dir, out_dir):
    """
    Extract ACTUAL exon sequences from GenBank feature locations.
    This is the correct approach — NOT midpoint splitting.

    For each intron-containing tRNA in GenBank:
        feat.location.parts[0] → exon1 sequence
        feat.location.parts[1] → exon2 sequence
    """
    log("\n[2/4] Building tRNA exon database (actual exon sequences)...")
    out_dir.mkdir(parents=True, exist_ok=True)

    train_df  = pd.read_csv(train_csv)
    rawdata   = Path(rawdata_dir)
    exon_recs = []
    processed = 0
    skipped   = 0

    for _, row in train_df.iterrows():
        acc      = row["accession"].split(".")[0]
        gb_files = list(rawdata.glob(f"{acc}*.gb"))
        if not gb_files:
            skipped += 1
            continue

        try:
            record = SeqIO.read(gb_files[0], "genbank")
            for feat in record.features:
                if feat.type != "tRNA": continue
                gene = feat.qualifiers.get("gene", [""])[0]
                if gene not in INTRON_TRNA: continue
                if len(feat.location.parts) < 2: continue

                # Extract each exon directly from genome coordinates
                for ei, part in enumerate(feat.location.parts):
                    exon_seq = str(part.extract(record.seq))
                    if len(exon_seq) < 20: continue
                    exon_recs.append(make_record(
                        exon_seq,
                        f"{acc}_{gene}_e{ei+1}",
                        gene,
                    ))
            processed += 1
        except Exception:
            skipped += 1
            continue

        if processed % 1000 == 0:
            log(f"    Processed: {processed:,} genomes, "
                f"{len(exon_recs):,} exon sequences")

    log(f"  Processed: {processed:,} genomes")
    log(f"  Skipped  : {skipped:,}")
    log(f"  Exon seqs: {len(exon_recs):,}")

    write_fasta_and_build(exon_recs, out_dir / "trna_exon.fasta")
    log(f"  ✅ Exon DB complete")

# ──────────────────────────────────────────────────────────────────────────────
# Step 3: rRNA database (full-length sequences only)
# ──────────────────────────────────────────────────────────────────────────────

def build_rrna_db(df, out_dir):
    """
    Build full-length rRNA databases.
    Strict length filter to avoid fragments.
    """
    log("\n[3/4] Building rRNA databases...")
    out_dir.mkdir(parents=True, exist_ok=True)

    rrna = df[df["gene_type"] == "rRNA"]

    for gene, (min_len, max_len) in RRNA_LENGTH.items():
        subset = rrna[
            (rrna["gene_name"] == gene) &
            (rrna["sequence"].notna()) &
            (rrna["length_bp"] >= min_len) &
            (rrna["length_bp"] <= max_len)
        ]
        recs = [
            make_record(
                row["sequence"],
                f"{row['accession']}_{gene}",
                gene,
            )
            for _, row in subset.head(3000).iterrows()
        ]
        fname = gene.replace(".", "") + ".fasta"
        write_fasta_and_build(recs, out_dir / fname)
        log(f"  {gene:<8}: {len(recs):>5} sequences "
            f"({min_len}-{max_len}bp)")

    log(f"  ✅ rRNA DB complete")

# ──────────────────────────────────────────────────────────────────────────────
# Step 4: Gene catalog
# ──────────────────────────────────────────────────────────────────────────────

def build_gene_catalog(df, out_path):
    """
    Build gene catalog with region, type, product, multi-exon info.
    """
    log("\n[4/4] Building gene catalog...")

    # Known regions from v1 experience
    IR_GENES = {
        "rpl2","rpl23","ndhB","rps7","ycf2",
        "orf70","trnI-CAU","trnI-GAU","trnA-UGC",
        "trnR-ACG","trnN-GUU","trnL-CAA","trnV-GAC",
        "rrn16","rrn23","rrn4.5","rrn5",
    }
    SSC_GENES = {
        "ndhF","rpl32","ccsA","ndhD","psaC","ndhE",
        "ndhG","ndhI","ndhA","ndhH","rps15","ycf1",
        "trnL-UAG",
    }

    # Multi-exon CDS from v1 benchmark
    MULTIEXON_CDS = {
        "rps16" : 2, "atpF"  : 2, "rpoC1" : 2,
        "ycf3"  : 3, "rps12" : 3, "clpP"  : 3,
        "petB"  : 2, "petD"  : 2, "rpl16" : 2,
        "rpl2"  : 2, "ndhB"  : 2, "ndhA"  : 2,
    }

    # Short first exons (v1 finding)
    SHORT_FIRST_EXON = {
        "petB": 6,
        "petD": 9,
        "rpl16": 9,
    }

    # Gene synonyms from v1
    SYNONYMS = {
        "clpP1" : "clpP",
        "clpP2" : "clpP",
        "psbN"  : "pbf1",
        "orf70a": "orf70",
        "orf70b": "orf70",
    }

    # Build catalog from training data
    catalog = {}
    for gene in df["gene_name"].unique():
        subset = df[df["gene_name"] == gene]
        if len(subset) == 0: continue

        gtype = subset["gene_type"].mode()[0]

        # Determine region
        if gene in IR_GENES or gtype == "rRNA":
            region = "IRb"
        elif gene in SSC_GENES:
            region = "SSC"
        else:
            region = "LSC"

        catalog[gene] = {
            "type"             : gtype,
            "region"           : region,
            "n_exons"          : MULTIEXON_CDS.get(gene, 1),
            "short_first_exon" : SHORT_FIRST_EXON.get(gene, 0),
            "synonym_of"       : None,
            "product"          : (
                                     subset["product"].dropna().mode().iloc[0]
                                     if "product" in subset.columns
                                     and not subset["product"].dropna().empty
                                     else gene
                                 ),
        }

    # Add reverse synonyms
    for alias, canonical in SYNONYMS.items():
        if canonical in catalog:
            catalog[alias] = {**catalog[canonical],
                              "synonym_of": canonical}

    with open(out_path, "w") as f:
        json.dump(catalog, f, indent=2)

    log(f"  Genes in catalog: {len(catalog):,}")
    log(f"  Multi-exon CDS  : {len(MULTIEXON_CDS)}")
    log(f"  Synonyms        : {len(SYNONYMS)}")
    log(f"  ✅ Gene catalog complete")

# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plastanno v2 — Build all databases"
    )
    parser.add_argument("--genes",   required=True,
                        help="Path to all_genes_full.tsv")
    parser.add_argument("--train",   required=True,
                        help="Path to splits/train.csv")
    parser.add_argument("--rawdata", required=True,
                        help="Path to rawdata directory (GenBank files)")
    parser.add_argument("--outdir",  required=True,
                        help="Output directory for databases")
    args = parser.parse_args()

    OUT = Path(args.outdir)
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Plastanno v2 — Database Builder")
    print("=" * 60)

    # Load main data
    log("\nLoading all_genes_full.tsv...")
    df = pd.read_csv(args.genes, sep="\t", low_memory=False)
    log(f"  Total entries: {len(df):,}")
    log(f"  Gene types   : {df['gene_type'].value_counts().to_dict()}")

    # Build databases
    build_trna_db(df,  OUT / "trna_db")
    build_exon_db(args.train, args.rawdata, OUT / "exon_db")
    build_rrna_db(df,  OUT / "rrna_db")
    build_gene_catalog(df, OUT / "gene_catalog.json")

    print("\n" + "=" * 60)
    print("✅ All databases built successfully!")
    print(f"   Output: {OUT}")
    print("=" * 60)

if __name__ == "__main__":
    main()
