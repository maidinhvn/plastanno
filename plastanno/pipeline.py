#!/usr/bin/env python3
"""
Plastanno v2 — Main pipeline
Hybrid annotation: Engine A (reference) + Engine B (model)
with reconciliation layer.
"""
import time
from pathlib import Path
from typing import Optional

from .core.feature  import Feature
from .core.reconcile import reconcile


# ── Paths ─────────────────────────────────────────────────────────────────────
from .paths import db_root, config_dir
BASE_DIR    = db_root()
BLAST_DB    = str(BASE_DIR / "blast_db"    / "genus_reps")
PROTEIN_DB  = str(BASE_DIR / "protein_db")
TRNA_DB_DIR = BASE_DIR / "trna_db"
EXON_DB     = str(BASE_DIR / "exon_db"    / "trna_exon")
RRNA_DB_DIR = BASE_DIR / "rrna_db"
HMM_DIR     = BASE_DIR / "hmm_db"         / "hmm_profiles"
GENE_CAT    = str(config_dir() / "gene_catalog.json")

RRNA_DBS = {
    "rrn16" : str(RRNA_DB_DIR / "rrn16"),
    "rrn23" : str(RRNA_DB_DIR / "rrn23"),
    "rrn4.5": str(RRNA_DB_DIR / "rrn45"),
    "rrn5"  : str(RRNA_DB_DIR / "rrn5"),
}


# ── Pipeline ──────────────────────────────────────────────────────────────────

def run(
    input_fasta : str,
    output_dir  : str,
    prefix      : Optional[str] = None,
    threads     : int = 4,
    no_plot     : bool = False,
    reference   : Optional[str] = None,
    use_trnascan: bool = False,
    verbose     : bool = False,
) -> dict:
    """
    Main annotation pipeline.

    Steps:
      1. Read FASTA
      2. IR boundary detection
      3. Closest relative finding
      4A. Engine A: Exonerate CDS (reference-based)
      4B. Engine B: HMM scan + tRNA + rRNA (model-based)
      5. Reconciliation
      6. Special cases
      7. Output
    """
    import json
    from Bio import SeqIO

    t_start = time.time()
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Read FASTA ────────────────────────────────────────────────────
    print("\n[1/7] Reading input FASTA...")
    record     = SeqIO.read(input_fasta, "fasta")
    genome_seq = str(record.seq)
    accession  = record.id
    genome_len = len(genome_seq)
    prefix     = prefix or accession
    print(f"      {accession} — {genome_len:,} bp")

    # Load gene catalog
    with open(GENE_CAT) as f:
        gene_catalog = json.load(f)

    # ── Step 2: IR boundary detection ────────────────────────────────────────
    print("\n[2/7] Detecting IR/LSC/SSC boundaries...")
    t2 = time.time()
    from .identify.ir_detector import detect_ir_boundaries
    ir_boundaries = detect_ir_boundaries(genome_seq)
    if ir_boundaries is None:
        ir_boundaries = {"LSC": (0, genome_len)}   # IR-lacking / single-copy plastome
        print("      No IR pair >= 10 kb found — treating as single-copy (LSC only)")
    for region, (s, e) in ir_boundaries.items():
        print(f"      {region}: {s:,} - {e:,} ({e-s:,} bp)")
    print(f"      ({time.time()-t2:.1f}s)")

    # ── Step 3: Closest relatives ────────────────────────────────────────────
    print("\n[3/7] Finding closest relatives...")
    t3 = time.time()
    from .identify.closest_rel import find_closest_relatives
    rel_result = find_closest_relatives(
        query_fasta = input_fasta,
        blast_db    = BLAST_DB,
        gene_catalog= gene_catalog,
        top_n       = 10,
        threads     = threads,
    )
    relatives = rel_result["relatives"]
    print(f"      {len(relatives)} relatives ({time.time()-t3:.1f}s)")
    if relatives:
        top = relatives[0]
        print(f"      Top: {top['accession']} "
              f"({top['pident']:.1f}%) "
              f"{top.get('genus','')} "
              f"{top.get('family','')}")

    # ── Step 4A: Engine A — Reference-based CDS ───────────────────────────────
    print("\n[4A/7] Engine A: Reference-based annotation...")
    t4a = time.time()
    from .identify.engine_a import run_engine_a, load_reference_proteins
    ref_proteins = None
    if reference:
        ref_proteins = load_reference_proteins(reference)
        print(f"      Using --reference {reference}: "
              f"{len(ref_proteins)} genes overlaid (auto fallback to DB)")
    features_a = run_engine_a(
        genome_seq   = genome_seq,
        ir_boundaries= ir_boundaries,
        relatives    = relatives,
        protein_db   = PROTEIN_DB,
        gene_catalog = gene_catalog,
        threads      = threads,
        ref_proteins = ref_proteins,
    )
    cds_a  = sum(1 for f in features_a if f.gene_type=="CDS")
    rrna_a = sum(1 for f in features_a if f.gene_type=="rRNA")
    print(f"      CDS: {cds_a} rRNA: {rrna_a} ({time.time()-t4a:.1f}s)")

    # ── Step 4B: Engine B — Model-based annotation ────────────────────────────
    print("\n[4B/7] Engine B: Model-based annotation...")
    t4b = time.time()
    from .identify.engine_b import run_engine_b
    features_b = run_engine_b(
        genome_seq   = genome_seq,
        ir_boundaries= ir_boundaries,
        relatives    = relatives,
        hmm_dir      = str(HMM_DIR),
        trna_db_dir  = TRNA_DB_DIR,
        exon_db      = EXON_DB,
        rrna_dbs     = RRNA_DBS,
        gene_catalog = gene_catalog,
        threads      = threads,
        use_trnascan = use_trnascan,
    )
    cds_b  = sum(1 for f in features_b if f.gene_type=="CDS")
    trna_b = sum(1 for f in features_b if f.gene_type=="tRNA")
    rrna_b = sum(1 for f in features_b if f.gene_type=="rRNA")
    print(f"      CDS: {cds_b} tRNA: {trna_b} rRNA: {rrna_b} "
          f"({time.time()-t4b:.1f}s)")

    # ── Step 5: Reconciliation ────────────────────────────────────────────────
    print("\n[5/7] Reconciliation...")
    t5 = time.time()
    annotations = reconcile(
        engine_a    = features_a,
        engine_b    = features_b,
        genome_seq  = genome_seq,
        gene_catalog= gene_catalog,
    )

    high   = sum(1 for a in annotations if a.flag=="HIGH")
    medium = sum(1 for a in annotations if a.flag=="MEDIUM")
    review = sum(1 for a in annotations if a.flag=="NEEDS_REVIEW")
    print(f"      Total: {len(annotations)} genes "
          f"(HIGH={high} MEDIUM={medium} REVIEW={review}) "
          f"({time.time()-t5:.1f}s)")

    # ── Step 6: Special cases ─────────────────────────────────────────────────
    print("\n[6/7] Handling special cases...")
    from .annotate.special_cases import run_all_special_cases
    annotations, changes = run_all_special_cases(
        annotations  = annotations,
        genome_seq   = genome_seq,
        ir_boundaries= ir_boundaries,
        gene_catalog = gene_catalog,
        protein_db   = PROTEIN_DB,
        trna_db_dir  = TRNA_DB_DIR,
    )
    for change in changes:
        print(f"      {change}")

    # ── Step 6b: Splice-site refinement (multi-exon CDS) ──────────────────────
    from .annotate.refine_splice import refine_all
    n_refined = refine_all(annotations, genome_seq, gene_catalog)
    if n_refined:
        print(f"      Splice-site refined: {n_refined} multi-exon CDS")

    # ── Step 7: Output ────────────────────────────────────────────────────────
    print("\n[7/7] Writing output files...")
    from .output.writers import write_all
    files = write_all(
        annotations  = annotations,
        genome_seq   = genome_seq,
        accession    = accession,
        genome_len   = genome_len,
        ir_boundaries= ir_boundaries,
        relatives    = relatives,
        out_dir      = out_dir,
        prefix       = prefix,
        no_plot      = no_plot,
    )
    for f in files:
        print(f"      ✅ {f}")

    elapsed = time.time() - t_start
    n_genes = len(annotations)

    print(f"\n{'='*60}")
    print(f"  Done in {elapsed:.1f}s — {n_genes} genes annotated")
    print(f"  HIGH={high} MEDIUM={medium} NEEDS_REVIEW={review}")
    print(f"{'='*60}\n")

    return {
        "annotations" : annotations,
        "ir_boundaries": ir_boundaries,
        "relatives"   : relatives,
        "elapsed"     : elapsed,
    }
