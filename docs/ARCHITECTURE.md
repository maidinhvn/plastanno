# Plastanno v2 Architecture

## Pipeline Overview
## Scoring Formula (Reconciliation)

C = 0.4·S_overlap + 0.2·S_ref + 0.2·S_model + 0.2·S_orf

S_overlap = Jaccard overlap A∩B / A∪B
S_ref     = pident/100 (Exonerate identity)
S_model   = normalized HMM bitscore
S_orf     = ORF validity (start/stop/length)

Confidence flags:
  C ≥ 0.8 → HIGH
  C ≥ 0.5 → MEDIUM
  C < 0.5 → NEEDS_REVIEW

## Key improvements over v1

1. Reconciliation layer (NEW)
   - v1 had NO reconciliation
   - v2: explicit scoring + provenance

2. Engine B model-based (IMPROVED)
   - v1: HMM only for naming
   - v2: HMM for detection too

3. Pre-built protein DB (NEW)
   - v1: on-the-fly from relatives (slow)
   - v2: pre-built per gene (fast)

4. RNA editing table (NEW)
   - v1: missing → false pseudogenes
   - v2: bảng tra start codon + editing sites

5. Provenance per feature (NEW)
   - v1: no tracking
   - v2: engine source + score components

## Databases

database/
├── blast_db/       ← genus reps (2,899 genomes)
├── protein_db/     ← per-gene AA sequences (Exonerate)
├── trna_db/        ← hierarchical (genus/family/global)
├── exon_db/        ← actual exon seqs (intron tRNAs)
├── rrna_db/        ← full-length rRNA (1400-1600bp)
├── hmm_db/         ← 81 HMM profiles
└── gene_catalog.json ← metadata + special cases
