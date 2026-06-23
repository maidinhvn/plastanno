# Database provenance and reproducibility

This note records how Plastanno's reference databases are assembled and how the
held-out evaluation is kept free of train–test leakage, so the results can be
reproduced and audited.

## Reference collection and split

Reference data is drawn from public RefSeq plastomes. The collection is
partitioned with a fixed seed, stratified by structural mode, into a development
set and a frozen held-out set (see `splits/split_manifest.json`, which records
the seed and the sha256 checksums of both sets). Tuning and calibration use the
development set.

## Databases

| Database | Content |
|---|---|
| `blast_db/genus_reps` | genus-representative genomes (2,899) for closest-relative search |
| `protein_db/` | per-gene CDS proteins (Engine A / Exonerate) |
| `hmm_db/` | profile HMMs (Engine B) |
| `trna_db/` | hierarchical tRNA (genus / family / global) |
| `exon_db/` | exon sequences for intron-containing tRNAs |
| `rrna_db/` | full-length rRNA per gene |
| `boundary_db/`, `exon_templates.json` | per-gene exon panel + length templates (splice refinement) |
| `gene_catalog.json` | curated per-gene metadata (region, exon count, expected length) |

Only the small configuration files (`gene_catalog.json`, `exon_templates.json`,
`boundary_db/`) are version-controlled. The large sequence and HMM databases are
distributed via Zenodo (DOI in the README) or can be rebuilt with
`scripts/build/build_all.py`.

## Leakage-free held-out evaluation

Generalisation is measured on a **leakage-free test set of 2,151 land-plant
plastomes whose own sequences are verified to be absent from every reference
database** (the per-gene proteins, the profile HMMs, the tRNA and tRNA-exon
databases, the rRNA databases, and the genus-representative set). It comprises:

- **498** plastomes held out from development and confirmed absent from all
  databases, and
- **1,653** RefSeq plastomes deposited after the databases were built, hence
  never present during database construction or tuning.

Because no test genome's sequence appears in any reference database, the reported
global F1 of 92.5% reflects performance on plastomes the tool has genuinely never
seen, rather than memorised reference data.
