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

## The held-out evaluation, and why its claim was withdrawn

An earlier version of this document stated that generalisation was measured on a
**leakage-free test set of 2,151 land-plant plastomes whose sequences are absent
from every reference database**, and reported a global F1 from it.

That premise does not hold. The reference databases were built **before** the
evaluation split was drawn, so most of that set was in fact represented in them.
The set was held out from *tuning*, which is a real and useful property, but it
was not held out from *database construction* — and it is the databases that
Engine A and Engine B search. The F1 measured on it therefore does not separate
generalisation from retrieval of material already in the databases, which is
exactly what the claim asserted it did.

Both the claim and the number are withdrawn rather than adjusted. What the set
still is:

- **498** plastomes held out from development, and
- **1,653** RefSeq plastomes deposited after the databases were built.

The second stratum is genuinely unseen by construction — a plastome deposited
later cannot be in an earlier database. The first is not, in general.

`Plastanno_dataset_inventory.xlsx` lists both the development collection and this
set, accession by accession, so the composition is inspectable regardless of what
is claimed about it. A corrected evaluation, on a pool built after the split
rather than before it, will be published with the manuscript.
