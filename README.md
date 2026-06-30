# Plastanno v2

**A hybrid, self-evaluating chloroplast-genome (plastome) annotator.**

[![CI](https://github.com/maidinhvn/plastanno/actions/workflows/ci.yml/badge.svg)](https://github.com/maidinhvn/plastanno/actions/workflows/ci.yml)
[![Reference databases (Zenodo)](https://zenodo.org/badge/DOI/10.5281/zenodo.20807994.svg)](https://doi.org/10.5281/zenodo.20807994)

Given a plastome FASTA, Plastanno predicts CDS, tRNA and rRNA features and writes
GenBank, GFF3, FASTA and report outputs — plus a circular plastome map. Its
defining design is a **hybrid of two independent annotation engines whose results
are merged by a scoring/reconciliation layer** that assigns every feature a
normalised confidence score, a review flag, and full provenance.

On a **leakage-free held-out set of 2,151 land-plant plastomes** (sequences not
represented in any reference database), Plastanno reaches a gene-by-gene
**F1 of 92.5%** (CDS 93.0, tRNA 90.3, rRNA 97.2). In a head-to-head on a common
subset it scores **92.0 vs 87.8 F1** against PGA, winning on 145 of 245 genomes.

---

## How it works

A 7-step pipeline (`plastanno/pipeline.py`):

1. **Read FASTA** — single plastome record.
2. **IR detection** (`identify/ir_detector.py`) — self-BLASTN; the longest
   minus-strand HSP ≥ 10 kb defines the inverted-repeat pair → `{LSC, IRb, SSC, IRa}`.
3. **Closest relatives** (`identify/closest_rel.py`) — BLAST against genus
   representatives; taxonomy selects the tRNA-DB tier.
4. **Two engines, in parallel:**
   - **Engine A — reference-based** (`identify/engine_a.py`): Exonerate
     `protein2genome` per gene (IR genes searched in both copies); rRNA via BLAST.
   - **Engine B — model-based** (`identify/engine_b.py`): CDS via 6-frame
     translation → `hmmsearch`; tRNA via ARAGORN + BLAST; rRNA via BLAST.
5. **Reconciliation** (`core/reconcile.py`) — bin (overlap) → adaptive confidence
   → ORF validation → locus selection; one annotation set with flags + provenance.
6. **Special cases** (`annotate/special_cases.py`) — CAU tRNA disambiguation,
   *rps12* trans-splicing, short first exons, multi-exon splice refinement,
   internal-stop QC.
7. **Output** (`output/writers.py`).

See `docs/figures/Fig1_pipeline.*` and `Fig2_reconciliation.*` for schematics.

## Installation

**Supported platforms:** Linux and macOS. Windows is not supported directly (the
external tools and the helper shell scripts assume a Unix environment) — use
**WSL2** as a workaround.

Plastanno needs Python ≥ 3.9 with `biopython`, `pandas` and `matplotlib` (the
last only for the circular map), plus four external tools on `PATH`: **BLAST+**,
**Exonerate**, **HMMER** (`hmmsearch`) and **ARAGORN** (optionally
**tRNAscan-SE** for the `--trnascan` option).

**Prerequisite:** a working `conda`. We recommend
[Miniforge](https://github.com/conda-forge/miniforge) (it defaults to the
conda-forge channel and avoids the Anaconda Terms-of-Service prompt noted below),
or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

### Recommended: install from Bioconda

The easiest way — conda pulls in `plastanno` together with every Python
dependency and external tool. No cloning, no `pip`:

```bash
conda config --set channel_priority strict   # recommended for bioconda
conda create -n plastanno -c conda-forge -c bioconda plastanno
conda activate plastanno

# REQUIRED before the first run: download the reference database (~266 MB)
plastanno fetch-db
```

To add the optional tRNAscan-SE source (used only by the `--trnascan` option):

```bash
conda install -n plastanno -c conda-forge -c bioconda trnascan-se
```

> List `conda-forge` before `bioconda` so it keeps the higher priority — this is
> the channel order Bioconda requires; with `channel_priority strict` it also
> speeds up the solver and avoids mixing incompatible builds.
>
> If `conda create` stops with a `CondaToSNonInteractiveError` (Anaconda Terms of
> Service not accepted), it is coming from the `defaults` channel — use Miniforge,
> which has no `defaults` channel, or append `--override-channels` to this one
> command to exclude it.

### From source (developers / unreleased code)

Use this only to run a development checkout. Here `pip install .` is run **from
inside the cloned repository**, so the leading `git clone && cd` matters:

```bash
git clone https://github.com/maidinhvn/plastanno.git
cd plastanno

# create the environment (Python deps + all external tools) from the pinned file
conda env create -f environment.yml
conda activate plastanno

# install the `plastanno` command from this checkout
pip install .

# REQUIRED before the first run: download the reference database (~266 MB)
plastanno fetch-db
```

`environment.yml` pins the channels and lists every dependency. To build the
environment by hand instead, the equivalent one-liner is:

```bash
conda create -n plastanno -c conda-forge -c bioconda \
    python=3.10 biopython pandas matplotlib platformdirs blast exonerate hmmer aragorn
```

(Same channel-order / Terms-of-Service notes as above apply.) You can also skip
`pip install .` and run the tool directly from the checkout with
`python3 plastanno.py …` (see [Usage](#usage)).

`plastanno fetch-db` downloads the database into a platform data directory
(`~/.local/share/plastanno/database` on Linux). Override the location with
`$PLASTANNO_DB`. Running from a source checkout without installing also works:
use `python3 plastanno.py …` and `bash scripts/get_database.sh` (which places the
database under `./database`).

Then verify everything is found (modules, databases, external tools):

```bash
bash diagnose.sh
```

The Exonerate executable is resolved from `$PLASTANNO_EXONERATE`, then `PATH`,
then a built-in default — set `PLASTANNO_EXONERATE` if yours is elsewhere.

## Usage

> **Before running:** the reference database must be present — run
> `plastanno fetch-db` once (see [Databases](#databases)). Without it, annotation
> stops at the "Finding closest relatives" step with a BLAST database error.

```bash
# Annotate one genome -> 6 files + a circular map (.png/.pdf/.svg) in out/
plastanno run genome.fasta --output out/ --threads 8

# Skip the map (e.g. for large batches / benchmarks)
plastanno run genome.fasta --output out/ --no-plot

# Annotate every *.fasta / *.fa in a directory
plastanno batch genomes/ --output out/ --threads 8
```

> Running from a source checkout without `pip install`? Use `python3 plastanno.py
> run …` (and `bash scripts/get_database.sh`) — the same commands work unchanged.

### Optional inputs

Both are off by default, so default output is unchanged:

```bash
# Overlay a closely-related annotated reference (GenBank): its per-gene proteins
# are tried first in Engine A, with automatic fallback to the built-in DB.
python3 plastanno.py run genome.fasta -o out/ --reference close_relative.gb

# Add tRNAscan-SE as an extra tRNA source (must be on PATH); contributes
# intronless tRNAs alongside ARAGORN/BLAST.
python3 plastanno.py run genome.fasta -o out/ --trnascan
```

### Worked example

Two land-plant plastomes are bundled under [`example/`](example/) so you can try
the tool immediately (you still need the external tools and the full `database/`):

- `example/NC_053537.1.fasta` — *Gynostemma yixingense* (Cucurbitales)
- `example/NC_008325.1.fasta` — *Daucus carota* (carrot, Apiaceae)

```bash
# Annotate one of the examples (writes 6 files + a circular map)
python3 plastanno.py run example/NC_053537.1.fasta --output out_example/NC_053537.1 --threads 4

# ...or run the bundled script, which does a single run and a batch run
bash example/run_example.sh
```

`example/NC_053537.1.fasta` should yield ~129 genes; open
`out_example/NC_053537.1/NC_053537.1.report` for the categorised gene-group
summary, and `..._map.png` for the circular map.

### Outputs

| File | Contents |
|------|----------|
| `<acc>.gb` | GenBank with provenance in each `/note` |
| `<acc>.gff3` | GFF3 with confidence flags |
| `<acc>.faa` / `.ffn` / `.frn` | protein / CDS-nucleotide / RNA FASTA |
| `<acc>.report` | QC summary: IR boundaries, confidence distribution, a categorised functional gene table, and features flagged for review |
| `<acc>_map.png/.pdf/.svg` | circular plastome map (unless `--no-plot`) |

## Databases

`database/` holds the reference data the engines use: genus-representative BLAST
databases, per-gene reference proteins, profile HMMs, a taxonomically tiered tRNA
database, an intron-exon database, full-length rRNA databases, a per-gene exon
panel and length templates, and `gene_catalog.json`.

**The large reference databases (~266 MB compressed download, ~700 MB once
extracted) are not in this Git repository** —
GitHub's per-file size limit makes them unsuitable for version control. Only the
small runtime configs (`gene_catalog.json`, `boundary_db/`, `exon_templates.json`)
are tracked. To run the tool after cloning, obtain the full `database/` in one of
three ways:

1. **Automatic (recommended):** run
   ```bash
   bash scripts/get_database.sh
   ```
   which downloads the archive from Zenodo, verifies its checksum, and unpacks it
   into `database/`.
2. **Manual:** download `plastanno-database.tar.gz` from Zenodo
   ([doi.org/10.5281/zenodo.20807994](https://doi.org/10.5281/zenodo.20807994))
   and `tar xzf plastanno-database.tar.gz` so its contents sit in `database/`.
3. **Rebuild:** `python3 scripts/build/build_all.py` (see that script's header
   for inputs).

For a ready-to-run copy that already bundles the databases, use the release
tarball instead of cloning.

## Benchmarking

Quality is measured by gene-by-gene comparison against reference GenBank files
(true positive = name match + both ends within ±tol bp + sequence similarity):

```bash
# Score one predicted .gb against a reference .gb
python3 scripts/benchmark/benchmark_gene_by_gene.py reference.gb predicted.gb --tol 60 --sim 0.6

# Aggregate F1 over a sample
python3 scripts/benchmark/multi_genome_bench.py --n 120 --workers 16
```

## Performance targets

CDS Sn > 92% / Pr > 95%, rRNA Sn/Pr > 97%, tRNA Sn > 88% / Pr > 90%,
runtime < 60 s per genome.

## Citation

If you use Plastanno, please cite the manuscript (in preparation; see `docs/`) and
the reference-database archive on Zenodo:
[doi.org/10.5281/zenodo.20807994](https://doi.org/10.5281/zenodo.20807994).

## License

See [LICENSE](LICENSE).
