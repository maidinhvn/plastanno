# Plastanno architecture

Current as of 3.1.0. Where this document and the code disagree, the code is
right — the previous version of this file described the pre-3.0.0 scoring and
was three months stale.

## Pipeline

`plastanno/pipeline.py::run()` orchestrates the steps; the `Step N` comments in
that file are the authoritative numbering.

| Step | Module | What it does |
|---|---|---|
| 1 | `pipeline` | read a single plastome record from FASTA |
| 2 | `identify/ir_detector` | self-BLASTN; longest minus-strand HSP ≥ 10 kb → `{LSC, IRb, SSC, IRa}`. Returns `None` for IR-lacking plastomes, and the whole genome is then one `LSC` region |
| 3 | `identify/closest_rel` | BLAST against `blast_db/genus_reps`; ranked neighbours for the tRNA search |
| 4A | `identify/engine_a` | Exonerate `protein2genome` per gene against `protein_db/`, both IR copies for IR genes; a gene whose catalog region gives no hit reaching 0.6× its expected length is also searched in IRb and IRa, and those hits carry an `[IR search]` note; rRNA via BLAST → `s_ref` |
| 4B | `identify/engine_b` | 6-frame translation → `hmmsearch` (CDS); ARAGORN + tRNAscan-SE + BLAST (tRNA); BLAST (rRNA) → `s_model` |
| 5 | `core/reconcile` | candidate selection — see below |
| 6 | `annotate/special_cases` | CAU disambiguation, *rps12* trans-splicing, short first exons, internal-stop QC |
| 6b | `annotate/refine_splice`, `identify/trna_hybrid`, `identify/intron_refine` | boundary refinement. `refine_all` runs, in order: exon-panel junctions of multi-exon CDS (a gene marked `keep_ends`, rpl2, keeps its ends); the terminal stop (trim a short read-through, complete a truncated 3' end; `[3' completed]`); the start codon (move a start that is not a start codon to the best-supported one nearby; `[start moved]`). Then intron-free tRNA ends and intron-bearing tRNA exons |
| 6c | `core/finalize`, `core/reconcile` | `finalize_qc`, `revoke_implausible_rescues`, CDS products from the catalog, then `submission_check`: a CDS that NCBI's validator would reject, or that carries one of the three repair notes, is NEEDS_REVIEW (`[submission]`) |
| 7 | `output/writers` | `.gb .gff3 .tbl .faa .ffn .frn .report .provenance.json .trna_alternatives.tsv` |

Step 6b runs **after** everything that decides which loci exist. That ordering is
load-bearing, not stylistic: `_select` clusters same-name features that overlap
at all, so moving a boundary earlier in the pipeline turns a boundary change into
an inventory change — two IR duplicates nudged apart stopped overlapping and both
survived.

## Candidate selection (`--mode pooled`, the default)

Both engines propose candidates; they are pooled by gene name and one is kept per
locus.

- **both engines found it** — keep Engine B's coordinates when its ORF check
  passes, or is at least as good as Engine A's; substitute A only otherwise.
- **one engine found it** — keep it. This rescue is what lifts inventory.
- **same name, too far apart to be one locus** — for an IR-duplicated gene those
  are the two copies and both are kept; otherwise ORF validity decides.

**Agreement between the engines earns no score bonus.** The H3 ablation found
that the cross-engine integration in the older `reconcile()` — AB matching, the
coordinate-donor rule, the agreement boost, conflict resolution — bought no
detectable improvement in strict CDS coordinate concordance over pooling the same
candidates through the same selector. What did help was Engine B's coordinates:
on the development set Engine B alone had the best strict CDS precision. So the
candidates were kept and the integration dropped. `--mode legacy` restores it.

### The score

Every feature is scored on single-engine signals only: `{ref, orf}` for a CDS
Engine A found, `{model, orf}` otherwise, as a weighted mean over the signals
actually available, renormalised to 1. There is no overlap term — in pooled mode
`s_overlap` is zero by construction.

    flag:  ≥ 0.8 HIGH   ·   ≥ 0.5 MEDIUM   ·   otherwise NEEDS_REVIEW

This is a **candidate-selection score**: it ranks candidates for one locus and
nothing more. It is deliberately not called a confidence. An earlier version of
this document gave a fixed-weight formula
`C = 0.4·S_overlap + 0.2·S_ref + 0.2·S_model + 0.2·S_orf`; that formula is gone,
and it capped any single-engine feature at 0.4 and forced NEEDS_REVIEW.

`_select` then union-finds features into per-locus clusters (same name + any
overlap → duplicate fragments; different names + reciprocal overlap > 0.5 and
length ratio > 0.5 → paralog cross-hit), keeps the best per cluster, and drops
CDS whose spliced length is below 0.6× expected. `rps12` passes through untouched
and is handled in step 6.

**A dropped locus is rescued (3.0.1).** When both engines found a gene, the call
not kept is held aside. If the kept call does not survive `_select` -- typically a
fragment removed by the length filter, while the other engine had the complete
gene -- the held-aside call is tried: `_select` is re-run on copies of the
survivors plus the candidate, and the candidate is added only if every survivor
comes back unchanged, so nothing that would have been output can move or vanish.
After the final QC, a rescued CDS that is flagged NEEDS_REVIEW or contains an
in-frame stop is withdrawn (`reconcile.revoke_implausible_rescues`). ndhA is the
gene this most often recovers.

## Provenance

Every annotation is a `Feature` (`core/feature.py`) carrying its source engine,
the component scores, the flag, exons, pseudogene status and free-form notes. The
writers propagate all of it into GenBank `/note` and GFF3 attributes, and into
`.provenance.json`.

## Databases

Counts are from the published 3.0.0 bundle.

    database/
    ├── blast_db/           genus representatives (BLAST nucleotide DB)
    ├── protein_db/         per-gene amino-acid sequences for Exonerate (403 files)
    ├── hmm_db/             profile HMMs for 88 genes
    ├── trna_db/            tiered tRNA (genus / family / global)
    ├── exon_db/            exon sequences for intron-bearing tRNA
    ├── rrna_db/            full-length rRNA
    ├── boundary_db/        exon panel and length templates
    ├── exon_templates.json
    └── gene_catalog.json   per-gene region, n_exons, expected_len

`all_profiles.hmm` holds **168 records for those 88 genes**: 80 of them appear
twice, byte-identical, and the 8 added in 3.0.0 appear once. Quote the gene
count, not the record count — an earlier version of this file said "168
profiles", which reads as 168 genes. hmmsearch therefore scans 80 profiles twice
on every genome; deduplicating the file is a free saving nobody has taken.

Two caveats that have each cost debugging time:

- **`gene_catalog.json` exists twice.** `paths.config_dir()` prefers the packaged
  `plastanno/data/gene_catalog.json` and falls back to `database/gene_catalog.json`
  only if that is absent. Edits to the latter are inert while the former exists.
- **The tRNA tiers do not engage.** `closest_rel` reads taxonomy from
  `database/splits/train.csv`, which is not shipped, so the genus and family tiers
  are skipped and everything uses `global`. Restoring the file was measured and
  rejected: the tiers did engage on 16 of 30 genomes, but tRNA identity was
  unchanged and the acceptor-stem score moved −0.008. Its absence is a decision.
