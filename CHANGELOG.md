# Changelog

All notable changes to Plastanno v2. Benchmarks are measured on the DEV split
(n=123 shared genomes) against reference GenBank annotations; the held-out set is
never used during development.

## [Unreleased]

### Fixed

- **The 3.0.0 entry miscounted its own test suite.** It said "Thirteen test
  files, 206 checks, none needing an external dataset". The tag has fourteen
  files and 365 checks, and one of them, `test_relatives_ordering.py`, does need
  the downloaded database — it passed `protein_db="database/protein_db"` and
  raised `IndexError` rather than saying so when the directory was absent. The
  count is corrected above and the test now skips with a message.
- **Circular map: names containing `_` were drawn wrongly.** The organism name
  is italicised with mathtext, and the name was interpolated into
  `$\mathit{...}$` without escaping, so mathtext read its markup characters as
  markup. When `--organism` was not given the accession took its place, and
  `NC_053537.1` was drawn as "NC" with a subscript zero followed by "53537.1" —
  in the title and the centre label, in all three output formats. All of
  `$ { } ^ _ # & % ~` and the backslash are now escaped.
- **`scripts/get_database.sh` downloaded the previous database.** It kept its own
  copy of the Zenodo URL and MD5, which was not updated when the 3.0.0 database
  was published. The superseded record still exists and its checksum still
  matched, so the script fetched the older bundle — the one without the eight
  added genes — and exited 0. It now reads both constants from
  `plastanno/fetch_db.py`, the single place they are written.

### Changed

- **Gene names follow the majority usage in GenBank plastomes: psbN, ycf3, ycf4,
  clpP.** The synonym table mapped psbN to pbf1, a newer name that 88-97% of
  GenBank plastome records do not use (23% of those submitted in 2024 do). It now
  maps pbf1 to psbN, and pafI and pafII -- which the 3.0.0 database carries as
  profiles of their own beside ycf3 and ycf4 -- to ycf3 and ycf4; clpP1 to clpP is
  unchanged. The rename happens where it always did, after reconciliation, so no
  coordinate moves: on 30 development genomes every feature was identical once
  the old names were mapped. A renamed call that overlaps an existing call of the
  same gene is dropped rather than written twice. Renaming before reconciliation
  was tried and rejected: it changed which call won, not only its name.
- **Documentation brought up to 3.0.0.** `README.md` still described the
  pre-3.0.0 architecture and listed tRNAscan-SE as optional, which it is not
  under the default `--trna-mode hybrid`; `numpy`, `scipy` and `platformdirs`
  were missing from the dependency list. `docs/ARCHITECTURE.md` still gave the
  fixed-weight scoring formula that 3.0.0 does not use and reported 81 HMM
  profiles where the shipped database has 168.
- **Performance claims withdrawn from the public README.** The held-out F1 and
  the head-to-head against one other tool have been removed rather than
  updated: the reference databases were built before the evaluation split, so
  that set was largely represented in them, and the comparison used a scorer
  that has since been replaced. Measured performance will be published with the
  manuscript after an independent review of the scorer and protocol.
- **Schematic figures redrawn** from the code. `Fig2_reconciliation` is replaced
  by `Fig2_selection`, since the cross-engine reconciliation it depicted is
  neither the default nor a claim still made. `Fig3`, `Fig4` and `Fig7` are
  removed for the reasons above — `Fig7` asserted a leakage-free evaluation set
  that was not leakage-free.
- CI installs the dependencies 3.0.0 actually needs, and runs the test suite.

## [3.0.0] — 2026-09-26

The pooled-architecture line reaches the public release — 149 commits in the
development repository between the 2.0.5 and 3.0.0 tags. The previous release
predates that line entirely, so this is the architecture the tool was refactored
into rather than a patch on top of it.

**A major version because upgrading changes results and requires more of the
environment.** Every changed default keeps a `legacy` setting that reproduces the
2.0.5 behaviour, so a run can be pinned to the old output if needed.

### Changed — defaults that alter annotation output

- `--mode pooled` — candidate pooling with a B-primary/A-rescue rule, replacing
  the cross-engine integration layer. `--mode legacy` restores it.
- `--trna-mode hybrid` — intron-free tRNA boundaries come from tRNAscan-SE while
  the set of loci, intron models and naming stay as before. Costs roughly 4x
  runtime. `--trna-mode legacy` restores ARAGORN/exon-DB boundaries and removes
  the tRNAscan-SE requirement.
- `--exon-mode aragorn` — ARAGORN's own structural intron call is preferred over
  the exon-database call at intron-bearing tRNA loci. `--exon-mode legacy`
  restores the previous ranking.
- `--intron-mode glocal` — **new.** For the seven intron-bearing tRNA (trnA-UGC,
  trnI-GAU, trnL-UAA, trnK-UUU, trnV-UAC, trnG-UCC, trnG-GCC), each reference
  exon is aligned end to end inside the existing call's window and donors vote
  per coordinate, so the boundary follows reference geometry instead of where
  local similarity tails off. It does not model the intron; the intron is the gap
  between the two exons. Inventory, naming and every other feature class are
  untouched. `--intron-mode legacy` restores ARAGORN's structural call.

### Changed — requirements

- **tRNAscan-SE is now required**, not optional. The default tRNA mode needs it
  and `pipeline.run` checks for the binary before any work starts, raising
  `TrnascanUnavailable` rather than failing later. Use `--trna-mode legacy` to
  run without it.
- **numpy and scipy are now required.** Reconciliation pairs the two engines'
  calls with a maximum-weight assignment (`scipy.optimize.linear_sum_assignment`)
  over a numpy cost matrix, on essentially every genome.
- `environment.yml` had described tRNAscan-SE as optional and omitted numpy and
  scipy; it, `pyproject.toml` and the conda recipe now agree.

### Fixed

- **`tRNA-???` aborted the whole genome.** ARAGORN writes that when it finds a
  tRNA structure it cannot assign an amino acid to — ordinary output, not
  corruption. The fail-loud parser treated it as an unreadable coordinate and
  raised, losing the entire annotation; three of sixty genomes in one benchmark
  sample died this way. One unassignable line now costs one candidate, with a
  warning. Genuinely malformed coordinates still raise.
- **`plastanno run --help` crashed.** A help string contained `84.1% legacy`,
  which argparse read as a format specifier (`% l` + `e`), so both `run --help`
  and `batch --help` raised instead of printing.
- **`fetch-db` could install where the tool would not look.** `fetch_db` and
  `paths` each resolved the user-data directory separately and the two copies had
  drifted: without `platformdirs`, `fetch-db` wrote to `~/.local/share/plastanno`
  while `db_root` skipped that location and returned the repo layout — so a user
  could download 266 MB and still be told the database was missing.
  `paths.user_data_parent` is now the single source of truth for both.

### Changed — reference database

The bundle gains eight genes the previous one lacked: chlB, chlL, chlN, lhbA,
pafI, pafII, rpl21 and ycf66, each with an HMM profile and a protein FASTA, with
`hmm_db/all_profiles.hmm` regenerated. Everything else is byte-identical to the
previous bundle. Published as Zenodo record 22960386 (concept DOI
10.5281/zenodo.20807994 continues to resolve to the newest version).

**Existing installations must re-fetch:** `plastanno fetch-db --force`.

`database_CHECKSUMS.sha256` is regenerated with 5634 entries. It now also covers
`boundary_db/*` and `exon_templates.json`, which ship in the repository but were
missing from the previous list — so `sha256sum -c` used to fail for every user.
The six `*.bak` entries, development snapshots that are not part of the package,
are dropped.

### Added — tests

Fourteen test files, 365 checks. All but one need no external data;
`test_relatives_ordering.py` drives Engine A against `database/protein_db`
and skips with a message when the database has not been downloaded:

    for t in tests/test_*.py; do python3 "$t" || echo "FAILED: $t"; done

The benchmark scoring machinery and its own tests are not in this repository;
they score against reference files that are not distributed. See
`tests/README.md`.

## [2.0.5] — 2026-09-11

Reported by a user: ycf1 was missing entirely from the output, and rps12 was
annotated as a 69 kb coding sequence. Both turned out to be symptoms of wider
defects.

### Fixed — annotation

- **ycf1 was deleted from 23 of the 111 DEV genomes that carry it.** The minimum
  length filter in `reconcile._select` drops any CDS below 0.6x its expected
  length; ycf1 is the most divergent plastid ORF, so Exonerate usually reports only
  a short conserved core (807 bp of a 5040 bp gene) and the IR-junction copy is
  genuinely truncated. ycf1 is now exempt from that filter, the length decision
  moves to `handle_ycf1` after ORF completion, and the completion window scales
  with the gene's own coding length. When only one copy is found the other is
  derived from the measured inverted-repeat symmetry (median error 3 bp) and marked
  as derived. *ycf1 F1 70.9 -> 79.8; ycf1 pseudogene sensitivity 0 -> 35/35.*
- **Frame-broken CDS.** Exonerate models frameshifts and reports them in its exon
  attributes, which the parser discarded; the resulting spans are not whole codons
  and both existing repair steps refuse to touch them. The indel count is now kept,
  and such a boundary is re-derived by scanning all three frames for a complete ORF.
  *46 -> 20 broken CDS (0.45% -> 0.19%); 19 of 26 rescues match the reference exactly.*
- **Missing 5' exon.** Short-exon recovery covered only petB/petD/rpl16 and searched
  first-exon lengths of 3-12 bp, so it could not recover an exon that is merely
  absent: rpoC1 (435 bp), atpF (145 bp), rps16 (40 bp) were annotated as exon 2
  alone. Recovery now fires whenever the catalog says the gene is multi-exon but one
  exon was produced, sizing its search from the reference alignment. *29 of 40 cases
  fixed.*

*DEV: F1 92.41 -> 92.73, CDS F1 93.02 -> 93.51, CDS TP 9574 -> 9631, FP 633 -> 590.
30 genomes improved, 0 regressed (sign test p < 1e-5). tRNA and rRNA unchanged.*

### Added — output

- **`.tbl` (NCBI 5-column feature table).** NCBI does not accept a GenBank flatfile
  as a submission; table2asn needs a feature table. Without one the annotation could
  not be submitted at all. Validated with NCBI's table2asn over six rounds:
  20 tool-attributable errors -> 3, and the three that remain are annotation
  boundaries, not format.

### Fixed — output format

- **GFF3 reported every intron-containing gene as one uninterrupted CDS**, ignoring
  exons. rps12 came out as ~70 kb; 28 genes were affected (trnK-UUU 2599 -> 73 bp,
  trnI-GAU 1015 -> 104 bp, clpP 2048 -> 591 bp). CDS lines are now one per exon with
  the spec-required phase.
- **Duplicate GFF3 feature IDs** (106 genomes): both rps12 copies share their 5'
  exon and so shared an ID, which a parser merges into one six-exon gene.
- **Gene features were written as `join()`** for all 23 intron-containing genes. A
  gene is one interval spanning its introns; only the coding feature is spliced.
  Trans-spliced rps12 keeps its join, with `/trans_splicing`.
- Pseudogenes follow the NCBI convention (a gene with `/pseudo`, no CDS) and stay
  out of the .faa/.ffn files; a CDS whose length is not a whole number of codons is
  written as partial with no invented translation; an RNA-edited initiator carries
  `/transl_except`.

### Fixed — benchmark

- **tRNA scoring under-counted by 128 features per run.** A tRNA may be named after
  its anticodon (trnA-UGC) or the codon it reads (trnA-GCA); the two are reverse
  complements of the same gene, and reference records use both. Requiring identical
  triplets scored correct calls as a false negative plus a false positive.
  *Measured tRNA F1 90.70 -> 93.53 — a correction to the measurement, not an
  improvement to the tool.* Numbers scored before and after this change are not
  directly comparable.

### Known limitations

- 20 CDS (0.19%) still have no readable frame, against 0.01% for PGA and 0.04% in
  the NCBI references. PGA reaches its figure by not annotating genes it cannot
  resolve (5.4 genes per genome in its own logs); Plastanno keeps them, flagged.
- 71.7% of remaining boundary errors are at the 5' end, median 63 bp.
- ycf15 and ycf68 are in the gene catalog but have no protein database, so neither
  engine can find them (92 false negatives on DEV).
- `NC_084420.1` exceeds the 300 s self-BLASTN timeout in IR detection (67,715 HSPs);
  pre-existing.
