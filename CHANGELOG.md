# Changelog

All notable changes to Plastanno. Benchmarks are measured on the DEV split
(n=123 shared genomes) against reference GenBank annotations; the held-out set is
never used during development.

## [3.1.0] — 2026-10-01

A minor release: the same commands, options and reference database (no `fetch-db`
needed after upgrading), with new behaviour that changes annotations. A CDS that
would fail NCBI validation is now flagged, with the reason, so it can be fixed before submission. A truncated 3' end is completed and a
short read-through trimmed. A start that is not a start codon is moved to the
best-supported one nearby. A gene whose catalog region gives no usable hit is also
searched in the inverted repeats. The exon panel now covers rpl2.

**Every repair is marked in the CDS's note and keeps the CDS NEEDS_REVIEW.** More CDS
are flagged for review than in 3.0.1, by design: a repair makes an ORF valid, but it
does not prove the gene model right.

How it was checked:
- Each change was validated on development genomes against a criterion fixed before
  the run.
- The release as a whole was checked against 3.0.1 on 200 development genomes not used
  before: angiosperms, non-angiosperms, orchids and ndh-less plastomes. It passed every
  locked criterion but one, which was accepted with its reason (see Known
  limitations).
- Fresh environments were built from `environment.yml` and from the bioconda recipe
  (Python 3.10 and 3.14). On the genomes checked, they gave the same annotation as the
  development environment.
- Measured accuracy will be published with the manuscript.

### Added

- **A submission check.** Each CDS, as the writers will write it, goes through
  table2asn's three CDS checks for genetic code 11: a valid start codon, a stop
  codon, no internal stop codon.
  - The RNA-editing exception and partial CDS are exempt, exactly as they are for the
    validator.
  - A failing CDS becomes NEEDS_REVIEW with a `[submission]` note naming the problem
    and the flag it had.
  - The `.report` lists every such CDS in a new "Submission check" section.
  - Flags and notes only: no coordinate, feature or confidence changes.
- **3' ends are completed and read-throughs trimmed.** The terminal-stop step used to
  look only 3 codons ahead, within a length guard. It now:
  - trims a CDS back to its first in-frame stop when that stop lies within its last 5
    codons;
  - takes a stop within 3 codons silently, within the old length guard, as 3.0.1 did;
  - otherwise completes a CDS with no internal stop to the first in-frame stop within
    100 codons, reading round the origin.

  Trims and completions carry a `[3' completed]` note.
- **A start that is not a start codon is moved.**
  - When a CDS starts on a codon that is neither a table-11 initiation codon nor one of
    the gene's special start codons, a last pass looks up to 10 codons either way for
    ATG or one of those codons.
  - It never crosses an in-frame stop and stays inside the first exon.
  - It takes the start the reference proteins' N-termini support best.
  - Starts that are valid but differ from a reference are not touched.
  - A moved CDS carries a `[start moved]` note.
- **A gene is searched in the inverted repeats when its catalog region gives no usable
  hit.**
  - Where the IR has expanded over most of the SSC, the SSC genes sit in the IR
    copies. ndhA, the SSC gene with an intron, was then lost whole.
  - When no hit in the catalog region reaches 0.6× the gene's expected length, Engine A
    also searches IRb and IRa, and the IR hits replace the region hits they overlap.
  - Such a CDS carries an `[IR search]` note.
- **rpl2 is in the exon panel.** Exonerate starts rpl2's exon 2 one codon late against
  almost every reference. The panel gains four rpl2 references, and the rpl2 exon
  template takes their exon lengths. rpl2 is marked `keep_ends`: the panel moves its
  junction and leaves its 5' and 3' ends to the start-codon and terminal-stop passes,
  because many lineages start rpl2 on an edited ACG.
- **A warning when the exon panel cannot be searched.** A missing or unreadable panel
  database used to skip junction refinement for every panel gene in silence. The first
  failure now prints one WARNING naming the cause.
- `submission_check` keeps every CDS with a `[3' completed]`, `[start moved]` or
  `[IR search]` note NEEDS_REVIEW, whatever its ORF. The `.report` lists these CDS
  apart from those that would fail validation.
- **Eight new test files, 181 checks:**
  - `test_submission_check.py`;
  - `test_terminal_stop.py`;
  - `test_start_rescue.py`;
  - `test_ir_fallback_search.py`;
  - `test_panel_keep_ends.py`;
  - `test_panel_warning.py`;
  - `test_renamed_gene_profile.py`;
  - `test_synonym_lookups.py`, which checks the reference-protein files only when the
    downloaded database is present.

### Changed

- **The tool is called Plastanno.** "v2" is dropped from its name, which it had
  outlived: the version is 3.x. The README title, `--help`, the database builder, the
  report header and the DEFINITION line written when no organism is given now say
  "Plastanno".
- **More CDS are flagged NEEDS_REVIEW**, by design: the repairs above and the
  submission check each keep a CDS under review.

### Fixed

- **psbN's reference proteins were not found after the rename.**
  - Step 6 renames pbf1 to psbN before anything else, but the protein database keeps
    psbN's proteins in `pbf1.fasta`.
  - So from 3.0.1 on, no psbN was ORF-completed, and a psbN call a few codons short of
    its stop stayed short.
  - The lookup now falls back to the file of a name that `SYNONYMS` maps to the gene.
  - A test walks every renamed gene's lookups: proteins, catalog, special start codons,
    exon panel and templates.
- **The README did not say which F1 the benchmark scripts report.** With the default
  ±60 bp, a call one base off counts as exact. The aggregate script also reports the
  F1 at exact coordinates. Neither is a published accuracy figure.

### Known limitations

- A pseudogene can be written as a CDS instead of a pseudogene. The CDS carries
  NEEDS_REVIEW and a `[submission]` note that names the defect (frameshifts, internal
  stop codons, no stop codon). The search in the inverted repeat can add such a CDS
  for a gene whose catalog region gave no usable hit.
- In some genomes outside the flowering plants, chlB or chlN is written twice at the
  same coordinates (also in 3.0.1). NCBI's validator reports it as a duplicate
  feature, a warning.
- In ferns, the start-codon repair can move an ACG start, which RNA editing turns
  into a real start, to a downstream ATG. The CDS keeps a `[start moved]` note and
  NEEDS_REVIEW.
- In Cycadaceae, a spurious psbM fragment can end up on the same coordinates as psbZ.
  It is flagged NEEDS_REVIEW.

## [3.0.1] — 2026-09-29

A patch release: the same commands and options, and the same reference database
(no `fetch-db` needed after upgrading), but annotations change where 3.0.0 was
wrong or incomplete. Gene names follow the majority usage in GenBank plastomes;
every CDS carries its protein name as `/product`; and a gene the pooled rule used
to drop -- most often ndhA -- is rescued. Each change was validated on development
genomes against a criterion fixed before the run; the locus rescue on 60 genomes
not used before.

### Fixed

- **The pooled rule could drop a whole gene** when the call it kept was a
  fragment: the selector's length filter then deleted it, although the other
  engine had a complete call. ndhA was the gene most often lost this way. The
  discarded call is now kept aside and tried if the kept one did not survive,
  in a probe of the real selector that must return every other feature
  unchanged; after the final QC, a rescued call flagged NEEDS_REVIEW or
  containing an in-frame stop is withdrawn. It was validated on development
  genomes not used before, against a criterion fixed in advance: no feature
  was lost or moved. Measured accuracy will be published with the manuscript,
  after the independent review of the scorer and protocol.
- **Two tests failed wherever the database was not unpacked in the working
  directory.** `test_ir_lacking_genome.py` and `test_relatives_ordering.py` passed
  the literal path `database/protein_db`, so in a fresh clone, a worktree, or after
  `plastanno fetch-db` (which installs under the user data directory) the first
  failed a check and the second skipped everything -- and told the user to run
  the very command that could not help. Both now find the database as the tool
  does, through `plastanno.paths.db_root()`.
- **The start-up banner said "Plastanno v2.0".** It now prints the installed
  version.
- **A CDS carried its gene symbol as /product.** Every CDS but rps12 was
  written as, for example, `/product="psbA"` rather than "photosystem II
  protein D1", in the `.gb`, the `.gff3` and the `.tbl` meant for NCBI
  submission, in every release since the first. GenBank records almost never
  use the bare symbol. The catalog held the protein names all along, but
  nothing read them. Each CDS now takes its protein name from the catalog
  once, before any writer runs, so the three formats cannot disagree; a
  product already set is kept. The catalog's chlN name is corrected from
  "photochlorophyllide" to "protochlorophyllide reductase ChlN subunit". tRNA
  and rRNA products were already correct and are unchanged.
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

### Removed

- **Dead code, found by a static and a traced audit** (`parked/dead_code/FINDINGS.md`
  in the development repository): eight definitions referenced nowhere and never
  executed (`coords.locus_span`, `coords._touching`, `feature.Exon`,
  `reconcile.jaccard_overlap`, `TRNAIdentity.same_family`,
  `engine_b._normalize_trna_name`, `ambiguity.describe`, `writers._spliced_len`);
  the empty `plastanno/ml` and `plastanno/utils` packages; and
  `scripts/viz/plastome_circular_map.py`, a copy of the package's map module that
  had fallen behind it (it lacked the title-escaping fix). The writer used to load
  that copy on any error importing the package module, which would have replaced
  a broken module with old code silently; it now loads the package module only,
  and a failure skips the map with a message. No output changes.

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
