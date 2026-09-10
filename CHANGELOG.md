# Changelog

All notable changes to Plastanno v2. Benchmarks are measured on the DEV split
(n=123 shared genomes) against reference GenBank annotations; the held-out set is
never used during development.

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
