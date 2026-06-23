# Plastanno v2 Design Document

## Lessons from v1

### Pipeline design
- Search BOTH IRs directly → no mirroring needed
- rRNA: detect in full genome → no IR duplicate
- tRNA: ARAGORN (de novo) + BLAST hierarchical DB + exon DB
- CDS: Exonerate with pre-built protein DB (not on-the-fly)

### Database design
- Build ALL databases with 1 script (build_all.py)
- tRNA hierarchical: genus (≥3 genomes) > family (≥5) > global
- tRNA exon DB: actual exon sequences from GenBank features
- rRNA: full-length only (1400-1600bp for rrn16 etc.)
- rRNA: 25000bp binning to avoid duplicates

### Gene handling
- Synonyms: clpP1→clpP, psbN→pbf1
- CAU disambiguation: trnI-CAU (IR) vs trnM-CAU (LSC)
- Multi-exon CDS: special cases for petB(6bp), petD(9bp)
- Trans-spliced: rps12 (3 exons across genome)
- Pseudogene: ycf1 at IRb/SSC junction

### ML integration (v2.0)
- ML guides Exonerate search regions (not standalone)
- Train after pipeline is complete and working
- Validate ML improvement with ablation study

## Benchmark methodology
- Gene-by-gene comparison (not count-based)
- Tolerance: ±60bp
- Separate metrics: CDS, rRNA, tRNA
- Compare with GeSeq and PGA on same dataset
