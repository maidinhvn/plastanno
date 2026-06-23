# GĐ4 Head-to-head provenance — Plastanno vs PGA

Self-contained, reproducible comparison of the **frozen Plastanno tool** against
**PGA** (Qu et al. 2019) on a common land-plant subset, both scored with the same
gene-by-gene metric.

## What was run
| | |
|---|---|
| Plastanno tool code | `plastanno/` at commit **b67fc4d** (boundary fix + internal-stop QC) |
| Targets | `splits/h2h_subset.txt` — 249 land-plant plastomes (245 paired after PGA) |
| PGA references | `splits/h2h_refmap.tsv` — each target's closest DEV reference(s); 73% same-genus, 26% same-family, 1% same-order |
| Plastanno command | `python3 plastanno.py run <acc>.fasta -o … -t 1 --no-plot` |
| PGA command | `perl PGA.pl -r <refs> -t <target> -o <out> -f circular` |
| Scorer (both tools) | `scripts/benchmark/score_h2h.py` → `benchmark_gene_by_gene.py` (±60 bp ends + 0.6 sequence similarity) |
| Reproduce | `bash scripts/benchmark/run_h2h.sh && python3 scripts/benchmark/score_h2h.py` |

Predictions: `bench_runs/h2h/plast_out/<acc>/<acc>.gb` (Plastanno) and
`bench_runs/h2h/pga_out/<acc>.gb` (PGA); per-genome paired F1 in
`docs/Plastanno_GD4_h2h_paired.csv`.

## Result (micro F1, n=245 paired)
| Tool | CDS | tRNA | rRNA | **Global** |
|---|---|---|---|---|
| **Plastanno** (b67fc4d) | 92.4 | 90.2 | 96.6 | **92.0** |
| PGA | 89.7 | 83.7 | 83.6 | **87.8** |

Per-genome: **Plastanno wins 145**, ties 33, PGA wins 67 (of 245).

## Note on the boundary fix
The multi-exon splice refinement (commit b67fc4d) changes exon boundaries only
within the metric's ±60 bp tolerance, so the gene-by-gene F1 is identical whether
scored on pre-fix or post-fix predictions. This run uses the **post-fix frozen
tool** so the reported numbers match the released code exactly.
