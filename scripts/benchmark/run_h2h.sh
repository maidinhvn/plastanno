#!/usr/bin/env bash
# ============================================================================
# GĐ4 head-to-head: Plastanno (frozen tool) vs PGA, self-contained & reproducible.
#
#   Tool version : git commit of plastanno/ at run time (record below)
#   Subset       : splits/h2h_subset.txt        (249 land-plant targets)
#   PGA refs     : splits/h2h_refmap.tsv         (per-target close DEV references;
#                  73% same-genus, 26% same-family, 1% same-order)
#   Scorer       : scripts/benchmark/score_h2h.py (same gene-by-gene scorer,
#                  +-60 bp ends + 0.6 similarity) used for BOTH tools.
#
# Reproduce:  bash scripts/benchmark/run_h2h.sh && python3 scripts/benchmark/score_h2h.py
# ============================================================================
set -u
# Ensure the plastanno conda env is active so external tools are on PATH (conda activate plastanno)
cd "$(dirname "$0")/../.."
EVAL=${PLASTANNO_DATA:-benchmark_data}/eval_v2_raw
RAW=${PLASTANNO_DATA:-benchmark_data}/rawdata
PGA="${PGA:-PGA.pl}"
H=bench_runs/h2h
mkdir -p "$H/plast_out" "$H/pga_out" "$H/pt"
echo "tool commit: $(git rev-parse HEAD)"  | tee "$H/PROVENANCE.txt"
echo "run date:    $(date -u '+%F %T UTC')" | tee -a "$H/PROVENANCE.txt"

# --- Plastanno (frozen tool), one process per genome, threads=1 ---
cat splits/h2h_subset.txt | xargs -P 24 -I{} bash -c \
  "python3 plastanno.py run $EVAL/{}.fasta -o $H/plast_out/{} -t 1 --no-plot >/dev/null 2>&1"

# --- PGA, per target with its closest DEV reference(s) (PGA's intended use) ---
run_pga () {
  local a="$1" d="$H/pt/$1"
  mkdir -p "$d/r" "$d/t" "$d/o"
  for r in $(awk -F'\t' -v a="$a" '$1==a{print $2}' splits/h2h_refmap.tsv | tr ',' ' '); do
    [ -f "$RAW/$r.gb" ] && cp "$RAW/$r.gb" "$d/r/"; done
  cp "$EVAL/$a.fasta" "$d/t/"
  perl "$PGA" -r "$d/r" -t "$d/t" -o "$d/o" -f circular >"$d/log" 2>&1
  [ -f "$d/o/$a.gb" ] && cp "$d/o/$a.gb" "$H/pga_out/$a.gb"
}
export -f run_pga; export RAW EVAL PGA H PATH
cat splits/h2h_subset.txt | xargs -P 24 -I{} bash -c 'run_pga "$@"' _ {}

echo "Plastanno: $(ls $H/plast_out/*/*.gb 2>/dev/null|wc -l) | PGA: $(ls $H/pga_out/*.gb 2>/dev/null|wc -l)"
echo "Now score:  python3 scripts/benchmark/score_h2h.py"
