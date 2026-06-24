#!/usr/bin/env bash
# ============================================================================
# Run-once held-out evaluation of the FROZEN tool (git commit 8e88671).
# Runs the 10 stratified-by-family chunks sequentially, each to its own workdir
# and log. This is the final, unbiased generalisation test:
#   * the tool is frozen — do NOT change/tune anything after this run;
#   * run exactly once; if it crashes for an operational reason (env/path), fix
#     that and rerun, but never tune accuracy in response to the results.
#
# Usage (survives logout, ~5 h total):
#   nohup bash scripts/run_heldout.sh > bench_runs/heldout/run_all.log 2>&1 &
#
# Resume after interruption: already-finished chunks have a results.json in
# bench_runs/heldout/chunk_<i>/ ; this script skips those, so just rerun it.
# ============================================================================
set -u
# Ensure the plastanno conda env is active so external tools are on PATH (conda activate plastanno)
cd "$(dirname "$0")/.."
mkdir -p bench_runs/heldout

for i in 00 01 02 03 04 05 06 07 08 09; do
    wd="bench_runs/heldout/chunk_$i"
    if [ -f "$wd/results.json" ]; then
        echo "=== chunk_$i already done (results.json present) — skipping ==="
        continue
    fi
    echo "=== chunk_$i START $(date '+%F %T') ==="
    python3 scripts/benchmark/multi_genome_bench.py --set heldout --final \
        --acc-file "splits/heldout_chunks/chunk_$i.txt" \
        --workers 16 --keep --workdir "$wd" \
        > "bench_runs/heldout/chunk_$i.log" 2>&1
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "!!! chunk_$i FAILED (exit $rc) — see bench_runs/heldout/chunk_$i.log; stopping."
        exit $rc
    fi
    echo "=== chunk_$i DONE  $(date '+%F %T') ==="
done

echo "=========================================================="
echo "ALL 10 HELD-OUT CHUNKS DONE  $(date '+%F %T')"
echo "Aggregate with: python3 scripts/benchmark/aggregate_heldout.py"
echo "=========================================================="
