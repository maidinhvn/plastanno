#!/usr/bin/env bash
# Worked example: annotate the bundled plastomes with Plastanno v2.
# Requires BLAST+, Exonerate, HMMER (hmmsearch) and ARAGORN on PATH, plus the
# full database/ (see README "Databases").
set -e
cd "$(dirname "$0")/.."        # project root (so database/ is found)

# 1) Single genome -> out_example/NC_053537.1/ (6 files + a circular map)
python3 plastanno.py run example/NC_053537.1.fasta \
        --output out_example/NC_053537.1 --threads 4

# 2) Batch: annotate every FASTA in example/ -> out_example/<acc>/
python3 plastanno.py batch example --output out_example --threads 4

echo
echo "Done. See out_example/<accession>/<accession>.report for the gene summary."
