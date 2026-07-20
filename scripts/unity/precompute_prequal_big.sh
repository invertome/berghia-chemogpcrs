#!/bin/bash
#SBATCH --job-name=precompute_prequal
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=logs/precompute_prequal-%j.out
#SBATCH --error=logs/precompute_prequal-%j.err
#SBATCH --account=pi_pkatz_umass_edu
#
# precompute_prequal_big.sh — Pre-compute PREQUAL output for the large
# all_berghia_refs alignment input (~2800 seqs, max length ~1700 aa).
# PREQUAL's pairwise-posterior pass is O(n^2) and exceeds the per-task
# wallclock budget of stage 04's array jobs. By precomputing here and
# placing the output at the path the filter stack expects, stage 04
# array tasks reuse the cached PREQUAL output via the size-check in
# run_alignment_filter_stack.
#
# Per-OG inputs are small (<100 seqs) and run PREQUAL inline in stage
# 04's array tasks — no precompute needed for those.
#
# Usage: sbatch scripts/unity/precompute_prequal_big.sh

set -eo pipefail
mkdir -p logs

source ~/.miniconda3/etc/profile.d/conda.sh
set +u
conda activate berghia-gpcr
set -u 2>/dev/null || true

BASE_DIR="/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05"
cd "$BASE_DIR"
source config.sh

INPUT="${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_filtered.fa"
TAG="all_berghia_refs"
WORKDIR="${RESULTS_DIR}/phylogenies/protein/_filter_stack_global"
OUTPUT="${WORKDIR}/${TAG}_prequal.fa"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: input not found: $INPUT" >&2
    exit 1
fi

mkdir -p "$WORKDIR"

if [ -s "$OUTPUT" ]; then
    n=$(grep -c '^>' "$OUTPUT")
    echo "PREQUAL output already exists at $OUTPUT ($n seqs); nothing to do."
    exit 0
fi

echo "=== PREQUAL precompute for $TAG ==="
echo "Input:  $INPUT ($(grep -c '^>' "$INPUT") seqs)"
echo "Output: $OUTPUT"
echo "Threshold: ${PREQUAL_FILTERTHRESH:-0.90}"
date

PREQUAL_FILTERTHRESH="${PREQUAL_FILTERTHRESH:-0.90}" \
    bash "${SCRIPTS_DIR}/run_prequal.sh" \
    --input="$INPUT" --output="$OUTPUT"

date
echo "=== Done ==="
n_out=$(grep -c '^>' "$OUTPUT" 2>/dev/null || true)
n_out=${n_out:-0}
echo "PREQUAL kept $n_out sequences"
ls -la "$OUTPUT"
