#!/bin/bash
#SBATCH --job-name=build_curated_chemo
#SBATCH --partition=cpu
#SBATCH --time=30
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/build_curated_chemo-%j.out
#SBATCH --error=logs/build_curated_chemo-%j.err
#SBATCH --account=pi_pkatz_umass_edu
#
# build_curated_chemo.sh — fetch-then-build wrapper for the curated
# invertebrate chemoreceptor HMM library.
#
# Sequence:
#   1. fetch_curated_chemoreceptors.py reads manifest.tsv, pulls
#      sequences from NCBI Entrez into per-family FASTAs.
#   2. build_curated_chemo_hmms.sh runs MAFFT-DASH E-INS-i + hmmbuild
#      per family (>=3 seqs/family).
#
# Idempotent: both steps skip already-cached / already-built outputs.
#
# Bead -k0g, 2026-05-19.

set -eo pipefail
mkdir -p logs

source ~/.miniconda3/etc/profile.d/conda.sh
set +u
conda activate berghia-gpcr
set -u 2>/dev/null || true

cd /scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
source config.sh

echo "=== fetch_curated_chemoreceptors.py ==="
python3 scripts/fetch_curated_chemoreceptors.py \
    --manifest references/curated_chemoreceptors_lophotrochozoa/manifest.tsv \
    --output-dir references/curated_chemoreceptors_lophotrochozoa

echo
echo "=== build_curated_chemo_hmms.sh ==="
REFS_DIR=references/curated_chemoreceptors_lophotrochozoa \
OUT_DIR="${RESULTS_DIR}/classification/hmms_curated_chemo" \
CPUS="${SLURM_CPUS_PER_TASK:-4}" \
bash scripts/build_curated_chemo_hmms.sh

echo
echo "=== consolidated curated_chemo.hmm ==="
ls -la "${RESULTS_DIR}/classification/hmms_curated_chemo/curated_chemo.hmm" 2>&1
grep -c '^NAME' "${RESULTS_DIR}/classification/hmms_curated_chemo/curated_chemo.hmm" 2>&1 \
    | xargs -I {} echo "HMMs in consolidated library: {}"

echo "BUILD_CURATED_CHEMO_COMPLETE"
