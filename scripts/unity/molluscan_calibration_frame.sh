#!/bin/bash
#SBATCH --job-name=moll_cal_frame
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Build the molluscan background-null sampling frame + the FASTAs to embed.
# Calibration/validation only -- never a production ranking.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
cd "$W"
OUT="$W/results/ranking/diagnostics/molluscan_calibration"
mkdir -p "$OUT" logs

python3 scripts/unity/molluscan_calibration_frame.py \
    --scan-dir      "$W/results/p5_phase1a_validation/scan" \
    --class-tsv     "$W/results/p5_phase1a_validation/classify/class_phase1a.tsv" \
    --proteome-dir  "$W/species_tree_data/ncbi_proteomes" \
    --taxonomy-tsv  "$OUT/molluscan_calibration_taxonomy.tsv" \
    --out-dir       "$OUT" \
    --n-postgate    "${N_POSTGATE:-6000}" \
    --seed          "${SEED:-20260720}"

echo "[frame] DONE $(date -Is)"
