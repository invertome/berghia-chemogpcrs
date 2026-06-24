#!/bin/bash
# run_landmark_placement_validation.sh — one-time pilot landmark-placement study (bead 521.7).
#
# Runs validate_landmark_placement.py over the C3 pilot data: places the dropped
# tier-2/3 out-group anchors back onto the clean (without-anchor) per-class trees
# with EPA-ng and runs three directional comparisons — functional-family read vs
# the HMM classifier (axis 1), placement position vs in-inference position (axis 2),
# and Berghia-restricted RF / LSE robustness (axis 3).
#
# Submit from the repo root:
#   sbatch scripts/unity/run_landmark_placement_validation.sh
#
#SBATCH --job-name=landmark_placement_validation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --output=logs/landmark_placement-%j.out
#SBATCH --error=logs/landmark_placement-%j.err

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"

python3 scripts/validate_landmark_placement.py \
    --calibration-dir results/p5_phase1a_validation/anchor_calibration \
    --anchor-tsv references/anchors/anchor_set.tsv \
    --anchor-fasta references/anchors/anchor_set.fasta \
    --class-berghia-tsv results/p5_phase1a_validation/classify/class_berghia.tsv \
    --out results/p5_phase1a_validation/landmark_placement \
    --classes "A B C F" \
    --threads "${SLURM_CPUS_PER_TASK:-4}"
