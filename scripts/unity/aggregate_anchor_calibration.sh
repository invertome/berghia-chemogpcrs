#!/bin/bash
# aggregate_anchor_calibration.sh — merge the C3 per-class verdicts into the
# final report once both the monolithic driver (class A, in verdicts/) and the
# B/C/F side-array (verdicts_bcf/) have finished.
#
# Submit AFTER both, gated on their terminal state:
#   sbatch --dependency=afterany:<monolith_jobid>:<bcf_array_jobid> \
#          scripts/unity/aggregate_anchor_calibration.sh
#
#SBATCH --job-name=anchor_cal_agg
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/anchor_cal_agg-%j.out
#SBATCH --error=logs/anchor_cal_agg-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"
# shellcheck source=/dev/null
source "${REPO_ROOT}/config.sh"
# shellcheck source=/dev/null
source "${REPO_ROOT}/functions.sh"

CAL="${RESULTS_DIR}/p5_phase1a_validation/anchor_calibration"
FINAL="${CAL}/verdicts_final"
mkdir -p "${FINAL}"

# Monolith verdicts first (class A, possibly a redundant B); side-array second so
# the definitely-complete B/C/F verdicts win on any tie.
cp -f "${CAL}/verdicts/"verdict_class_*.json     "${FINAL}/" 2>/dev/null || true
cp -f "${CAL}/verdicts_bcf/"verdict_class_*.json "${FINAL}/" 2>/dev/null || true

n=$(find "${FINAL}" -name 'verdict_class_*.json' | wc -l)
log "anchor_cal_agg: aggregating ${n} per-class verdict(s) from ${FINAL}"

python3 scripts/aggregate_anchor_verdicts.py \
    --in-dir "${FINAL}" \
    --out-json "${CAL}/anchor_calibration_report.json" \
    --out-md "${CAL}/anchor_calibration_report.md"

log "anchor_cal_agg: report -> ${CAL}/anchor_calibration_report.md"
