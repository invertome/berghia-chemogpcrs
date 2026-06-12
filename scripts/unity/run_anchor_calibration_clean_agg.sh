#!/bin/bash
# run_anchor_calibration_clean_agg.sh — aggregate the DASH-free C3 verdicts.
#
# Submit gated on the 4-task array:
#   sbatch --dependency=afterany:<clean_array_jobid> \
#          scripts/unity/run_anchor_calibration_clean_agg.sh
#
#SBATCH --job-name=anchor_cal_clean_agg
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/anchor_cal_clean_agg-%j.out
#SBATCH --error=logs/anchor_cal_clean_agg-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"
# shellcheck source=/dev/null
source "${REPO_ROOT}/config.sh"

CAL="${RESULTS_DIR}/p5_phase1a_validation/anchor_calibration"
n=$(find "${CAL}/verdicts_clean" -name 'verdict_class_*.json' 2>/dev/null | wc -l)
echo "[clean_agg] aggregating ${n} DASH-free verdict(s) from verdicts_clean/"

python3 scripts/aggregate_anchor_verdicts.py \
    --in-dir "${CAL}/verdicts_clean" \
    --out-json "${CAL}/anchor_calibration_report.json" \
    --out-md "${CAL}/anchor_calibration_report.md"
echo "===== anchor_calibration_report.md (DASH-free) ====="
cat "${CAL}/anchor_calibration_report.md"
