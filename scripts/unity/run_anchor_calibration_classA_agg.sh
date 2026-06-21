#!/bin/bash
# run_anchor_calibration_classA_agg.sh — compute the class-A verdict (once both
# split trees exist) and regenerate the 4-class calibration report.
#
# The B/C/F verdicts were produced by the 48h clean array; class A is built by
# run_anchor_calibration_classA_split.sh (2-task array, with/without). This step
# runs evaluate_anchor_divergence.py for class A, then re-aggregates ALL verdicts
# in verdicts_clean/ into anchor_calibration_report.{json,md}.
#
# Submit gated on the class-A split array:
#   sbatch --dependency=afterok:<split_array_jobid> \
#          scripts/unity/run_anchor_calibration_classA_agg.sh
#
#SBATCH --job-name=anchor_cal_A_agg
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/anchor_cal_A_agg-%j.out
#SBATCH --error=logs/anchor_cal_A_agg-%j.err
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
TREES_DIR="${CAL}/trees_${CAL_SUBDIR:-clean}"
VERD_DIR="${CAL}/verdicts_${CAL_SUBDIR:-clean}"
# Report basename: keep the default name for the canonical 'clean' run; suffix it
# for alternate runs (e.g. CAL_SUBDIR=nocloak -> anchor_calibration_report_nocloak).
REPORT_BASE="${CAL}/anchor_calibration_report"
[ "${CAL_SUBDIR:-clean}" != "clean" ] && REPORT_BASE="${REPORT_BASE}_${CAL_SUBDIR}"

# --- Class A verdict (needs BOTH split trees) ------------------------------
with_tree="${TREES_DIR}/with_class_A.treefile"
without_tree="${TREES_DIR}/without_class_A.treefile"
if [ -s "$with_tree" ] && [ -s "$without_tree" ]; then
    echo "[classA_agg] computing class-A verdict"
    python3 scripts/evaluate_anchor_divergence.py \
        --class A \
        --tree-without "$without_tree" \
        --tree-with "$with_tree" \
        --pool-members "${CAL}/pools/pool_members_class_A.tsv" \
        --out "${VERD_DIR}/verdict_class_A.json"
else
    echo "[classA_agg] WARN: class-A tree(s) missing (with=$([ -s "$with_tree" ] && echo ok || echo MISSING), without=$([ -s "$without_tree" ] && echo ok || echo MISSING)); verdict NOT computed"
fi

# --- Re-aggregate all verdicts into the 4-class report ---------------------
n=$(find "${VERD_DIR}" -name 'verdict_class_*.json' 2>/dev/null | wc -l)
echo "[classA_agg] aggregating ${n} verdict(s) from $(basename "${VERD_DIR}")/"
python3 scripts/aggregate_anchor_verdicts.py \
    --in-dir "${VERD_DIR}" \
    --out-json "${REPORT_BASE}.json" \
    --out-md "${REPORT_BASE}.md"
echo "===== $(basename "${REPORT_BASE}").md (4-class) ====="
cat "${REPORT_BASE}.md"
