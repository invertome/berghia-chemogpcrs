#!/bin/bash
# rerun_anchor_eval_bcf.sh — salvage the C3 calibration for classes B/C/F.
#
# The B/C/F trees (with_/without_) were built successfully by the side-array
# (job 60712714) but evaluate_anchor_divergence.py crashed on IQ-TREE's 3-part
# aBayes/SH-aLRT/UFBoot support (fixed in load_tree, commit 7685660). The trees
# still exist, so just re-run the evaluation (cheap) + re-aggregate. No tree
# rebuild. Class A is handled separately (its tree-build OOM'd at 64G).
#
# Submit:  sbatch scripts/unity/rerun_anchor_eval_bcf.sh
#
#SBATCH --job-name=anchor_eval_bcf
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/anchor_eval_bcf-%j.out
#SBATCH --error=logs/anchor_eval_bcf-%j.err
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

for CLS in B C F; do
    wt="${CAL}/trees_bcf/with_class_${CLS}.treefile"
    ot="${CAL}/trees_bcf/without_class_${CLS}.treefile"
    if [ ! -s "$wt" ] || [ ! -s "$ot" ]; then
        echo "[eval_bcf] class ${CLS}: missing tree(s), skipping"; continue
    fi
    echo "[eval_bcf] class ${CLS}: evaluating existing trees"
    python3 scripts/evaluate_anchor_divergence.py \
        --class "${CLS}" \
        --tree-without "$ot" \
        --tree-with "$wt" \
        --pool-members "${CAL}/pools/pool_members_class_${CLS}.tsv" \
        --out "${CAL}/verdicts_bcf/verdict_class_${CLS}.json"
done

# Re-aggregate (monolith class-A verdicts in verdicts/ + B/C/F in verdicts_bcf/).
FINAL="${CAL}/verdicts_final"
mkdir -p "${FINAL}"
cp -f "${CAL}/verdicts/"verdict_class_*.json     "${FINAL}/" 2>/dev/null || true
cp -f "${CAL}/verdicts_bcf/"verdict_class_*.json "${FINAL}/" 2>/dev/null || true
n=$(find "${FINAL}" -name 'verdict_class_*.json' | wc -l)
echo "[eval_bcf] aggregating ${n} verdict(s)"
python3 scripts/aggregate_anchor_verdicts.py \
    --in-dir "${FINAL}" \
    --out-json "${CAL}/anchor_calibration_report.json" \
    --out-md "${CAL}/anchor_calibration_report.md"
echo "===== anchor_calibration_report.md ====="
cat "${CAL}/anchor_calibration_report.md"
