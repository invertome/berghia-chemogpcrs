#!/bin/bash
#SBATCH --job-name=moll_cal_negctl
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# SAME-DIMENSION negative control for the geometry checker.
#
# The first negative control failed on dimensionality (2560 vs 128), which only
# proves the checker catches a gross space mismatch. This one differs from
# production in EXACTLY ONE step -- the L2 normalization inside
# ONNXModel.predict, a default argument no caller passes and therefore the step
# a reimplementation would silently drop -- and yields 128-d vectors either way,
# so dimensionality cannot rescue the check.
#
# The pathcheck MUST FAIL here. A pass means the norms / cosine instruments do
# not discriminate and a same-dimension divergence could reach a production null
# undetected.
#
# CPU only. Reuses the existing ESM-2 3B probe embeddings; no GPU forward pass.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate esmc
set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
PR=/scratch3/workspace/jperezmoreno_umass_edu-jorge/proteinclip
cd "$W"
E="$W/results/ranking/embeddings"
OUT="$W/results/ranking/diagnostics/molluscan_calibration/pathcheck"
mkdir -p "$OUT" logs

[ -s "$OUT/probe_esm2_3b.npz" ] || {
    echo "FATAL: $OUT/probe_esm2_3b.npz missing; run the pathcheck first" >&2
    exit 2; }

python3 scripts/unity/molluscan_calibration_negcontrol_norm.py \
    --base-npz          "$OUT/probe_esm2_3b.npz" \
    --out-npz           "$OUT/probe_proteinclip3b_NONORM.npz" \
    --proteinclip-repo  "$PR" \
    --head 36

echo
echo "##### same-dimension negative control: proteinclip3b WITHOUT input L2 norm #####"
if python3 scripts/unity/molluscan_calibration_pathcheck.py \
        --produced-npz   "$OUT/probe_proteinclip3b_NONORM.npz" \
        --production-npz "$E/candidates_proteinclip3b_classA.npz" \
        --ref-npz        "$E/reference_proteinclip3b_PROD.npz" \
        --ref-labels     "$W/references/anchors/anchor_set_PROD.tsv" \
        --label          "SAME-DIM NEGATIVE CONTROL (apply_norm=False)" \
        --out-json       "$OUT/pathcheck_negcontrol_nonorm.json"; then
    echo >&2
    echo "FATAL: the same-dimension negative control PASSED." >&2
    echo "The checker cannot detect a dropped L2 normalization, so it does NOT" >&2
    echo "protect against the realistic same-dimension divergence. The molluscan" >&2
    echo "null's geometry guarantee is weaker than claimed and must be revisited." >&2
    exit 3
fi
echo
echo "same-dimension negative control correctly FAILED (as required)"
echo "[negctl] DONE $(date -Is)"
