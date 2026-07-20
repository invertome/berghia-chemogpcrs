#!/bin/bash
#SBATCH --job-name=moll_cal_pathproj
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# STAGE 1 of the geometry proof (CPU only): does the projection step my null
# pipeline uses reproduce the PRODUCTION proteinclip3b npz bit for bit?
#
# Takes the production ESM-2 3B candidate embeddings, pushes them through the
# exact projection command scripts/unity/molluscan_calibration_merge.sh runs,
# and compares the result against the production proteinclip3b npz. This
# isolates step 2+3 of the derivation (L2 norm + proteinclip_esm2_36.onnx) from
# step 1 (the GPU forward pass), which is checked separately by
# molluscan_calibration_pathcheck.sh.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate esmc
set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
PR=/scratch3/workspace/jperezmoreno_umass_edu-jorge/proteinclip
cd "$W"
E="$W/results/ranking/embeddings"
CAL="$W/results/ranking/diagnostics/molluscan_calibration"
OUT="$CAL/pathcheck"
mkdir -p "$OUT" logs

# identical invocation to molluscan_calibration_merge.sh and to
# scripts/unity/run_proteinclip3b.sh, which produced the production npz
PROTEINCLIP_REPO="$PR" PC_HEAD=36 \
    BASE_NPZ="$E/candidates_esm2_3b_classA.npz" \
    OUT_NPZ="$OUT/reproduced_proteinclip3b_from_production_esm2.npz" \
    python3 scratch_proteinclip_project_param.py

python3 scripts/unity/molluscan_calibration_pathcheck.py \
    --produced-npz   "$OUT/reproduced_proteinclip3b_from_production_esm2.npz" \
    --production-npz "$E/candidates_proteinclip3b_classA.npz" \
    --ref-npz        "$E/reference_proteinclip3b_PROD.npz" \
    --ref-labels     "$W/references/anchors/anchor_set_PROD.tsv" \
    --label          "projection-only (production ESM-2 3B -> proteinclip3b)" \
    --out-json       "$OUT/pathcheck_projection.json"

echo "[pathproj] DONE $(date -Is)"
