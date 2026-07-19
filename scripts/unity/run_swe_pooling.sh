#!/bin/bash
# run_swe_pooling.sh — Sliced-Wasserstein (SWE) pooling of PER-RESIDUE protein
# embeddings into the embedding bake-off's per-model npz naming convention.
#
# WHY this script exists / where it plugs in:
#   scripts/swe_pooling.py (A5, epic v4bs) is an embedding-GENERATION transform,
#   not a stage-07 ranking producer. It replaces mean-pooling of a PLM's (L, D)
#   per-residue output with a distribution-aware, permutation-invariant,
#   length-robust pool. So its integration point is the EMBEDDING / BAKE-OFF
#   path, not any numbered pipeline stage: this wrapper simply emits
#       ${EMBEDDINGS_DIR}/candidates_${SWE_TAG}_classA.npz
#       ${EMBEDDINGS_DIR}/reference_${SWE_TAG}_PROD.npz
#   so a `swe`-pooled variant enters the existing bake-off harness as an
#   ordinary entry and is scored head-to-head against the mean-pooled models.
#   This script does NO ranking and promotes nothing to a voter.
#
# SELF-GATED: the per-residue npz inputs are produced GPU-side and do not exist
# yet. Each arm is skipped (with a "stays dormant" note) when its input is
# missing; the job exits 0 either way, so this is safe to run at any time.
#
# Usage:
#   sbatch scripts/unity/run_swe_pooling.sh
#   SWE_TAG=swe128 SWE_N_SLICES=128 sbatch scripts/unity/run_swe_pooling.sh
#   SWE_CAND_RESIDUE_NPZ=/path/cand_residue.npz \
#     SWE_REF_RESIDUE_NPZ=/path/ref_residue.npz \
#     sbatch scripts/unity/run_swe_pooling.sh
#
#SBATCH --job-name=swe_pooling
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO}"
source "${REPO}/config.sh"
source "${REPO}/functions.sh"

# --- Params (override via env) -----------------------------------------------
# EMBEDDINGS_DIR comes from config.sh (default ${RESULTS_DIR}/ranking/embeddings)
# and honours an env override.
SWE_TAG="${SWE_TAG:-swe}"                 # bake-off model tag for the pooled npz
SWE_N_SLICES="${SWE_N_SLICES:-64}"        # random projection directions
SWE_N_QUANTILES="${SWE_N_QUANTILES:-8}"   # quantiles read off each sorted slice
SWE_SEED="${SWE_SEED:-0}"                 # seeds the projection directions

# INPUTS — per-residue npz (protein_id -> (L, D) array), written by a GPU-side
# embedding run. These defaults name where such a run is expected to deposit
# them alongside the pooled bake-off npz; neither exists yet.
SWE_CAND_RESIDUE_NPZ="${SWE_CAND_RESIDUE_NPZ:-${EMBEDDINGS_DIR}/candidates_residue_classA.npz}"
SWE_REF_RESIDUE_NPZ="${SWE_REF_RESIDUE_NPZ:-${EMBEDDINGS_DIR}/reference_residue_PROD.npz}"

# OUTPUTS — the bake-off harness' per-model naming convention.
CAND_OUT_NPZ="${EMBEDDINGS_DIR}/candidates_${SWE_TAG}_classA.npz"
REF_OUT_NPZ="${EMBEDDINGS_DIR}/reference_${SWE_TAG}_PROD.npz"

mkdir -p "${EMBEDDINGS_DIR}" "${LOGS_DIR}"

log "swe_pooling: tag=${SWE_TAG} n_slices=${SWE_N_SLICES} n_quantiles=${SWE_N_QUANTILES} seed=${SWE_SEED}"

# --- Helpers -----------------------------------------------------------------
# Number of arrays in an npz, or "-" when unreadable/absent.
npz_count() {
    NPZ_PATH="$1" python3 - <<'PYEOF' 2>/dev/null || echo "-"
import os
import numpy as np
print(len(np.load(os.environ["NPZ_PATH"]).files))
PYEOF
}

# Pool one arm; skip (exit-0 semantics) when its per-residue input is missing.
pool_arm() {
    local label="$1" in_npz="$2" out_npz="$3"
    if [ ! -s "${in_npz}" ]; then
        log --level=WARN "swe_pooling: ${label} skipped — stays dormant (no per-residue npz at ${in_npz})"
        return 1
    fi
    log "swe_pooling: pooling ${label}: ${in_npz} -> ${out_npz}"
    python3 "${REPO}/scripts/swe_pooling.py" \
        --in-npz "${in_npz}" \
        --out-npz "${out_npz}" \
        --n-slices "${SWE_N_SLICES}" \
        --n-quantiles "${SWE_N_QUANTILES}" \
        --seed "${SWE_SEED}"
    return 0
}

# --- Run both arms independently ---------------------------------------------
cand_done=0
ref_done=0
pool_arm "candidates" "${SWE_CAND_RESIDUE_NPZ}" "${CAND_OUT_NPZ}" && cand_done=1
pool_arm "references" "${SWE_REF_RESIDUE_NPZ}" "${REF_OUT_NPZ}" && ref_done=1

if [ "${cand_done}" -eq 0 ] && [ "${ref_done}" -eq 0 ]; then
    log --level=WARN "swe_pooling: no per-residue inputs present — nothing pooled, channel stays dormant"
fi

# --- Provenance summary -------------------------------------------------------
echo "--- swe_pooling provenance ---"
echo "  tag         : ${SWE_TAG}"
echo "  n_slices    : ${SWE_N_SLICES}"
echo "  n_quantiles : ${SWE_N_QUANTILES}"
echo "  seed        : ${SWE_SEED}"
echo "  pooled dim  : $(( SWE_N_SLICES * SWE_N_QUANTILES ))"
if [ "${cand_done}" -eq 1 ]; then
    echo "  candidates  : ${CAND_OUT_NPZ} (n=$(npz_count "${CAND_OUT_NPZ}"))"
else
    echo "  candidates  : SKIPPED (missing ${SWE_CAND_RESIDUE_NPZ})"
fi
if [ "${ref_done}" -eq 1 ]; then
    echo "  references  : ${REF_OUT_NPZ} (n=$(npz_count "${REF_OUT_NPZ}"))"
else
    echo "  references  : SKIPPED (missing ${SWE_REF_RESIDUE_NPZ})"
fi

log "swe_pooling: DONE (candidates=${cand_done} references=${ref_done})"
exit 0
