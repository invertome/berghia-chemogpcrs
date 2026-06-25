#!/bin/bash
# run_per_class_pools.sh — production per-class reference-pool build (the "P2"
# step stage 04 consumes from ${PER_CLASS_POOL_DIR}).
#
# Reads ANCHOR_OUTGROUP_CLASSES (config.sh) and builds refs_class_{A,B,C,F}.fa
# with per-class anchor tiers: classes listed in ANCHOR_OUTGROUP_CLASSES get
# out-group anchor tiers (1,2,3); the rest get tier-1 (in-group mollusc) only.
# With ANCHOR_OUTGROUP_CLASSES="" (the C3 decision) every class is tier-1 only.
#
# Idempotent: the builder skips if pools already exist (FORCE=1 to rebuild).
# The new-flag pass is EXPLICIT (never env-defaulted) so the calibration
# scripts' with-all-anchors pool is unaffected.
#
# Usage:  sbatch scripts/unity/run_per_class_pools.sh
#
#SBATCH --job-name=per_class_pools
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/per_class_pools-%j.out
#SBATCH --error=logs/per_class_pools-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"
source "${REPO_ROOT}/config.sh"
source "${REPO_ROOT}/functions.sh"

# --- Resolve production inputs (mirror stage 04 / the calibration script) ----
SCAN_GLOB="${SCAN_DERIVED_CANDIDATES_DIR}/*.chemo_candidates.fa"
CLASSIFY_DIR="$(dirname "${BERGHIA_CLASS_TSV}")"
CLASS_PHASE1A="${CLASSIFY_DIR}/class_phase1a.tsv"
BERGHIA_FA="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"
ANCHOR_FASTA="${REPO_ROOT}/references/anchors/anchor_set.fasta"
ANCHOR_TSV="${REPO_ROOT}/references/anchors/anchor_set.tsv"
CPUS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${PER_CLASS_POOL_DIR}" logs

# --- Pre-flight: required inputs must exist ----------------------------------
for f in "${CLASS_PHASE1A}" "${BERGHIA_CLASS_TSV}" "${BERGHIA_FA}" \
         "${ANCHOR_FASTA}" "${ANCHOR_TSV}"; do
    [ -s "$f" ] || { echo "ERROR: missing required input: $f" >&2; exit 1; }
done
# shellcheck disable=SC2086
if ! ls ${SCAN_GLOB} >/dev/null 2>&1; then
    echo "ERROR: no scan candidates match ${SCAN_GLOB}" >&2
    exit 1
fi

echo "[$(date '+%H:%M:%S')] per_class_pools: ANCHOR_OUTGROUP_CLASSES='${ANCHOR_OUTGROUP_CLASSES}'"
echo "  scan glob : ${SCAN_GLOB}"
echo "  class tsv : ${CLASS_PHASE1A}"
echo "  out dir   : ${PER_CLASS_POOL_DIR}"

# --- Build the per-class pools -----------------------------------------------
# REF_CLUSTER_IDENTITY / REF_FLOOR_PER_SPECIES / REF_CAP_PER_TAXON come from
# config.sh via env (the builder reads them as flag defaults).
FORCE_FLAG=""
[ "${FORCE:-0}" = "1" ] && FORCE_FLAG="--force"

# shellcheck disable=SC2086
python3 "${REPO_ROOT}/scripts/build_per_class_reference_pools.py" \
    --scan-fasta-glob "${SCAN_GLOB}" \
    --class-tsv "${CLASS_PHASE1A}" \
    --berghia-fasta "${BERGHIA_FA}" \
    --berghia-class-tsv "${BERGHIA_CLASS_TSV}" \
    --anchor-fasta "${ANCHOR_FASTA}" \
    --anchor-tsv "${ANCHOR_TSV}" \
    --anchor-outgroup-classes "${ANCHOR_OUTGROUP_CLASSES}" \
    --out-dir "${PER_CLASS_POOL_DIR}" \
    --threads "${CPUS}" \
    ${FORCE_FLAG}

echo "[$(date '+%H:%M:%S')] per_class_pools: DONE -> ${PER_CLASS_POOL_DIR}/refs_class_{A,B,C,F}.fa"
