#!/bin/bash
# test_models_gpcrtm_c20.sh — HPC-only substitution-model retest.
#
# Bead -5b0. Re-runs ModelFinder on a stratified subsample of the global
# Berghia+refs alignment with GPCRtm (Rios 2015, BMC Bioinf 16:206) and
# the LG profile mixture series (LG+C20+R10, LG+C60+R10, LG4X+R10) added
# to the candidate set. Reports BIC; the user picks the winner if ΔBIC > 10.
#
# Expected runtime: 24-72 hours on 16 CPUs / 64 GB RAM for ~2400 taxa.
# This SLURM script is NOT auto-submitted by the pipeline; run only when
# you want to validate / improve the substitution model.

#SBATCH --job-name=model_retest
#SBATCH --output=logs/test_models_gpcrtm_c20_%j.out
#SBATCH --error=logs/test_models_gpcrtm_c20_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$PROJECT_DIR/config.sh"
source "$PROJECT_DIR/functions.sh"

ALIGNMENT="${1:-${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa}"
OUT_DIR="${RESULTS_DIR}/model_retest"
mkdir -p "$OUT_DIR"

if [ ! -f "$ALIGNMENT" ]; then
    echo "ERROR: alignment not found at $ALIGNMENT" >&2
    exit 1
fi

# GPCRtm matrix file: must be downloaded separately. Source:
# https://www.iqtree.org/doc/Substitution-Models#protein-models
# (or ask the maintainers for the matrix file)
GPCRTM_MATRIX="${GPCRTM_MATRIX:-${PROJECT_DIR}/references/GPCRtm.txt}"
EXTRA_MADD=""
if [ -f "$GPCRTM_MATRIX" ]; then
    EXTRA_MADD="GPCRtm,LG+C20+R10,LG+C60+R10,LG4X+R10"
    GPCRTM_ARG="-mfile $GPCRTM_MATRIX"
else
    echo "INFO: GPCRtm matrix not at $GPCRTM_MATRIX; testing only LG+C20/C60/LG4X profile mixtures." >&2
    EXTRA_MADD="LG+C20+R10,LG+C60+R10,LG4X+R10"
    GPCRTM_ARG=""
fi

log "Running ModelFinder retest (-mset LG,VT,WAG,GPCRtm -madd $EXTRA_MADD)"

# Use ModelFinder only (no tree inference) — saves wall time.
# Full set including profile mixtures.
$IQTREE -s "$ALIGNMENT" \
    -m TESTONLY \
    -mset "${IQTREE_MODEL_SET},GPCRtm" \
    -madd "$EXTRA_MADD" \
    $GPCRTM_ARG \
    -seed "${IQTREE_SEED}" \
    -T "${SLURM_CPUS_PER_TASK:-${CPUS}}" \
    --prefix "${OUT_DIR}/model_retest"

log "ModelFinder retest complete. See ${OUT_DIR}/model_retest.iqtree for BIC table."
log "If a profile-mixture model wins by ΔBIC > 10, rebuild the global tree under that model."
