#!/bin/bash
#SBATCH --job-name=align_filters_install
#SBATCH --partition=cpu
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# install_alignment_filters.sh — install PREQUAL, TAPER (Julia), CLOAK on Unity.
#
# Replaces the abandoned HmmCleaner install. Pipeline plan after this:
#   PREQUAL (input residue mask)
#     -> MAFFT canonical + 4 MAFFT variants + FAMSA (6-aln ensemble)
#     -> CLOAK consensus mask
#     -> TAPER residue-outlier
#     -> ClipKit column trim
#     -> IQ-TREE 3
#
# Tools:
#   PREQUAL — Whelan, Irisarri & Burki 2018 Bioinformatics 34:3929 — bioconda.
#   TAPER   — Zhang et al. 2021 MEE — Julia, source github chaoszhang/TAPER.
#   CLOAK   — Chatur (Wheeler lab) bioRxiv Dec 2025 — Python, source github
#             phylowheeler/CLOAK.
#
# All three install paths are independent; failures in one don't abort the
# others. End-of-script verification reports each tool's status.

set -eo pipefail
mkdir -p logs

WORKDIR="/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05"
TOOLS_DIR="${WORKDIR}/tools"
mkdir -p "${TOOLS_DIR}"

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate berghia-gpcr

echo "===================================================================="
echo "[$(date -u +%FT%TZ)] Installing PREQUAL + TAPER + CLOAK"
echo "===================================================================="
echo "Active env: ${CONDA_DEFAULT_ENV}"
echo "CONDA_PREFIX: ${CONDA_PREFIX}"
echo ""

PREQUAL_OK=0
TAPER_OK=0
CLOAK_OK=0

# -----------------------------------------------------------------------------
# 1. PREQUAL — bioconda
# -----------------------------------------------------------------------------
echo "--- 1/3: PREQUAL (bioconda) ---"
if command -v prequal >/dev/null 2>&1; then
    echo "PREQUAL already installed at $(command -v prequal); skipping mamba."
    PREQUAL_OK=1
else
    if mamba install -y -n berghia-gpcr -c bioconda -c conda-forge prequal; then
        if command -v prequal >/dev/null 2>&1; then
            echo "PREQUAL installed at $(command -v prequal)"
            PREQUAL_OK=1
        fi
    fi
fi
[[ ${PREQUAL_OK} -eq 1 ]] || echo "WARN: PREQUAL install failed; continuing."
echo ""

# -----------------------------------------------------------------------------
# 2. TAPER — Julia + clone repo
# -----------------------------------------------------------------------------
echo "--- 2/3: TAPER (Julia source clone) ---"
TAPER_DIR="${TOOLS_DIR}/TAPER"
if [[ ! -d "${TAPER_DIR}" ]]; then
    git clone https://github.com/chaoszhang/TAPER.git "${TAPER_DIR}" || \
        echo "WARN: git clone TAPER failed."
fi
if [[ ! -f "${TAPER_DIR}/correction_multi.jl" ]]; then
    echo "ERROR: TAPER repo missing correction_multi.jl after clone." >&2
else
    if ! command -v julia >/dev/null 2>&1; then
        echo "Installing julia via conda-forge..."
        mamba install -y -n berghia-gpcr -c conda-forge julia || \
            echo "WARN: julia install failed."
    fi
    if command -v julia >/dev/null 2>&1; then
        echo "julia at $(command -v julia) ($(julia --version))"
        if julia "${TAPER_DIR}/correction_multi.jl" -h >/dev/null 2>&1; then
            echo "TAPER --help: OK"
            TAPER_OK=1
        else
            echo "WARN: TAPER --help failed; may need first-run package compile."
        fi
    fi
fi
echo ""

# -----------------------------------------------------------------------------
# 3. CLOAK — Python single-file clone
# -----------------------------------------------------------------------------
echo "--- 3/3: CLOAK (Python source clone) ---"
CLOAK_DIR="${TOOLS_DIR}/CLOAK"
if [[ ! -d "${CLOAK_DIR}" ]]; then
    git clone https://github.com/phylowheeler/CLOAK.git "${CLOAK_DIR}" || \
        echo "WARN: git clone CLOAK failed."
fi
if [[ -f "${CLOAK_DIR}/cloak.py" ]]; then
    echo "CLOAK at ${CLOAK_DIR}/cloak.py"
    # cloak.py only imports os + sys (verified). No pip deps needed.
    if python3 -c "import os, sys" 2>/dev/null; then
        CLOAK_OK=1
    fi
fi
echo ""

# -----------------------------------------------------------------------------
# Verification
# -----------------------------------------------------------------------------
echo "===================================================================="
echo "[$(date -u +%FT%TZ)] Install summary"
echo "===================================================================="
printf "  PREQUAL : %s\n" "$([[ ${PREQUAL_OK} -eq 1 ]] && echo OK || echo FAILED)"
printf "  TAPER   : %s\n" "$([[ ${TAPER_OK}   -eq 1 ]] && echo OK || echo FAILED)"
printf "  CLOAK   : %s\n" "$([[ ${CLOAK_OK}   -eq 1 ]] && echo OK || echo FAILED)"
echo ""
echo "Tool paths to set in config.local.sh:"
[[ ${PREQUAL_OK} -eq 1 ]] && echo "  export PREQUAL=$(command -v prequal)"
[[ ${TAPER_OK}   -eq 1 ]] && echo "  export TAPER=${TAPER_DIR}/correction_multi.jl"
[[ ${TAPER_OK}   -eq 1 ]] && echo "  export JULIA=$(command -v julia)"
[[ ${CLOAK_OK}   -eq 1 ]] && echo "  export CLOAK=${CLOAK_DIR}/cloak.py"

# Exit non-zero only if ALL three failed (a partial install is still useful)
if [[ ${PREQUAL_OK} -eq 0 && ${TAPER_OK} -eq 0 && ${CLOAK_OK} -eq 0 ]]; then
    exit 1
fi
exit 0
