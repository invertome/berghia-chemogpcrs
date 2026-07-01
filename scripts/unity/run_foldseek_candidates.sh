#!/bin/bash
# run_foldseek_candidates.sh — Foldseek structural-evidence search (Task 4 of
# the ML/PLM chemoreceptor ranking plan,
# docs/plans/2026-07-01-ml-plm-chemoreceptor-ranking.md).
#
# Runs `foldseek easy-search` of each candidate's AlphaFold model (stage 08
# output) against three reference structure databases: PDB, AFDB50, GPCRdb.
# The resulting hits.tsv files feed scripts/structural_evidence.py, which
# turns them into two HONEST signals per candidate -- never a positive
# "looks like a known chemoreceptor" score:
#   - novelty/recall:  no confident structural hit anywhere (struct_novelty)
#   - exclusion:       confident hit to a known NON-chemoreceptor GPCR
#                       family (struct_nonchemo_corrob)
#
# CPU-only: Foldseek structural search (not structure prediction) does not
# need a GPU.
#
# Submit from the repo root:
#   sbatch scripts/unity/run_foldseek_candidates.sh
#
#SBATCH --job-name=foldseek_candidates
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --output=logs/foldseek_candidates-%j.out
#SBATCH --error=logs/foldseek_candidates-%j.err

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"
source "${REPO_ROOT}/config.sh"

# Candidate AlphaFold models (stage 08 output; one structure file per
# candidate, mmCIF or PDB -- foldseek ingests both natively).
AF_DIR="${AF_DIR:-${RESULTS_DIR}/structural_analysis/alphafold}"

# Reference structure databases, pre-built with `foldseek createdb` /
# `foldseek databases`. Sane defaults live under references/foldseek/;
# override per env if built elsewhere on Unity.
PDB_FOLDSEEK_DB="${PDB_FOLDSEEK_DB:-${REPO_ROOT}/references/foldseek/pdb}"
AFDB50_FOLDSEEK_DB="${AFDB50_FOLDSEEK_DB:-${REPO_ROOT}/references/foldseek/afdb50}"
GPCRDB_FOLDSEEK_DB="${GPCRDB_FOLDSEEK_DB:-${REPO_ROOT}/references/foldseek/gpcrdb_2025}"

OUT_DIR="${OUT_DIR:-${RESULTS_DIR}/foldseek/candidates}"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "${OUT_DIR}" "${TMP_DIR}" logs

FOLDSEEK_BIN="$(command -v "${FOLDSEEK:-foldseek}" || true)"
if [ -z "${FOLDSEEK_BIN}" ]; then
    echo "ERROR: foldseek not found (FOLDSEEK=${FOLDSEEK:-foldseek}). " \
         "Install via scripts/unity/install_structural_tools.sh." >&2
    exit 2
fi
echo "[foldseek-candidates] using foldseek: ${FOLDSEEK_BIN}"
echo "[foldseek-candidates] query AF models: ${AF_DIR}"

if [ ! -d "${AF_DIR}" ]; then
    echo "ERROR: candidate AlphaFold model dir not found: ${AF_DIR}" >&2
    echo "       (run stage 08_structural_analysis.sh first)" >&2
    exit 2
fi

declare -A DB_PATHS=(
    [PDB]="${PDB_FOLDSEEK_DB}"
    [AFDB50]="${AFDB50_FOLDSEEK_DB}"
    [GPCRdb]="${GPCRDB_FOLDSEEK_DB}"
)

for DB in PDB AFDB50 GPCRdb; do
    TARGET_DB="${DB_PATHS[$DB]}"
    if [ ! -e "${TARGET_DB}" ]; then
        echo "WARNING: ${DB} foldseek DB not found at ${TARGET_DB} -- skipping." \
             "Build via: foldseek createdb <structures/> ${TARGET_DB}" >&2
        continue
    fi
    OUT="${OUT_DIR}/${DB}.tsv"
    echo "[foldseek-candidates] searching vs ${DB}: ${TARGET_DB} -> ${OUT}"
    "${FOLDSEEK_BIN}" easy-search "${AF_DIR}" "${TARGET_DB}" "${OUT}" "${TMP_DIR}" \
        --format-output "query,target,fident,alntmscore,evalue" \
        --threads "${SLURM_CPUS_PER_TASK:-8}"
done

echo "[foldseek-candidates] DONE -> ${OUT_DIR}/{PDB,AFDB50,GPCRdb}.tsv"
