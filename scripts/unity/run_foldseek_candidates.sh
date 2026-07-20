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

# Candidate AlphaFold models (stage 08 output; mmCIF or PDB -- foldseek ingests
# both natively).
#
# Bead 5ubd: stage 08 nests AF3 output as
#     alphafold/<candidate_id>/<af3_job_name>/<af3_job_name>_model.cif
# (08_structural_analysis.sh:298), so AF_DIR itself holds only per-candidate
# DIRECTORIES. Handing it straight to `foldseek easy-search` searches nothing
# usable, and any query id foldseek did derive would be the AF3 job name rather
# than the candidate id -- build_structural_channel.py uses the foldseek `query`
# field VERBATIM as the channel join key, so the whole structural voter would
# join against nothing. Stage the models into a flat query dir instead.
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

# --- Bead 5ubd: flat query dir, one symlink per candidate --------------------
# Each model is linked in as <candidate_id>.<ext> so foldseek's query id carries
# the candidate id. The candidate id is the FIRST path component under AF_DIR --
# the `cand_dir` stage 08 created (08:222) and reads back at 08:303 -- never the
# file stem, which is AF3's internal job name.
#
# Deliberately NOT results/structural_analysis/all_pdb/, even though it is
# already flat: stage 08 pools GPCRdb REFERENCE structures into that directory
# alongside the predictions (08:274). Searching it would emit channel rows keyed
# on reference PDB accessions -- non-candidates scored as candidates -- and the
# flattened names are AF3 job names, so the candidate id is gone either way.
QUERY_DIR="${TMP_DIR}/query_structures"
rm -rf -- "${QUERY_DIR}"
mkdir -p "${QUERY_DIR}"

N_MODELS=0
while IFS= read -r _model; do
    [ -n "${_model}" ] || continue
    _rel="${_model#"${AF_DIR%/}"/}"
    case "${_rel}" in
        */*) _cand_id="${_rel%%/*}" ;;                 # <cand_id>/.../<file>
        *)   _cand_id="$(basename "${_rel%.*}")" ;;    # already flat
    esac
    _ext="${_model##*.}"
    ln -sfn "$(readlink -f "${_model}")" "${QUERY_DIR}/${_cand_id}.${_ext}"
    N_MODELS=$((N_MODELS + 1))
done < <(find "${AF_DIR}" \( -name '*_model.cif' -o -name '*_model.pdb' \) | sort)

if [ "${N_MODELS}" -eq 0 ]; then
    echo "ERROR: no candidate structures found under ${AF_DIR}" >&2
    echo "       (searched recursively for *_model.cif / *_model.pdb)" >&2
    echo "       Stage 08 writes alphafold/<candidate_id>/<af3_job>/<af3_job>_model.cif;" >&2
    echo "       run 08_structural_analysis.sh first, or point AF_DIR at the tree that holds them." >&2
    echo "       Refusing to run a search over zero queries -- that is what made the" >&2
    echo "       structural ranking axis silently dormant (bead 5ubd)." >&2
    exit 3
fi
echo "[foldseek-candidates] staged ${N_MODELS} candidate structures -> ${QUERY_DIR}"

declare -A DB_PATHS=(
    [PDB]="${PDB_FOLDSEEK_DB}"
    [AFDB50]="${AFDB50_FOLDSEEK_DB}"
    [GPCRdb]="${GPCRDB_FOLDSEEK_DB}"
)

N_DBS_SEARCHED=0
for DB in PDB AFDB50 GPCRdb; do
    TARGET_DB="${DB_PATHS[$DB]}"
    if [ ! -e "${TARGET_DB}" ]; then
        echo "WARNING: ${DB} foldseek DB not found at ${TARGET_DB} -- skipping." \
             "Build via: foldseek createdb <structures/> ${TARGET_DB}" >&2
        continue
    fi
    OUT="${OUT_DIR}/${DB}.tsv"
    echo "[foldseek-candidates] searching vs ${DB}: ${TARGET_DB} -> ${OUT}"
    "${FOLDSEEK_BIN}" easy-search "${QUERY_DIR}" "${TARGET_DB}" "${OUT}" "${TMP_DIR}" \
        --format-output "query,target,fident,alntmscore,evalue" \
        --threads "${SLURM_CPUS_PER_TASK:-8}"
    N_DBS_SEARCHED=$((N_DBS_SEARCHED + 1))
done

# Bead 5ubd: every DB being absent means stage 07 finds no hit TSVs and logs the
# structural channel as "stays dormant". That must not look like success.
if [ "${N_DBS_SEARCHED}" -eq 0 ]; then
    echo "ERROR: none of the PDB/AFDB50/GPCRdb foldseek DBs were found, so no" >&2
    echo "       search ran. Build at least one (foldseek databases / createdb)" >&2
    echo "       or set PDB_FOLDSEEK_DB / AFDB50_FOLDSEEK_DB / GPCRDB_FOLDSEEK_DB." >&2
    exit 3
fi

echo "[foldseek-candidates] DONE: ${N_MODELS} candidates vs ${N_DBS_SEARCHED} DB(s) -> ${OUT_DIR}"
