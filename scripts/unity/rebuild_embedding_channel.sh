#!/bin/bash
#SBATCH --job-name=emb_channel_rebuild
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# rebuild_embedding_channel.sh — rebuild the production embedding/novelty channel
# with the corrected reference label set. Bead berghia-chemogpcrs-nsew (P0).
#
# THIS IS A SCORING REBUILD, NOT A GPU JOB. All four .npz already exist under
# ${EMBEDDINGS_DIR}; nothing here runs a model. The work is matrix algebra over
# ~790 candidate x ~953 reference vectors at dim ~1250 — the largest object is a
# ~1250x1250 Ledoit-Wolf covariance, i.e. tens of MB. cpu/2 cores/16G is
# therefore already generous, and stays inside a standard node's memory so the
# job is schedulable on the widest possible pool. --time=2h is ~25x the expected
# few-minute runtime, per the standing 2-3x-minimum rule.
#
# WHAT WAS BROKEN. config.sh defaulted EMB_REF_LABELS to the June
# references/anchors/anchor_set.tsv (206 rows) while the reference embeddings are
# reference_<tag>_PROD.npz, built from the 1094-row PROD set whose class-A FASTA
# carries 953 composite ANCHOR_<class>_<tier>_<acc> headers. load_ref_labels
# rebuilds that composite key from BOTH files, so the key FORMAT was right on
# both sides and the join succeeded silently while resolving 116/953 = 12% of the
# references. 837 reference vectors were discarded with no error. Family
# prototypes collapsed (aminergic 37 vs 212, opsin 17 vs 105, peptide 39 vs 386,
# lipid 1 vs 70) and chemokine / nucleotide / glycoprotein-hormone got NO
# prototype at all — so a candidate that genuinely IS one of those had nothing to
# be near and scored MAXIMALLY NOVEL, ranking to the top of the wet-lab
# shortlist. That is the exact inversion of the exclusion axis's purpose.
#
# Fix: config.sh now defaults EMB_REF_LABELS to anchor_set_PROD.tsv (953/953).
# This wrapper re-derives the channel and PROVES the fix via three hard guards
# (see scripts/unity/embedding_channel_preflight.py) — reference-label overlap,
# deconfound key overlap, and degenerate zero-vector prototypes. Any guard that
# fires aborts the job non-zero rather than emitting a plausible artifact.
#
# Usage:  sbatch scripts/unity/rebuild_embedding_channel.sh
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr

cd "${SLURM_SUBMIT_DIR:-$PWD}"
# config.sh is sourced BEFORE `set -u` (it is not nounset-safe, and conda's
# deactivate hooks reference CONDA_BACKUP_* unguarded).
source ./config.sh
set -u

# SELF-SUFFICIENT DEFAULTS. This checkout is ~10 commits behind main and its
# config.sh predates the whole EMB_* block — it defines RESULTS_DIR and
# REFERENCE_DIR and nothing else this job needs. Syncing config.sh is NOT safe
# right now: two long A1 jobs are live and source it. So the embedding settings
# are declared here, matching main's config.sh values exactly, and every one
# stays env-overridable. When the checkout is synced these become no-ops.
EMB_SCORER="${EMB_SCORER:-consensus}"
EMBEDDINGS_DIR="${EMBEDDINGS_DIR:-${RESULTS_DIR}/ranking/embeddings}"
EMB_CONSENSUS_MODELS="${EMB_CONSENSUS_MODELS:-proteinclip3b protrek}"
EMB_CANDIDATE_FASTA="${EMB_CANDIDATE_FASTA:-${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia_classA.fa}"
# THE FIX (config.sh:554 on main): the PROD label set, not the 206-row June one.
EMB_REF_LABELS="${EMB_REF_LABELS:-${REFERENCE_DIR}/anchors/anchor_set_PROD.tsv}"
EMB_CANDIDATE_IDENTITY_TSV="${EMB_CANDIDATE_IDENTITY_TSV:-${EMBEDDINGS_DIR}/candidate_ref_identity_PROD.tsv}"
EMB_DIAG_DIR="${EMB_DIAG_DIR:-${RESULTS_DIR}/ranking/diagnostics}"

OUT_DIR="${RESULTS_DIR}/ranking/channels"
DIAG_DIR="${EMB_DIAG_DIR}"
mkdir -p "${OUT_DIR}" "${DIAG_DIR}" logs

echo "[rebuild] host=$(hostname) job=${SLURM_JOB_ID:-local} $(date -Is)"
echo "[rebuild] EMB_SCORER      = ${EMB_SCORER}"
echo "[rebuild] EMB_REF_LABELS  = ${EMB_REF_LABELS}"
echo "[rebuild] CANDIDATE_FASTA = ${EMB_CANDIDATE_FASTA}"
echo "[rebuild] EMBEDDINGS_DIR  = ${EMBEDDINGS_DIR}"

# Refuse the known-wrong label file outright, so a stale environment override
# cannot silently reintroduce the defect this job exists to fix.
case "${EMB_REF_LABELS}" in
  *anchor_set.tsv)
    echo "FATAL: EMB_REF_LABELS points at the 206-row June anchor_set.tsv." >&2
    echo "       That is the defect (12% reference overlap). Use anchor_set_PROD.tsv." >&2
    exit 2 ;;
esac
[ -f "${EMB_REF_LABELS}" ]       || { echo "FATAL: missing ${EMB_REF_LABELS}" >&2; exit 2; }
[ -f "${EMB_CANDIDATE_FASTA}" ]  || { echo "FATAL: missing ${EMB_CANDIDATE_FASTA}" >&2; exit 2; }

# Build the model specs exactly as 07_candidate_ranking.sh does.
SPECS=()
for TAG in ${EMB_CONSENSUS_MODELS}; do
    CAND="${EMBEDDINGS_DIR}/candidates_${TAG}_classA.npz"
    REF="${EMBEDDINGS_DIR}/reference_${TAG}_PROD.npz"
    case "${CAND}${REF}" in
      *_FN.npz*) echo "FATAL: refusing a _FN.npz input (corrupt sequences)" >&2; exit 2 ;;
    esac
    [ -f "${CAND}" ] || { echo "FATAL: missing ${CAND}" >&2; exit 2; }
    [ -f "${REF}" ]  || { echo "FATAL: missing ${REF}" >&2; exit 2; }
    # The identity TSV is the OPTIONAL 4th spec field (an extra confound, itself
    # dropped by the independence gate below 50% coverage). Append it only when
    # it exists: embedding_candidate_diagnostics._read_identity returns {} for an
    # EMPTY path but raises FileNotFoundError for a MISSING one, so passing a
    # path that isn't there kills the whole consensus producer. Observed on job
    # 61993481 — candidate_ref_identity_PROD.tsv has never been generated in this
    # checkout. NOTE: 07_candidate_ranking.sh appends this field
    # UNCONDITIONALLY, so in production the same crash is swallowed by its
    # `|| log --level=WARN` and the channel silently goes dormant.
    if [ -f "${EMB_CANDIDATE_IDENTITY_TSV}" ]; then
        SPECS+=("${TAG}:${CAND}:${REF}:${EMB_CANDIDATE_IDENTITY_TSV}")
    else
        echo "[rebuild] note: identity TSV absent (${EMB_CANDIDATE_IDENTITY_TSV}) -- omitting the optional confound field for ${TAG}"
        SPECS+=("${TAG}:${CAND}:${REF}")
    fi
done
[ "${#SPECS[@]}" -ge 2 ] || { echo "FATAL: need >=2 models, got ${#SPECS[@]}" >&2; exit 2; }
echo "[rebuild] model specs: ${SPECS[*]}"

# ---------------------------------------------------------------- PREFLIGHT --
# Hard guards. Non-zero exit here aborts the job BEFORE any artifact is written,
# so a voided run can never be mistaken for a result.
echo
echo "=== PREFLIGHT (3 hard guards) ==="
python3 scripts/unity/embedding_channel_preflight.py \
    --ref-labels "${EMB_REF_LABELS}" \
    --candidate-fasta "${EMB_CANDIDATE_FASTA}" \
    --models "${SPECS[@]}" \
    --out-json "${OUT_DIR}/embedding_preflight.json"

# --------------------------------------------------------------- PRODUCTION --
# The consensus channel is what stage 07 writes; arguments mirror
# 07_candidate_ranking.sh's invocation verbatim (combiner rra, deconfound seq_len).
echo
echo "=== CONSENSUS CHANNEL (production: RRA over ${EMB_CONSENSUS_MODELS}) ==="
python3 scripts/fusion_consensus.py \
    --models "${SPECS[@]}" \
    --ref-labels "${EMB_REF_LABELS}" \
    --candidate-fasta "${EMB_CANDIDATE_FASTA}" \
    --combiner rra --deconfound seq_len \
    --diag-dir "${DIAG_DIR}" \
    --out "${OUT_DIR}/embedding_channel.tsv"

# The consensus schema carries emb_nonchemo_family but NOT emb_nonchemo_sim
# (CONSENSUS_TSV_COLUMNS). The cosine scorer is the only path that emits a
# similarity, and the stale 888-row calibration figures came from it — so run it
# too, purely as a diagnostic, to make the saturation question answerable.
# This does NOT overwrite the production channel.
echo
echo "=== COSINE EXCLUSION CHANNEL (diagnostic: emb_nonchemo_sim) ==="
FIRST_TAG="${EMB_CONSENSUS_MODELS%% *}"
python3 scripts/build_embedding_channel.py \
    --candidate-npz "${EMBEDDINGS_DIR}/candidates_${FIRST_TAG}_classA.npz" \
    --ref-npz "${EMBEDDINGS_DIR}/reference_${FIRST_TAG}_PROD.npz" \
    --ref-labels "${EMB_REF_LABELS}" \
    --out "${OUT_DIR}/embedding_channel_cosine_${FIRST_TAG}.tsv"

# ------------------------------------------------------------------ REPORT ---
echo
echo "=== REPORT ==="
python3 scripts/unity/embedding_channel_report.py \
    --consensus-tsv "${OUT_DIR}/embedding_channel.tsv" \
    --cosine-tsv "${OUT_DIR}/embedding_channel_cosine_${FIRST_TAG}.tsv" \
    --top 30

echo
echo "[rebuild] DONE $(date -Is)"
