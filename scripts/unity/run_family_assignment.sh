#!/bin/bash
#SBATCH --job-name=emb_family_assign
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# run_family_assignment.sh — calibrated non-chemoreceptor family assignment.
# Bead berghia-chemogpcrs-nsew, deliverable 4.
#
# Answers "are any Berghia candidates actually bioamine / opsin / peptide
# receptors rather than chemoreceptors?" by asking whether a candidate falls
# inside the envelope a family's OWN reference members occupy — not merely which
# prototype happens to be nearest (argmin labels all 790 as something).
#
# Reads the .npz directly, so it does NOT depend on rebuild_embedding_channel.sh
# having finished; the two can run concurrently. Runs once per consensus model
# so a call that only holds for one embedding space is visible as such.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source ./config.sh
set -u

EMBEDDINGS_DIR="${EMBEDDINGS_DIR:-${RESULTS_DIR}/ranking/embeddings}"
EMB_CONSENSUS_MODELS="${EMB_CONSENSUS_MODELS:-proteinclip3b protrek}"
EMB_REF_LABELS="${EMB_REF_LABELS:-${REFERENCE_DIR}/anchors/anchor_set_PROD.tsv}"
# seq_len sources for the MANDATORY length deconfounding (matches production's
# fusion_consensus.py --deconfound seq_len). References and candidates both
# need lengths so they share one residual scale.
EMB_CANDIDATE_FASTA="${EMB_CANDIDATE_FASTA:-${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia_classA.fa}"
EMB_ANCHOR_FASTA="${EMB_ANCHOR_FASTA:-${REFERENCE_DIR}/anchors/anchor_set_PROD.fasta}"
OUT_DIR="${RESULTS_DIR}/ranking/diagnostics"
mkdir -p "${OUT_DIR}" logs

case "${EMB_REF_LABELS}" in
  *anchor_set.tsv) echo "FATAL: 206-row label set — that is the defect." >&2; exit 2 ;;
esac

for _f in "${EMB_CANDIDATE_FASTA}" "${EMB_ANCHOR_FASTA}"; do
    [ -s "${_f}" ] || { echo "FATAL: missing seq_len source ${_f}; deconfounding cannot run." >&2; exit 2; }
done

for TAG in ${EMB_CONSENSUS_MODELS}; do
    echo
    echo "############ ${TAG} ############"
    python3 scripts/unity/embedding_family_assignment.py \
        --candidate-npz "${EMBEDDINGS_DIR}/candidates_${TAG}_classA.npz" \
        --ref-npz "${EMBEDDINGS_DIR}/reference_${TAG}_PROD.npz" \
        --ref-labels "${EMB_REF_LABELS}" \
        --candidate-fasta "${EMB_CANDIDATE_FASTA}" \
        --anchor-fasta "${EMB_ANCHOR_FASTA}" \
        --deconfound envelope \
        --out-tsv "${OUT_DIR}/family_assignment_${TAG}.tsv"
done

echo
echo "[assign] DONE $(date -Is)"
