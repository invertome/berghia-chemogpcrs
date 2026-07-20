#!/bin/bash
#SBATCH --job-name=moll_cal_analyze
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Molluscan background null + rescore of the 790 + gate-enrichment caveat +
# phylogenetic-scope check (deliverables 1-3b), then the family-coherence
# measurement (deliverable 4), for both production models.
#
# CALIBRATION / VALIDATION ONLY. Nothing here is a production ranking and no
# family definition, threshold, or channel is changed by running it.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
cd "$W"
CAL="$W/results/ranking/diagnostics/molluscan_calibration"
E="$W/results/ranking/embeddings"
REFL="$W/references/anchors/anchor_set_PROD.tsv"
ANCHOR_FA="$W/references/anchors/anchor_set_PROD.fasta"
CAND_FA="$W/results/chemogpcrs/chemogpcrs_berghia_classA.fa"
mkdir -p "$CAL" logs

# add detection-evidence strata (adds columns only; the seq_id set is unchanged,
# so the FASTAs and the embeddings already built from them stay valid)
python3 scripts/unity/molluscan_calibration_augment_sample.py \
    --sample-tsv "$CAL/molluscan_calibration_sample.tsv" \
    --class-tsv  "$W/results/p5_phase1a_validation/classify/class_phase1a.tsv"

for TAG in ${MODELS:-proteinclip3b protrek}; do
    NULLNPZ="$CAL/molluscan_null_${TAG}.npz"
    [ -s "$NULLNPZ" ] || { echo "FATAL: missing $NULLNPZ" >&2; exit 2; }

    echo
    echo "##################### ${TAG} :: null + rescore #####################"
    python3 scripts/unity/molluscan_calibration_null.py \
        --model "$TAG" \
        --ref-npz       "$E/reference_${TAG}_PROD.npz" \
        --candidate-npz "$E/candidates_${TAG}_classA.npz" \
        --null-npz      "$NULLNPZ" \
        --ref-labels    "$REFL" \
        --candidate-fasta "$CAND_FA" \
        --anchor-fasta    "$ANCHOR_FA" \
        --null-fasta      "$CAL/molluscan_null_all.fa" \
        --sample-tsv      "$CAL/molluscan_calibration_sample.tsv" \
        --out-prefix      "$CAL/molluscan_null_${TAG}"

    echo
    echo "################## ${TAG} :: family coherence ######################"
    python3 scripts/unity/molluscan_calibration_family_coherence.py \
        --model "$TAG" \
        --ref-npz       "$E/reference_${TAG}_PROD.npz" \
        --candidate-npz "$E/candidates_${TAG}_classA.npz" \
        --ref-labels    "$REFL" \
        --subfamily-tsv "$CAL/molluscan_calibration_subfamilies.tsv" \
        --candidate-fasta "$CAND_FA" \
        --anchor-fasta    "$ANCHOR_FA" \
        --null-npz        "$NULLNPZ" \
        --null-fasta      "$CAL/molluscan_null_all.fa" \
        --out-prefix      "$CAL/family_coherence_${TAG}"
done

echo
echo "[analyze] DONE $(date -Is)"
