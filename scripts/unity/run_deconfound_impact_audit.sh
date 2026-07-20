#!/bin/bash
#SBATCH --job-name=deconf_audit
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Re-derive family assignment + envelope membership with the SAME length
# deconfounding the production novelty channel uses (fusion_consensus.py
# --deconfound seq_len). See embedding_deconfound_impact_audit.py.

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
cd "$W"

OUT="$W/results/ranking/diagnostics"
mkdir -p "$OUT"

python3 scripts/unity/embedding_deconfound_impact_audit.py \
    --emb-dir "$W/results/ranking/embeddings" \
    --ref-labels "$W/references/anchors/anchor_set_PROD.tsv" \
    --candidate-fasta "$W/results/chemogpcrs/chemogpcrs_berghia_classA.fa" \
    --anchor-fasta "$W/references/anchors/anchor_set_PROD.fasta" \
    --diag-dir "$OUT" \
    --model-a proteinclip3b \
    --model-b protrek \
    --out-prefix "$OUT/deconfound_impact_audit"
