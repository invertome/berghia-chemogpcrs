#!/bin/bash
#SBATCH --job-name=taar_analysis
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --account=pi_pkatz_umass_edu
#
# TAAR embedding job D (bead h0y0): the Q1 placement + Q2 add-impact analysis.
# Submitted with --dependency=afterok on the two embed arms, so it runs only
# after both proteinclip3b_all.npz and protrek_all.npz exist. Reuses the
# production scorer embedding_evidence.py; writes VERDICT.md (ADD/HOLD rule).
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"; conda activate berghia-gpcr; set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
cd "$W"
export PYTHONPATH="$W/scripts:${PYTHONPATH:-}"
T=results/ranking/taar_embedding
LAB=references/anchors/anchor_set_PROD.tsv
mkdir -p "$T/proteinclip3b" "$T/protrek" "$T/consensus" logs

for m in proteinclip3b protrek; do
    [ -s "$T/${m}_all.npz" ] || { echo "FATAL: missing $T/${m}_all.npz" >&2; exit 2; }
done
[ -f "$LAB" ] || { echo "FATAL: missing anchor labels $LAB" >&2; exit 2; }

echo "[D] per-model proteinclip3b $(date -Is)"
python3 scripts/analysis/taar_embedding_placement.py \
    --emb-npz "$T/proteinclip3b_all.npz" --model proteinclip3b \
    --anchor-labels "$LAB" --out-dir "$T/proteinclip3b"

echo "[D] per-model protrek $(date -Is)"
python3 scripts/analysis/taar_embedding_placement.py \
    --emb-npz "$T/protrek_all.npz" --model protrek \
    --anchor-labels "$LAB" --out-dir "$T/protrek"

echo "[D] consensus $(date -Is)"
python3 scripts/analysis/taar_embedding_placement.py \
    --consensus-from "$T/proteinclip3b" "$T/protrek" --out-dir "$T/consensus"

echo; echo "=== VERDICT (consensus) ==="
cat "$T/consensus/VERDICT.md" 2>/dev/null || echo "(no consensus VERDICT.md written)"
echo "[D] DONE $(date -Is)"
