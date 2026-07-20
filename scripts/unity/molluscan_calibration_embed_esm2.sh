#!/bin/bash
#SBATCH --job-name=moll_cal_esm2
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=a100|a40|l40s|rtx8000
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=0-7%4
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#
# ESM-2 3B base embeddings for the molluscan background null.
#
# Same base + same embedder the production proteinclip3b channel is built from
# (facebook/esm2_t36_3B_UR50D via scratch_hf_auto_embed.py, MODE=esm,
# mean-pooled, MAXLEN=1024) -- comparing embeddings produced any other way
# would be meaningless, so this reuses the production embedder verbatim rather
# than reimplementing it.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate esmc
set -u

# The 11 GB ESM-2 3B checkpoint lives in the shared scratch HF cache, NOT in
# $HOME/.cache. Without this export every array task would try to re-download
# it (matches run_candidate_embed_all.sh, which produced the production npz).
export HF_HOME=/scratch3/workspace/jperezmoreno_umass_edu-jorge/hf_cache
export HF_HUB_DISABLE_TELEMETRY=1
export TORCH_HOME=/scratch3/workspace/jperezmoreno_umass_edu-jorge/torch_hub

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
cd "$W"
CAL="$W/results/ranking/diagnostics/molluscan_calibration"
SH="$CAL/shards"
NSHARD="${NSHARD:-8}"
IN="${IN_FASTA:-$CAL/molluscan_null_all.fa}"
TAG="${TAG:-esm2_3b}"
mkdir -p "$SH" logs

i="${SLURM_ARRAY_TASK_ID:-0}"
FA="$SH/${TAG}_shard${i}.fa"
OUT="$SH/${TAG}_shard${i}.npz"

python3 scripts/unity/molluscan_calibration_shard.py \
    --input "$IN" --out "$FA" --nshards "$NSHARD" --shard "$i"

if [ -s "$OUT" ]; then
    echo "[embed] shard $i already done: $OUT"
    exit 0
fi

MODE=esm \
MODEL_ID=facebook/esm2_t36_3B_UR50D \
TOKENIZER_ID=facebook/esm2_t36_3B_UR50D \
INPUT_FAA="$FA" \
OUT_NPZ="$OUT" \
python3 scratch_hf_auto_embed.py

echo "[embed] shard $i DONE $(date -Is)"
