#!/bin/bash
#SBATCH --job-name=moll_cal_protrek
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
# ProTrek_650M embeddings for the molluscan background null, via the SAME
# extractor the production protrek channel was built from
# (ProTrek/scratch_protrek_embed.py, get_protein_repr).
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate protrek
set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
PT=/scratch3/workspace/jperezmoreno_umass_edu-jorge/ProTrek
CAL="$W/results/ranking/diagnostics/molluscan_calibration"
SH="$CAL/shards"
NSHARD="${NSHARD:-8}"
IN="${IN_FASTA:-$CAL/molluscan_null_all.fa}"
mkdir -p "$SH" "$W/logs"

i="${SLURM_ARRAY_TASK_ID:-0}"
FA="$SH/protrek_shard${i}.fa"
OUT="$SH/protrek_shard${i}.npz"

python3 "$W/scripts/unity/molluscan_calibration_shard.py" \
    --input "$IN" --out "$FA" --nshards "$NSHARD" --shard "$i"

if [ -s "$OUT" ]; then
    echo "[embed] shard $i already done: $OUT"
    exit 0
fi

cd "$PT"
INPUT_FAA="$FA" OUT_NPZ="$OUT" python3 scratch_protrek_embed.py

echo "[embed] shard $i DONE $(date -Is)"
