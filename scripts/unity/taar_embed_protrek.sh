#!/bin/bash
#SBATCH --job-name=taar_protrek
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=a100|a40|l40s|rtx8000
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --account=pi_pkatz_umass_edu
#
# TAAR embedding job B (bead h0y0): ProTrek_650M embed of the combined 2395-seq
# FASTA -> protrek_all.npz. Runs concurrently with the proteinclip3b arm (job A).
# Reuses ProTrek/scratch_protrek_embed.py, the exact production extractor.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"; conda activate protrek; set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
PT=/scratch3/workspace/jperezmoreno_umass_edu-jorge/ProTrek
T="$W/results/ranking/taar_embedding"
COMB="$T/combined.fa"
OUT="$T/protrek_all.npz"
mkdir -p "$T" "$W/logs"
[ -s "$COMB" ] || { echo "FATAL: missing $COMB (run build_taar_embedding_fasta.py first)" >&2; exit 2; }
N=$(grep -c '^>' "$COMB"); [ "$N" -eq 2395 ] || { echo "FATAL: combined.fa has $N seqs, expected 2395" >&2; exit 2; }

echo "[B] ProTrek embed $(date -Is)"
cd "$PT"
[ -s "$OUT" ] || INPUT_FAA="$COMB" OUT_NPZ="$OUT" python3 scratch_protrek_embed.py

# verify coverage
python3 - "$COMB" "$OUT" <<'PY'
import sys, numpy as np
ids = {l[1:].split()[0] for l in open(sys.argv[1]) if l.startswith(">")}
got = set(np.load(sys.argv[2]).files)
miss = ids - got
print(f"[B] protrek coverage {len(got & ids)}/{len(ids)}")
if miss: sys.exit(f"FATAL: {len(miss)} ids missing from protrek_all.npz e.g. {sorted(miss)[:5]}")
PY
echo "[B] DONE $(date -Is)"
