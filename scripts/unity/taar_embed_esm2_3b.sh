#!/bin/bash
#SBATCH --job-name=taar_esm2_3b
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=a100|a40|l40s|rtx8000
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --account=pi_pkatz_umass_edu
#
# TAAR embedding job A (bead h0y0): ESM2-3B embed of the combined 2395-seq FASTA,
# then the ProteinCLIP projection (PC_HEAD=36) -> proteinclip3b_all.npz. This is
# the proteinclip3b arm; runs concurrently with the protrek arm (job B).
# Reuses the EXACT production extractors so the vectors match the existing
# reference_proteinclip3b_PROD.npz manifold.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"; conda activate esmc; set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
PC=/scratch3/workspace/jperezmoreno_umass_edu-jorge/proteinclip
export HF_HOME=/scratch3/workspace/jperezmoreno_umass_edu-jorge/hf_cache HF_HUB_DISABLE_TELEMETRY=1
cd "$W"
T=results/ranking/taar_embedding
COMB="$T/combined.fa"
mkdir -p "$T" logs
[ -s "$COMB" ] || { echo "FATAL: missing $COMB (run build_taar_embedding_fasta.py first)" >&2; exit 2; }
N=$(grep -c '^>' "$COMB"); [ "$N" -eq 2395 ] || { echo "FATAL: combined.fa has $N seqs, expected 2395" >&2; exit 2; }

# ensure the proteinclip projector head is importable (idempotent; as run_proteinclip*.sh does)
( cd "$PC" && pip install -q -e . --no-deps 2>/dev/null; pip install -q onnxruntime 2>/dev/null || true )

echo "[A] ESM2-3B embed $(date -Is)"
ESM_OUT="$T/esm2_3b_all.npz"
[ -s "$ESM_OUT" ] || MODE=esm MODEL_ID=facebook/esm2_t36_3B_UR50D TOKENIZER_ID=facebook/esm2_t36_3B_UR50D \
    INPUT_FAA="$COMB" OUT_NPZ="$ESM_OUT" python3 scratch_hf_auto_embed.py

echo "[A] ProteinCLIP projection (PC_HEAD=36) $(date -Is)"
PC_OUT="$T/proteinclip3b_all.npz"
[ -s "$PC_OUT" ] || PROTEINCLIP_REPO="$PC" PC_HEAD=36 BASE_NPZ="$ESM_OUT" OUT_NPZ="$PC_OUT" \
    python3 scratch_proteinclip_project_param.py

# verify coverage: every one of the 2395 ids embedded, no drops
python3 - "$COMB" "$PC_OUT" <<'PY'
import sys, numpy as np
ids = {l[1:].split()[0] for l in open(sys.argv[1]) if l.startswith(">")}
got = set(np.load(sys.argv[2]).files)
miss = ids - got
print(f"[A] proteinclip3b coverage {len(got & ids)}/{len(ids)}")
if miss: sys.exit(f"FATAL: {len(miss)} ids missing from proteinclip3b_all.npz e.g. {sorted(miss)[:5]}")
PY
echo "[A] DONE $(date -Is)"
