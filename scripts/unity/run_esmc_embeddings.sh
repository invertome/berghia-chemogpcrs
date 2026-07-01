#!/bin/bash
# run_esmc_embeddings.sh — ESM-C (300M) forward pass over a candidate .faa,
# writing per-id mean-pooled embeddings to a .npz for scripts/embedding_evidence.py
# (Rank Task 5, bead berghia-chemogpcrs-875).
#
# HARD RULE this feeds: the .npz produced here is later compared ONLY against
# non-chemoreceptor family centroids (exclusion) and a broad class-A-GPCR
# centroid (recall) — see scripts/embedding_evidence.py's module docstring.
# This script itself does no ranking; it only emits embeddings.
#
# Every embedding derived from this run should be treated as leakage-uncertain:
# Berghia (and nudibranchs generally) are ~absent from ESM-C's pretraining
# corpus, so scripts/embedding_evidence.py attaches emb_leakage_flag=True to
# every candidate unconditionally.
#
# Prerequisites: the `berghia-gpcr` conda env needs the `esm` package
# (EvolutionaryScale's SDK: `pip install esm`) plus a torch build with CUDA
# support. Neither is a default dependency of the rest of the pipeline, so
# install once with:  conda activate berghia-gpcr && pip install esm
#
# Usage:
#   sbatch scripts/unity/run_esmc_embeddings.sh
#   CANDIDATE_FAA=/path/to/other_candidates.faa OUT_NPZ=/path/to/out.npz \
#       sbatch scripts/unity/run_esmc_embeddings.sh
#
#SBATCH --job-name=esmc_embed
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/esmc_embed-%j.out
#SBATCH --error=logs/esmc_embed-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr                      # needs `esm` (pip) + CUDA torch
set -u

# --- Paths + params (override via env) ---------------------------------------
REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"

CANDIDATE_FAA="${CANDIDATE_FAA:-${REPO_ROOT}/results/chemogpcrs/chemogpcrs_berghia.fa}"
OUT_NPZ="${OUT_NPZ:-${REPO_ROOT}/results/ranking/embeddings/candidates_esmc300m.npz}"
ESMC_MODEL="${ESMC_MODEL:-esmc_300m}"             # esmc_300m | esmc_600m | esmc_6b
DEVICE="${DEVICE:-cuda}"

mkdir -p "$(dirname "${OUT_NPZ}")" logs

# --- Preflight ---------------------------------------------------------------
[ -s "${CANDIDATE_FAA}" ] || { echo "ERROR: missing candidate FASTA: ${CANDIDATE_FAA}" >&2; exit 1; }

echo "[$(date '+%H:%M:%S')] esmc_embed: model=${ESMC_MODEL} device=${DEVICE}"
echo "  input : ${CANDIDATE_FAA}"
echo "  output: ${OUT_NPZ}"

# --- Forward pass over every candidate sequence -------------------------------
# Pure inline driver: reads the FASTA, embeds each sequence with ESM-C, mean-pools
# over the sequence-length axis, and writes one .npz with one array per id.
# This is the ONLY place in the Rank-Task-5 code path that imports torch/esm —
# scripts/embedding_evidence.py (the scored/tested module) never does.
CANDIDATE_FAA="${CANDIDATE_FAA}" OUT_NPZ="${OUT_NPZ}" ESMC_MODEL="${ESMC_MODEL}" DEVICE="${DEVICE}" \
python3 - <<'PYEOF'
import os
import sys

import numpy as np
import torch
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein

candidate_faa = os.environ["CANDIDATE_FAA"]
out_npz = os.environ["OUT_NPZ"]
model_name = os.environ["ESMC_MODEL"]
device = os.environ["DEVICE"]


def read_fasta(path):
    """Minimal FASTA reader: id = first whitespace-delimited token of the
    header; trailing '*' stop codons are stripped (recover_cds_from_assemblies.py
    does the same on translated CDS)."""
    seqs = {}
    seq_id = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if seq_id is not None:
                    seqs[seq_id] = "".join(chunks).rstrip("*")
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if seq_id is not None:
            seqs[seq_id] = "".join(chunks).rstrip("*")
    return seqs


sequences = read_fasta(candidate_faa)
print(f"[esmc_embed] {len(sequences)} candidate sequences to embed", file=sys.stderr)

model = ESMC.from_pretrained(model_name).to(device)
model.eval()

vectors = {}
with torch.no_grad():
    for i, (seq_id, seq) in enumerate(sequences.items(), start=1):
        protein = ESMProtein(sequence=seq)
        tensor = model.encode(protein)
        per_residue = model.forward(tensor)          # (1, L, hidden_dim)
        mean_pooled = per_residue.mean(dim=1).squeeze(0)  # (hidden_dim,)
        vectors[seq_id] = mean_pooled.float().cpu().numpy()
        if i % 50 == 0:
            print(f"[esmc_embed] {i}/{len(sequences)} embedded", file=sys.stderr)

np.savez(out_npz, **vectors)
print(f"[esmc_embed] wrote {len(vectors)} embeddings -> {out_npz}", file=sys.stderr)
PYEOF

echo "[$(date '+%H:%M:%S')] esmc_embed: DONE -> ${OUT_NPZ}"
