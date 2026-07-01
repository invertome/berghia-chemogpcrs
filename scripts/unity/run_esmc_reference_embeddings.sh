#!/bin/bash
# run_esmc_reference_embeddings.sh — ESM-C (300M) forward pass over the
# curated GPCR family reference set (references/anchors/anchor_set.fasta),
# writing per-accession mean-pooled embeddings to a .npz for
# scripts/build_embedding_channel.py (Glue task G2, bead berghia-chemogpcrs-875).
#
# This is the REFERENCE-side counterpart to scripts/unity/run_esmc_embeddings.sh
# (which embeds the candidate .faa). build_embedding_channel.py's
# build_family_centroids() consumes THIS script's output .npz (--ref-npz)
# together with references/anchors/anchor_set.tsv (--ref-labels, the
# accession->family/class labels) to build the non-chemoreceptor family
# centroids + the broad class-A-GPCR centroid; scripts/embedding_evidence.py
# then compares each candidate embedding (from run_esmc_embeddings.sh)
# against those centroids. The candidate and reference .npz files MUST come
# from the SAME ESM-C checkpoint (ESMC_MODEL) -- comparing embeddings from
# different model sizes/checkpoints is not meaningful.
#
# HARD RULE this feeds: the .npz produced here is later compared ONLY as
# non-chemoreceptor family centroids (exclusion) and a broad class-A-GPCR
# centroid (recall) -- see scripts/embedding_evidence.py's module docstring.
# references/anchors/anchor_set.fasta is the curated NON-chemoreceptor GPCR
# reference set by construction (build_anchor_set.py / curate_gpcr_references.py);
# this script itself does no ranking or classification, it only emits
# embeddings.
#
# Every embedding derived from this run should be treated as leakage-uncertain:
# Berghia (and nudibranchs generally) are ~absent from ESM-C's pretraining
# corpus, so scripts/embedding_evidence.py attaches emb_leakage_flag=True to
# every candidate unconditionally -- this reference-side script does not
# change that caveat, it only supplies the centroids the candidate side is
# compared against.
#
# Prerequisites: the `berghia-gpcr` conda env needs the `esm` package
# (EvolutionaryScale's SDK: `pip install esm`) plus a torch build with CUDA
# support -- the same prerequisite as run_esmc_embeddings.sh. Install once with:
#   conda activate berghia-gpcr && pip install esm
#
# Usage:
#   sbatch scripts/unity/run_esmc_reference_embeddings.sh
#   REF_FAA=/path/to/other_reference.fasta OUT_NPZ=/path/to/reference_out.npz \
#       sbatch scripts/unity/run_esmc_reference_embeddings.sh
#
#SBATCH --job-name=esmc_embed_ref
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/esmc_embed_ref-%j.out
#SBATCH --error=logs/esmc_embed_ref-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr                      # needs `esm` (pip) + CUDA torch
set -u

# --- Paths + params (override via env) ---------------------------------------
REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"

REF_FAA="${REF_FAA:-${REPO_ROOT}/references/anchors/anchor_set.fasta}"
OUT_NPZ="${OUT_NPZ:-${REPO_ROOT}/results/ranking/embeddings/reference_esmc300m.npz}"
ESMC_MODEL="${ESMC_MODEL:-esmc_300m}"             # esmc_300m | esmc_600m | esmc_6b
DEVICE="${DEVICE:-cuda}"

mkdir -p "$(dirname "${OUT_NPZ}")" logs

# --- Preflight ---------------------------------------------------------------
[ -s "${REF_FAA}" ] || { echo "ERROR: missing reference FASTA: ${REF_FAA}" >&2; exit 1; }

echo "[$(date '+%H:%M:%S')] esmc_embed_ref: model=${ESMC_MODEL} device=${DEVICE}"
echo "  input : ${REF_FAA}"
echo "  output: ${OUT_NPZ}"

# --- Forward pass over every reference sequence -------------------------------
# Pure inline driver: reads the FASTA, embeds each sequence with ESM-C, mean-pools
# over the sequence-length axis, and writes one .npz with one array per
# accession. This is the ONLY place in the reference-embedding path that
# imports torch/esm -- scripts/build_embedding_channel.py (the tested glue)
# and scripts/embedding_evidence.py (the scored/tested module) never do.
REF_FAA="${REF_FAA}" OUT_NPZ="${OUT_NPZ}" ESMC_MODEL="${ESMC_MODEL}" DEVICE="${DEVICE}" \
python3 - <<'PYEOF'
import os
import sys

import numpy as np
import torch
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig

ref_faa = os.environ["REF_FAA"]
out_npz = os.environ["OUT_NPZ"]
model_name = os.environ["ESMC_MODEL"]
device = os.environ["DEVICE"]


def read_fasta(path):
    """Minimal FASTA reader: id = first whitespace-delimited token of the
    header (matches anchor_set.tsv's `accession` column, e.g.
    'ANCHOR_A_1_P31356'); trailing '*' stop codons are stripped
    (recover_cds_from_assemblies.py does the same on translated CDS)."""
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


sequences = read_fasta(ref_faa)
print(f"[esmc_embed_ref] {len(sequences)} reference sequences to embed", file=sys.stderr)

model = ESMC.from_pretrained(model_name).to(device)
model.eval()

vectors = {}
with torch.no_grad():
    for i, (seq_id, seq) in enumerate(sequences.items(), start=1):
        protein = ESMProtein(sequence=seq)
        tensor = model.encode(protein)               # -> ESMProteinTensor
        # ESM-C SDK: get embeddings via logits(..., LogitsConfig(return_embeddings=True)),
        # which batches encode->forward and unwraps the hidden states. Do NOT call
        # model.forward(tensor).mean(): forward() returns an ESMCOutput dataclass
        # (no .mean()) and expects a raw token tensor, not the ESMProteinTensor
        # container. The `esm` API can drift across versions — confirm this call
        # against the Unity env's installed `esm` version on the first run.
        # (Mirrors the corrected call in run_esmc_embeddings.sh verbatim, since
        # the candidate and reference sides must stay API-identical.)
        logits_output = model.logits(tensor, LogitsConfig(return_embeddings=True))
        per_residue = logits_output.embeddings       # (1, L, hidden_dim) tensor
        mean_pooled = per_residue.mean(dim=1).squeeze(0)  # (hidden_dim,)
        vectors[seq_id] = mean_pooled.float().cpu().numpy()
        if i % 50 == 0:
            print(f"[esmc_embed_ref] {i}/{len(sequences)} embedded", file=sys.stderr)

np.savez(out_npz, **vectors)
print(f"[esmc_embed_ref] wrote {len(vectors)} embeddings -> {out_npz}", file=sys.stderr)
PYEOF

echo "[$(date '+%H:%M:%S')] esmc_embed_ref: DONE -> ${OUT_NPZ}"
