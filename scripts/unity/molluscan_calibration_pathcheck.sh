#!/bin/bash
#SBATCH --job-name=moll_cal_pathchk
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=a100|a40|l40s|rtx8000
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# STAGE 2 of the geometry proof (GPU): full round trip.
#
# Takes a deterministic subset of the 790 class-A candidates -- sequences whose
# PRODUCTION embeddings already exist -- and pushes them through the EXACT
# pipeline the molluscan null uses, for BOTH models:
#
#   proteinclip3b : scratch_hf_auto_embed.py (MODE=esm, esm2_t36_3B_UR50D)
#                   -> scratch_proteinclip_project_param.py (PC_HEAD=36)
#   protrek       : ProTrek/scratch_protrek_embed.py (get_protein_repr)
#
# then compares against the production npz at both the vector level and the
# family-call level. If this passes, the molluscan background is provably in
# the same geometry as the candidates and references it calibrates; if it
# fails, no null is worth computing and the caller stops here.
#
# This is the check `inside_vert_prodbasis` cannot perform: that one re-derives
# the vertebrate readout from the EXISTING production npz and never exercises
# the code that produces the NEW molluscan embeddings.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate esmc
set -u

export HF_HOME=/scratch3/workspace/jperezmoreno_umass_edu-jorge/hf_cache
export HF_HUB_DISABLE_TELEMETRY=1
export TORCH_HOME=/scratch3/workspace/jperezmoreno_umass_edu-jorge/torch_hub

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
PR=/scratch3/workspace/jperezmoreno_umass_edu-jorge/proteinclip
PT=/scratch3/workspace/jperezmoreno_umass_edu-jorge/ProTrek
cd "$W"
E="$W/results/ranking/embeddings"
CAL="$W/results/ranking/diagnostics/molluscan_calibration"
OUT="$CAL/pathcheck"
CAND_FA="$W/results/chemogpcrs/chemogpcrs_berghia_classA.fa"
N="${PATHCHECK_N:-64}"
mkdir -p "$OUT" logs

# Deterministic subset: every Nth record of the production candidate FASTA, so
# the probe spans the whole file (and hence the whole length range) rather than
# clustering at the head. Reuses the null pipeline's own shard extractor.
STRIDE=$(python3 -c "
import sys
n=sum(1 for l in open('$CAND_FA') if l.startswith('>'))
print(max(1, n // $N))")
python3 scripts/unity/molluscan_calibration_shard.py \
    --input "$CAND_FA" --out "$OUT/probe.fa" --nshards "$STRIDE" --shard 0
echo "[pathcheck] probe: $(grep -c '^>' "$OUT/probe.fa") sequences (stride $STRIDE)"

# ---- proteinclip3b path -------------------------------------------------
if [ ! -s "$OUT/probe_esm2_3b.npz" ]; then
    MODE=esm \
    MODEL_ID=facebook/esm2_t36_3B_UR50D \
    TOKENIZER_ID=facebook/esm2_t36_3B_UR50D \
    INPUT_FAA="$OUT/probe.fa" \
    OUT_NPZ="$OUT/probe_esm2_3b.npz" \
    python3 scratch_hf_auto_embed.py
fi
if [ ! -s "$OUT/probe_proteinclip3b.npz" ]; then
    PROTEINCLIP_REPO="$PR" PC_HEAD=36 \
        BASE_NPZ="$OUT/probe_esm2_3b.npz" \
        OUT_NPZ="$OUT/probe_proteinclip3b.npz" \
        python3 scratch_proteinclip_project_param.py
fi

# ---- protrek path --------------------------------------------------------
if [ ! -s "$OUT/probe_protrek.npz" ]; then
    (
        set -eo pipefail
        source "$HOME/.miniconda3/etc/profile.d/conda.sh"
        conda activate protrek
        set -u
        cd "$PT"
        INPUT_FAA="$OUT/probe.fa" OUT_NPZ="$OUT/probe_protrek.npz" \
            python3 scratch_protrek_embed.py
    )
fi

cd "$W"

# ---- NEGATIVE CONTROL: the check must FAIL on raw ESM-2 space -------------
# If this comparison passed, the checker would be incapable of detecting the
# very failure it exists for, so a PASS here is itself a fatal result.
echo
echo "##### negative control: raw ESM-2 3B vs production proteinclip3b #####"
if python3 scripts/unity/molluscan_calibration_pathcheck.py \
        --produced-npz   "$OUT/probe_esm2_3b.npz" \
        --production-npz "$E/candidates_proteinclip3b_classA.npz" \
        --label          "NEGATIVE CONTROL (raw ESM-2 3B vs proteinclip3b)" \
        --out-json       "$OUT/pathcheck_negative_control.json"; then
    echo "FATAL: the negative control PASSED. The checker cannot distinguish" >&2
    echo "raw ESM-2 space from ProteinCLIP space and proves nothing." >&2
    exit 3
fi
echo "negative control correctly FAILED (as required)"

# ---- the real checks ------------------------------------------------------
echo
echo "##################### proteinclip3b round trip #####################"
python3 scripts/unity/molluscan_calibration_pathcheck.py \
    --produced-npz   "$OUT/probe_proteinclip3b.npz" \
    --production-npz "$E/candidates_proteinclip3b_classA.npz" \
    --ref-npz        "$E/reference_proteinclip3b_PROD.npz" \
    --ref-labels     "$W/references/anchors/anchor_set_PROD.tsv" \
    --label          "proteinclip3b full path" \
    --out-json       "$OUT/pathcheck_proteinclip3b.json"

echo
echo "######################### protrek round trip #######################"
python3 scripts/unity/molluscan_calibration_pathcheck.py \
    --produced-npz   "$OUT/probe_protrek.npz" \
    --production-npz "$E/candidates_protrek_classA.npz" \
    --ref-npz        "$E/reference_protrek_PROD.npz" \
    --ref-labels     "$W/references/anchors/anchor_set_PROD.tsv" \
    --label          "protrek full path" \
    --out-json       "$OUT/pathcheck_protrek.json"

echo
echo "[pathcheck] DONE $(date -Is)"
