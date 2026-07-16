#!/bin/bash
#SBATCH --job-name=fusion_bake
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fusion_bake-%j.out
#SBATCH --error=logs/fusion_bake-%j.err
#SBATCH --account=pi_pkatz_umass_edu
#
# Build fused reference embeddings for named cross-paradigm pairs, then score each
# with the UNCHANGED scratch_lofo_bakeoff.py (bead cw3.16, Module B). Scoring logic
# is never re-implemented here — this wrapper only fuses then delegates.
#
# Usage: PAIRS="protrek:saprot protrek:carp" sbatch scripts/unity/run_fusion_bakeoff.sh
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"; conda activate esmc; set -u
cd /scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
E=results/ranking/embeddings
# PAIRS supplied via env, e.g. PAIRS="protrek:saprot protrek:carp"
for pair in $PAIRS; do
  a=${pair%%:*}; b=${pair##*:}; tag=${a}X${b}
  python3 scripts/fusion_embed.py --pair "$a:$E/reference_${a}_PROD.npz" \
      "$b:$E/reference_${b}_PROD.npz" --out-prefix "$E/reference_${tag}_PROD"
  for pr in "FULL-clean:references/anchors/anchor_set_PROD.tsv" \
            "VERIFIED-clean:references/anchors/anchor_set_PROD_verified.tsv"; do
    echo "########## $tag :: ${pr%%:*} ##########"
    python3 scratch_lofo_bakeoff.py "$E/reference_${tag}_PROD.npz" "-" "${pr##*:}" \
        || echo ">>> $tag ${pr%%:*} FAILED"
  done
done
echo "=== fusion_bake DONE ==="
