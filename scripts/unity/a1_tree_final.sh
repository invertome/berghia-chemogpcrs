#!/bin/bash
#SBATCH --job-name=a1_tree_final
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=14-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# A1 dedicated class-A tree — FINAL IQ-TREE 3 (run after a1_tree_pilot.sh picks
# the filter condition). Builds the tree A1's tree_distance_to_refs.py consumes.
#
# ALN selects the trimmed alignment chosen by the pilot (default: canonical +
# smart-gap ClipKit, the minimal-filtering condition — preserves the ECL/ICL
# divergence A1's novelty signal depends on). Override ALN to use a
# PREQUAL/TAPER condition if the pilot verdict favours it.
#
# Model selection is restricted to nuclear AA matrices (-msub nuclear), UFBoot +
# SH-aLRT + TBE, deterministic seed — matching the pipeline's per-class trees.
set -eo pipefail

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
OUTDIR="${OUTDIR:-${REPO}/results/phylogenies/protein/class_A_a1}"
ALN="${ALN:-${OUTDIR}/canonical_trim.fa}"
PREFIX="${PREFIX:-${OUTDIR}/class_A_a1}"
THREADS="${THREADS:-32}"
SEED="${IQTREE_SEED:-20260718}"

source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
cd "$REPO"
source config.sh
source functions.sh

[ -s "$ALN" ] || { echo "[a1_final] ERROR: alignment $ALN missing (run a1_tree_pilot.sh first)"; exit 1; }
echo "[a1_final] IQ-TREE 3 on $ALN ($(grep -c '^>' "$ALN") seqs)"

# PROFILE MIXTURE MODELS. Divergent class-A GPCR paralogs have strong
# site-specific composition (TM helices vs ECL/ICL loops) that single-matrix
# models fit poorly, so ModelFinder explores C-series profile mixtures alongside
# the standard nuclear matrices and selects by BIC. -mset restricts the base
# matrices; --madd ADDS the mixture candidates, so the two compose.
#   IQTREE_MIXTURE_MODE=madd (default) | pmsf (two-stage, for intractable scale) | off
IQTREE_MIXTURE_MODE="${IQTREE_MIXTURE_MODE:-madd}"
IQTREE_MADD="${IQTREE_MADD:-C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,LG4M,LG4X}"
IQTREE_PMSF_MODEL="${IQTREE_PMSF_MODEL:-LG+C60+F+G}"

if [ "$IQTREE_MIXTURE_MODE" = "pmsf" ]; then
    echo "[a1_final] stage 1/2: guide tree ..."
    iqtree3 -s "$ALN" --prefix "${PREFIX}_guide" \
      -m MFP -msub nuclear -mset "${IQTREE_MODEL_SET:-LG,WAG,JTT,Q.pfam}" \
      -seed "$SEED" -T "$THREADS" --redo
    echo "[a1_final] stage 2/2: PMSF ${IQTREE_PMSF_MODEL} ..."
    iqtree3 -s "$ALN" --prefix "$PREFIX" \
      -m "$IQTREE_PMSF_MODEL" -ft "${PREFIX}_guide.treefile" \
      -B 1000 -alrt 1000 --tbe -seed "$SEED" -T "$THREADS" --redo
elif [ "$IQTREE_MIXTURE_MODE" = "off" ]; then
    iqtree3 -s "$ALN" --prefix "$PREFIX" \
      -m MFP -msub nuclear -mset "${IQTREE_MODEL_SET:-LG,WAG,JTT,Q.pfam}" \
      -B 1000 -alrt 1000 --tbe -seed "$SEED" -T "$THREADS" --redo
else
    echo "[a1_final] exploring mixtures: ${IQTREE_MADD}"
    iqtree3 -s "$ALN" --prefix "$PREFIX" \
      -m MFP -msub nuclear -mset "${IQTREE_MODEL_SET:-LG,WAG,JTT,Q.pfam}" \
      --madd "$IQTREE_MADD" \
      -B 1000 -alrt 1000 --tbe -seed "$SEED" -T "$THREADS" --redo
fi

grep -m1 "Best-fit model" "${PREFIX}.log" 2>/dev/null || true
echo "[a1_final] DONE -> ${PREFIX}.treefile"
echo "[a1_final] point A1 at it: export EMB_CLASSA_TREE=${PREFIX}.treefile ; RUN_EMB_RESIDUAL_NOVELTY=1"
