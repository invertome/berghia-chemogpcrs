#!/bin/bash
#SBATCH --job-name=a1_backbone
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=14-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# A1 reference BACKBONE build — SPECIES-INDEPENDENT, built ONCE.
#
# This produces the reusable placement backbone for A1's phylogeny-residualized
# novelty (novelty as EXCESS divergence beyond phylogenetic expectation, which
# needs each candidate's phylogenetic distance to its nearest characterized
# reference).
#
# WHY A BACKBONE INSTEAD OF A DE-NOVO TREE: the earlier A1 tree was built from
# candidates + anchors TOGETHER. That works as a one-off but does not
# generalize — every new species would force a full tree rebuild. Here the tree
# is inferred from the characterized class-A anchors ALONE, so it contains no
# species-specific information. Downstream, each species' candidates are simply
# PLACED onto this fixed backbone with EPA-ng, and the candidate->nearest-anchor
# distance is read off the placement.
#
# CONSEQUENCE: this script is run ONCE and its output is reused for every input
# species. Rebuild it ONLY when the reference/anchor set itself changes (new or
# revised characterized anchors) — NEVER per input species.
#
# INPUT: the class-A characterized anchors (anchor_set_PROD_classA.fasta).
# Candidates are DELIBERATELY excluded — including them would make the backbone
# species-dependent and defeat the entire design.
#
# ALIGNMENT: MAFFT-DASH E-INS-i with --originalseqonly. --originalseqonly is
# REQUIRED with --dash: DASH contributes structural priors as a GUIDE only, and
# its retrieved homologs must NOT be retained in the output alignment (they are
# not part of the anchor set and would pollute the backbone leaf set).
#
# Env (all have defaults): REPO, ANCHOR_FASTA, BACKBONE_DIR, THREADS, SEED,
#   RUN_BACKBONE_PREQUAL, RUN_BACKBONE_TAPER, RUN_BACKBONE_CLIPKIT.
set -eo pipefail

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
ANCHOR_FASTA="${ANCHOR_FASTA:-${REPO}/references/anchors/derived/anchor_set_PROD_classA.fasta}"
BACKBONE_DIR="${BACKBONE_DIR:-${REPO}/results/phylogenies/protein/a1_backbone}"
THREADS="${THREADS:-32}"
SEED="${SEED:-20260718}"

# --- Alignment-cleanup gates -------------------------------------------------
# These MIRROR the filter conditions being evaluated by the A1 filter-stack
# pilot (a1_tree_pilot.sh, 2x2 PREQUAL x TAPER with ClipKit held constant).
# That pilot is still in flight, so nothing is hard-coded here: PREQUAL and
# TAPER default OFF (minimal filtering, which preserves the ECL/ICL divergence
# the A1 novelty signal lives in) and ClipKit smart-gap defaults ON.
# ACTION ITEM: once the pilot's filter_comparison.json verdict lands, set these
# defaults to whatever condition it selects, so the backbone is built under the
# same filtering regime the pilot validated.
#
# CLOAK is intentionally NOT offered: it is a consensus-across-many-MSAs column
# masker that requires the K=5 alignment ensemble, so it does not apply to a
# single E-INS-i alignment.
RUN_BACKBONE_PREQUAL="${RUN_BACKBONE_PREQUAL:-0}"
RUN_BACKBONE_TAPER="${RUN_BACKBONE_TAPER:-0}"
RUN_BACKBONE_CLIPKIT="${RUN_BACKBONE_CLIPKIT:-1}"

# Activate conda BEFORE set -u (CONDA_BACKUP_* unbound-var trap on deactivate hooks)
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
cd "$REPO"
source config.sh
source functions.sh

# --- Self-gate: the anchor set is REQUIRED ----------------------------------
if [ ! -s "$ANCHOR_FASTA" ]; then
  echo "[a1_backbone] ERROR: anchor FASTA missing or empty: $ANCHOR_FASTA" >&2
  echo "[a1_backbone] ERROR: this input is required; set ANCHOR_FASTA or build the class-A anchor set first." >&2
  exit 1
fi

mkdir -p "$BACKBONE_DIR"
N_ANCHORS=$(grep -c '^>' "$ANCHOR_FASTA")
echo "[a1_backbone] anchors = ${N_ANCHORS} <- $ANCHOR_FASTA"

# Duplicate-leaf-id guard, INLINED on purpose: the deployed checkout's
# functions.sh predates assert_no_duplicate_fasta_ids, so calling it would
# break the job on the cluster.
DUPS=$(grep '^>' "$ANCHOR_FASTA" | sed 's/^>//; s/[[:space:]].*//' | sort | uniq -d)
if [ -n "$DUPS" ]; then
  echo "[a1_backbone] ERROR: duplicate anchor ids in $ANCHOR_FASTA:" >&2
  echo "$DUPS" >&2
  exit 1
fi

# --- Optional PREQUAL (pre-alignment residue masking) ------------------------
ALIGN_INPUT="$ANCHOR_FASTA"
if [ "$RUN_BACKBONE_PREQUAL" = "1" ]; then
  echo "[a1_backbone] PREQUAL (residue-level pre-alignment mask) ..."
  bash "${SCRIPTS_DIR}/run_prequal.sh" \
    --input="$ANCHOR_FASTA" \
    --output="${BACKBONE_DIR}/anchors_prequal.fa" \
    2>&1 | tee "${BACKBONE_DIR}/prequal.log"
  ALIGN_INPUT="${BACKBONE_DIR}/anchors_prequal.fa"
else
  echo "[a1_backbone] PREQUAL skipped (RUN_BACKBONE_PREQUAL=0)"
fi

# --- MAFFT-DASH E-INS-i ------------------------------------------------------
BACKBONE_ALN="${BACKBONE_DIR}/backbone_aln.fa"
echo "[a1_backbone] MAFFT-DASH E-INS-i (--originalseqonly) on $ALIGN_INPUT ..."
mafft --dash --originalseqonly --genafpair --maxiterate 1000 --thread "$THREADS" \
  "$ALIGN_INPUT" > "$BACKBONE_ALN" 2> "${BACKBONE_DIR}/mafft.log"
echo "[a1_backbone] aligned: $(grep -c '^>' "$BACKBONE_ALN") seqs -> $BACKBONE_ALN"

# --- Optional TAPER (per-sequence residue-outlier mask, post-alignment) ------
if [ "$RUN_BACKBONE_TAPER" = "1" ]; then
  echo "[a1_backbone] TAPER (per-sequence residue-outlier mask) ..."
  bash "${SCRIPTS_DIR}/run_taper.sh" \
    --input="$BACKBONE_ALN" \
    --output="${BACKBONE_DIR}/backbone_aln_taper.fa" \
    2>&1 | tee "${BACKBONE_DIR}/taper.log"
  BACKBONE_ALN="${BACKBONE_DIR}/backbone_aln_taper.fa"
else
  echo "[a1_backbone] TAPER skipped (RUN_BACKBONE_TAPER=0)"
fi

# --- Optional ClipKit (column trim) -----------------------------------------
# The TRIMMED alignment is the one shipped in the package: EPA-ng placement must
# use exactly the alignment the backbone tree was inferred from.
REF_ALN="${BACKBONE_DIR}/backbone_reference_alignment.fa"
if [ "$RUN_BACKBONE_CLIPKIT" = "1" ]; then
  echo "[a1_backbone] ClipKit smart-gap column trim ..."
  clipkit "$BACKBONE_ALN" -m smart-gap -o "$REF_ALN" 2>&1 | tee "${BACKBONE_DIR}/clipkit.log"
else
  echo "[a1_backbone] ClipKit skipped (RUN_BACKBONE_CLIPKIT=0)"
  cp "$BACKBONE_ALN" "$REF_ALN"
fi
# Sum the first record's residue lines: taking only the first line reports the
# FASTA wrap width (60), not the alignment width, in the provenance summary.
NCOL=$(awk '/^>/{if(seen) exit; seen=1; next} seen{n+=length($0)} END{print n+0}' "$REF_ALN")
echo "[a1_backbone] reference alignment: $(grep -c '^>' "$REF_ALN") seqs, ~${NCOL} cols -> $REF_ALN"

# --- Backbone tree: IQ-TREE 3 with PROFILE MIXTURE MODELS --------------------
# Divergent class-A GPCR paralogs have strong site-specific amino-acid
# composition (TM helices vs ECL/ICL loops), which single-matrix models fit
# poorly. ModelFinder therefore explores profile mixtures alongside the standard
# nuclear matrices and selects by BIC.
#   IQTREE_MIXTURE_MODE=madd (default) -> ModelFinder tests --madd candidates
#   IQTREE_MIXTURE_MODE=pmsf           -> two-stage PMSF (guide tree, then a
#        fixed profile mixture). PMSF is the field-standard way to apply C-series
#        mixtures at large taxon counts where full mixture ML is intractable
#        (Wang et al. 2018); use it if the madd run cannot finish in walltime.
#   IQTREE_MIXTURE_MODE=off            -> previous behaviour (no mixtures)
IQTREE_MIXTURE_MODE="${IQTREE_MIXTURE_MODE:-madd}"
IQTREE_MADD="${IQTREE_MADD:-C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,LG4M,LG4X}"
IQTREE_PMSF_MODEL="${IQTREE_PMSF_MODEL:-LG+C60+F+G}"

PREFIX="${BACKBONE_DIR}/a1_backbone"
if [ "$IQTREE_MIXTURE_MODE" = "pmsf" ]; then
    echo "[a1_backbone] IQ-TREE 3 stage 1/2: guide tree ..."
    iqtree3 -s "$REF_ALN" --prefix "${PREFIX}_guide" \
      -m MFP -msub nuclear -seed "$SEED" -T "$THREADS"
    echo "[a1_backbone] IQ-TREE 3 stage 2/2: PMSF ${IQTREE_PMSF_MODEL} ..."
    iqtree3 -s "$REF_ALN" --prefix "$PREFIX" \
      -m "$IQTREE_PMSF_MODEL" -ft "${PREFIX}_guide.treefile" \
      -B 1000 -alrt 1000 --tbe -seed "$SEED" -T "$THREADS"
elif [ "$IQTREE_MIXTURE_MODE" = "off" ]; then
    echo "[a1_backbone] IQ-TREE 3 (no mixtures) ..."
    iqtree3 -s "$REF_ALN" --prefix "$PREFIX" \
      -m MFP -msub nuclear -B 1000 -alrt 1000 --tbe -seed "$SEED" -T "$THREADS"
else
    echo "[a1_backbone] IQ-TREE 3 exploring mixtures: ${IQTREE_MADD}"
    iqtree3 -s "$REF_ALN" --prefix "$PREFIX" \
      -m MFP -msub nuclear --madd "$IQTREE_MADD" \
      -B 1000 -alrt 1000 --tbe -seed "$SEED" -T "$THREADS"
fi
BACKBONE_TREE="${PREFIX}.treefile"
grep -m1 "Best-fit model" "${PREFIX}.log" 2>/dev/null || true

# --- Anchor leaf ids (first whitespace token of each header) -----------------
# Consumed downstream to identify which backbone leaves are characterized
# references when computing candidate->nearest-reference distances.
ANCHOR_IDS="${BACKBONE_DIR}/backbone_anchor_ids.txt"
grep '^>' "$REF_ALN" | sed 's/^>//; s/[[:space:]].*//' > "$ANCHOR_IDS"

# --- Provenance summary ------------------------------------------------------
echo
echo "[a1_backbone] ================ BACKBONE PACKAGE ================"
echo "[a1_backbone] built from      : $ANCHOR_FASTA (${N_ANCHORS} class-A anchors)"
echo "[a1_backbone] filters         : PREQUAL=${RUN_BACKBONE_PREQUAL} TAPER=${RUN_BACKBONE_TAPER} CLIPKIT=${RUN_BACKBONE_CLIPKIT}"
echo "[a1_backbone] threads/seed    : ${THREADS} / ${SEED}"
echo "[a1_backbone] reference align : $REF_ALN"
echo "[a1_backbone] backbone tree   : $BACKBONE_TREE"
echo "[a1_backbone] anchor leaf ids : $ANCHOR_IDS ($(wc -l < "$ANCHOR_IDS") ids)"
echo "[a1_backbone] =================================================="
echo "[a1_backbone] SPECIES-INDEPENDENT: reuse this package for every input species."
echo "[a1_backbone] Rebuild ONLY when the characterized anchor set changes."
echo "[a1_backbone] Next: place a species' candidates onto it with EPA-ng."
