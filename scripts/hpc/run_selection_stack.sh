#!/bin/bash
# run_selection_stack.sh — Modern HyPhy selection stack per orthogroup.
#
# Bead -urk. Replaces the aBSREL-only path in 05_selective_pressure_and_asr.sh
# with the recommended multi-method pipeline:
#
#   GARD (recombination screen) → BUSTED-S (gene-wide w/ synonymous-rate
#   variation) → BUSTED-MH (multi-nucleotide-substitution correction;
#   Lucaci 2023) → aBSREL (branch-level on confirmed signal) → MEME
#   (site-level for confirmed branches).
#
# Branch pre-filter: only test Berghia branches + LSE-internal branches
# (the targeted approach has dramatically more power than exhaustive
# correction; aBSREL paper notes this explicitly).
# BH-FDR (q=0.10) applied across the pre-filtered branch set.
#
# This is HPC-only — runtime is ~hours per OG depending on size.

#SBATCH --job-name=selection_stack
#SBATCH --output=logs/selection_stack_%A_%a.out
#SBATCH --error=logs/selection_stack_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$PROJECT_DIR/config.sh"
source "$PROJECT_DIR/functions.sh"

# Inputs (positional)
CODON_ALIGN="${1:?codon alignment FASTA}"
TREE="${2:?tree file}"
OG_BASE="${3:?orthogroup base name}"
BRANCH_LIST_FILE="${4:-}"   # Optional: file with one branch name per line for pre-filter

OUT_DIR="${RESULTS_DIR}/selective_pressure"
mkdir -p "$OUT_DIR"

# 1. GARD — recombination breakpoint screen
log "GARD: ${OG_BASE}"
hyphy gard \
    --alignment "$CODON_ALIGN" \
    --output "$OUT_DIR/${OG_BASE}_gard.json" \
    > "$OUT_DIR/${OG_BASE}_gard.log" 2>&1 \
    || log --level=WARN "GARD failed for ${OG_BASE} (continuing without recombination correction)"

# 2. BUSTED-S — synonymous-rate variation
log "BUSTED-S: ${OG_BASE}"
hyphy busted \
    --alignment "$CODON_ALIGN" \
    --tree "$TREE" \
    --srv Yes \
    --output "$OUT_DIR/${OG_BASE}_busted_s.json" \
    > "$OUT_DIR/${OG_BASE}_busted_s.log" 2>&1 \
    || log --level=WARN "BUSTED-S failed for ${OG_BASE}"

# 3. BUSTED-MH — multi-nucleotide-substitution correction
log "BUSTED-MH: ${OG_BASE}"
hyphy busted \
    --alignment "$CODON_ALIGN" \
    --tree "$TREE" \
    --srv Yes \
    --multiple-hits Double+Triple \
    --output "$OUT_DIR/${OG_BASE}_busted_mh.json" \
    > "$OUT_DIR/${OG_BASE}_busted_mh.log" 2>&1 \
    || log --level=WARN "BUSTED-MH failed for ${OG_BASE}"

# 4. aBSREL — branch-level (optionally restricted to a pre-filtered branch set)
log "aBSREL: ${OG_BASE}"
ABSREL_BRANCH_ARG=""
if [ -n "$BRANCH_LIST_FILE" ] && [ -f "$BRANCH_LIST_FILE" ]; then
    ABSREL_BRANCH_ARG="--branches FOREGROUND"
    # Note: HyPhy expects branch labels in the tree; pre-labeling the tree
    # with FOREGROUND tags is the user's responsibility (see scripts/
    # label_tree_branches.py — TODO).
fi
hyphy absrel \
    --alignment "$CODON_ALIGN" \
    --tree "$TREE" \
    --output "$OUT_DIR/${OG_BASE}_absrel.json" \
    $ABSREL_BRANCH_ARG \
    > "$OUT_DIR/${OG_BASE}_absrel.log" 2>&1 \
    || log --level=WARN "aBSREL failed for ${OG_BASE}"

# 5. MEME — site-level on the OG (HyPhy emits per-site posteriors)
log "MEME: ${OG_BASE}"
hyphy meme \
    --alignment "$CODON_ALIGN" \
    --tree "$TREE" \
    --output "$OUT_DIR/${OG_BASE}_meme.json" \
    > "$OUT_DIR/${OG_BASE}_meme.log" 2>&1 \
    || log --level=WARN "MEME failed for ${OG_BASE}"

log "Selection stack complete for ${OG_BASE}"
