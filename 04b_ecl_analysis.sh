#!/bin/bash
# 04b_ecl_analysis.sh
# Purpose: Analyze extracellular loop (ECL) divergence vs transmembrane (TM) conservation.
# Candidates with divergent ECLs but conserved TM domains are prioritized — this pattern
# suggests ligand-binding site diversification, a hallmark of chemoreceptor evolution.
#
# Inputs:
#   - Canonical MAFFT alignment (class_A_aligned_canonical.fa) —
#     un-cleaned, retains variable ECL columns. Stage 04's filter stack
#     (PREQUAL+CLOAK+TAPER) would mask exactly the alignment-uncertain ECL
#     positions we want to study here, so this analysis must read the
#     pre-filter canonical alignment instead. Falls back to the filter-
#     stacked output with a warning if the canonical sibling is missing
#     (e.g. for legacy runs that pre-date the May 2026 filter-stack swap).
#   - DeepTMHMM/Kyte-Doolittle topology predictions (predicted_topologies.3line)
#
# Outputs:
#   - ${RESULTS_DIR}/ecl_analysis/ecl_divergence.csv
#
# Dependencies: Step 04 (phylogenetic analysis, for alignment), Step 02 (for TM predictions)
#
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=ecl_analysis
#SBATCH --output=${LOGS_DIR}/04b_ecl_analysis_%j.out
#SBATCH --error=${LOGS_DIR}/04b_ecl_analysis_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/ecl_analysis" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

log "=== Step 04b: ECL Divergence Analysis ==="

# --- Check prerequisites ---

# Step 04 completion (we need the alignment).
# Prefer the canonical sibling (un-cleaned MAFFT alignment, ECL columns
# preserved); fall back to the filter-stacked output for legacy runs.
CANONICAL="${RESULTS_DIR}/phylogenies/protein/class_A/class_A_aligned_canonical.fa"
FILTERED="${RESULTS_DIR}/phylogenies/protein/class_A/class_A_aligned.fa"
if [ -s "$CANONICAL" ]; then
    ALIGNMENT="$CANONICAL"
    log "Using canonical (un-CLOAKed) alignment for ECL analysis: $(basename "$ALIGNMENT")"
elif [ -s "$FILTERED" ]; then
    ALIGNMENT="$FILTERED"
    log --level=WARN "Canonical sibling missing; falling back to filter-stacked alignment: $(basename "$ALIGNMENT")"
    log --level=WARN "ECL columns may be CLOAK-masked, biasing the divergence ratio. Re-run stage 04 to regenerate."
else
    log "Warning: Neither canonical nor filter-stacked alignment found at $CANONICAL or $FILTERED"
    log "Skipping ECL analysis — run step 04 first to generate the alignment."
    exit 0
fi

# DeepTMHMM predictions (3-line format required for per-residue topology)
DEEPTMHMM_DIR="${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia"
TOPO_FILE="${DEEPTMHMM_DIR}/predicted_topologies.3line"

if [ ! -f "$TOPO_FILE" ]; then
    log --level=ERROR "Required topology file not found: $TOPO_FILE"
    log --level=ERROR "Stage 02 (02_chemogpcrs_identification.sh) must be run first to produce"
    log --level=ERROR "DeepTMHMM per-residue topology predictions. 04b's SBATCH allocation"
    log --level=ERROR "(2h / 8G / 1 CPU) is sized for ECL analysis only — running DeepTMHMM"
    log --level=ERROR "(ProtT5 transformer) inline here would exceed the wallclock limit and"
    log --level=ERROR "produce nothing. Re-run stage 02 to generate: $TOPO_FILE"
    exit 1
fi

TOPO_COUNT=$(grep -c '^>' "$TOPO_FILE" 2>/dev/null || true)
TOPO_COUNT=${TOPO_COUNT:-0}
log "Found topology predictions for $TOPO_COUNT sequences"

# --- Run ECL divergence analysis ---
log "Analyzing ECL divergence (alignment: $(basename "$ALIGNMENT"), min ECL length: ${MIN_ECL_LENGTH:-5})..."

run_command "ecl_divergence_analysis" \
    python3 "${SCRIPTS_DIR}/analyze_ecl.py" \
        --alignments "$ALIGNMENT" \
        --deeptmhmm "$DEEPTMHMM_DIR" \
        --min-ecl-length "${MIN_ECL_LENGTH:-5}" \
        --output "${RESULTS_DIR}/ecl_analysis/ecl_divergence.csv"

# --- Verify output ---
ECL_OUTPUT="${RESULTS_DIR}/ecl_analysis/ecl_divergence.csv"
if [ -f "$ECL_OUTPUT" ]; then
    GENE_COUNT=$(tail -n +2 "$ECL_OUTPUT" | wc -l)
    log "ECL analysis complete: $GENE_COUNT genes scored"

    if [ "$GENE_COUNT" -gt 0 ]; then
        HIGH_RATIO=$(awk -F',' 'NR>1 && $10>=2.0 {count++} END {print count+0}' "$ECL_OUTPUT")
        log "  Genes with high ECL/TM ratio (>=2.0): $HIGH_RATIO"
    fi
else
    log "Warning: ECL analysis produced no output"
fi

# --- Completion flag ---
create_checkpoint "ecl_analysis" "ECL divergence analysis complete"
log "=== Step 04b complete ==="
