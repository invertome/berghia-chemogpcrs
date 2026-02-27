#!/bin/bash
# 04b_ecl_analysis.sh
# Purpose: Analyze extracellular loop (ECL) divergence vs transmembrane (TM) conservation.
# Candidates with divergent ECLs but conserved TM domains are prioritized — this pattern
# suggests ligand-binding site diversification, a hallmark of chemoreceptor evolution.
#
# Inputs:
#   - Pre-trimming MAFFT alignment (all_berghia_refs_aligned.fa) — retains variable ECL columns
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

# Step 04 completion (we need the alignment)
ALIGNMENT="${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa"
if [ ! -s "$ALIGNMENT" ]; then
    log "Warning: Pre-trimming alignment not found or empty: $ALIGNMENT"
    log "Skipping ECL analysis — run step 04 first to generate the alignment."
    exit 0
fi

# DeepTMHMM predictions (3-line format required for per-residue topology)
DEEPTMHMM_DIR="${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia"
TOPO_FILE="${DEEPTMHMM_DIR}/predicted_topologies.3line"

if [ ! -f "$TOPO_FILE" ]; then
    log "No 3-line topology file found at $TOPO_FILE"
    log "Running DeepTMHMM wrapper to generate topology predictions..."

    BERGHIA_FASTA="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"
    if [ ! -s "$BERGHIA_FASTA" ]; then
        log "Warning: Berghia FASTA not found: $BERGHIA_FASTA"
        log "Skipping ECL analysis — need candidate sequences for TM prediction."
        exit 0
    fi

    run_command "deeptmhmm_topology" \
        bash "${SCRIPTS_DIR}/run_deeptmhmm.sh" -f "$BERGHIA_FASTA" -o "$DEEPTMHMM_DIR"

    if [ ! -f "$TOPO_FILE" ]; then
        log "Error: DeepTMHMM failed to produce topology predictions"
        exit 1
    fi
fi

TOPO_COUNT=$(grep -c '^>' "$TOPO_FILE" 2>/dev/null || echo 0)
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
