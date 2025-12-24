#!/bin/bash
# 07_candidate_ranking.sh
# Purpose: Rank GPCR candidates based on phylogeny, selection, expression, and synteny data.
# Inputs: GPCR IDs from step 02, expression data, results from steps 04-06
# Outputs: Ranked candidates CSV in ${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv, plots
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=candidate_ranking
#SBATCH --output=${LOGS_DIR}/07_candidate_ranking_%j.out
#SBATCH --error=${LOGS_DIR}/07_candidate_ranking_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directory
mkdir -p "${RESULTS_DIR}/ranking" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check required dependencies
check_file "${RESULTS_DIR}/step_completed_extract_berghia.txt" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile"

# Check optional dependencies (warn but don't fail)
check_file --warn-only "${RESULTS_DIR}/selective_pressure/absrel_results.csv"
if [ ! -f "${RESULTS_DIR}/selective_pressure/absrel_results.csv" ]; then
    log "Note: aBSREL results not found - dN/dS scoring will be skipped"
fi

log "Starting candidate ranking."

# Create expression output directory
mkdir -p "${RESULTS_DIR}/expression" "${RESULTS_DIR}/gproteins" "${RESULTS_DIR}/ecl_analysis"

# --- Phase 1: Process expression data if available ---
if [ -d "${SALMON_QUANT_DIR}" ]; then
    log "Processing expression data from ${SALMON_QUANT_DIR}..."
    python3 "${SCRIPTS_DIR}/process_expression.py" \
        --quant-dir "${SALMON_QUANT_DIR}" \
        --chemosensory-tissues "${CHEMOSENSORY_TISSUES}" \
        --other-tissues "${NON_CHEMOSENSORY_TISSUES}" \
        --min-tpm "${MIN_TPM_THRESHOLD}" \
        --tau-threshold "${TAU_THRESHOLD}" \
        --output "${RESULTS_DIR}/expression/expression_summary.csv" \
        || log "Warning: Expression processing failed or no data available"
else
    log "Note: No expression data directory found at ${SALMON_QUANT_DIR}"
fi

# --- Phase 2: Classify G-proteins if references available ---
if [ -f "${GPROTEIN_REF_FASTA}" ]; then
    log "Classifying G-proteins in transcriptomes..."
    python3 "${SCRIPTS_DIR}/classify_gproteins.py" \
        --transcriptomes "${TRANSCRIPTOME_DIR}/*.aa" \
        --reference "${GPROTEIN_REF_FASTA}" \
        --classes "${GPROTEIN_REF_CLASSES}" \
        --mollusc-reference "${GPROTEIN_MOLLUSC_FASTA}" \
        --mollusc-classes "${GPROTEIN_MOLLUSC_CLASSES}" \
        --evalue "${GPROTEIN_EVALUE}" \
        --output "${RESULTS_DIR}/gproteins/classified_gproteins.csv" \
        || log "Warning: G-protein classification failed or no references available"

    # Calculate co-expression if expression data available
    if [ -f "${RESULTS_DIR}/expression/expression_summary.csv" ]; then
        log "Analyzing G-protein co-expression..."
        python3 "${SCRIPTS_DIR}/coexpression_analysis.py" \
            --expression "${RESULTS_DIR}/expression/expression_summary.csv" \
            --gproteins "${RESULTS_DIR}/gproteins/classified_gproteins.csv" \
            --candidates "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" \
            --chemosensory-tissues "${CHEMOSENSORY_TISSUES}" \
            --output "${RESULTS_DIR}/gproteins/gprotein_coexpression.csv" \
            || log "Warning: G-protein co-expression analysis failed"
    fi
else
    log "Note: No G-protein reference found at ${GPROTEIN_REF_FASTA}"
fi

# --- Phase 3: Analyze ECL divergence ---
if [ -d "${RESULTS_DIR}/deeptmhmm" ] && [ -d "${RESULTS_DIR}/phylogenies" ]; then
    log "Analyzing extracellular loop divergence..."
    python3 "${SCRIPTS_DIR}/analyze_ecl.py" \
        --alignments "${RESULTS_DIR}/phylogenies/protein/aligned_*.fasta" \
        --deeptmhmm "${RESULTS_DIR}/deeptmhmm/" \
        --min-ecl-length "${MIN_ECL_LENGTH}" \
        --output "${RESULTS_DIR}/ecl_analysis/ecl_divergence.csv" \
        || log "Warning: ECL divergence analysis failed"
else
    log "Note: DeepTMHMM or phylogeny results not found for ECL analysis"
fi

# Extract all GPCR IDs
awk '/^>/ {print substr($1,2)}' "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" > "${RESULTS_DIR}/ranking/candidate_ids.txt"

# Run ranking script with improved algorithm
# Note: Now takes synteny directory (not file) for quantitative scoring
run_command "rank_candidates" python3 "${SCRIPTS_DIR}/rank_candidates.py" \
    "${RESULTS_DIR}/ranking/candidate_ids.txt" \
    "${EXPRESSION_DATA}" \
    "${RESULTS_DIR}/phylogenies/protein" \
    "${RESULTS_DIR}/selective_pressure" \
    "${RESULTS_DIR}/synteny" \
    "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv"

# Generate plots
python3 "${SCRIPTS_DIR}/plot_ranking.py" "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" "${RESULTS_DIR}/ranking/ranking_plot" || log "Warning: Ranking plot failed"

# Create completion flag
touch "${RESULTS_DIR}/step_completed_07.txt"
log "Candidate ranking completed."
