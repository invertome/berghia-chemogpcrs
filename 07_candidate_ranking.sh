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
