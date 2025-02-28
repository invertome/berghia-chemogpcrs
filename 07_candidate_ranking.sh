#!/bin/bash
# 07_candidate_ranking.sh
# Purpose: Rank GPCR candidates using integrated data (phylogeny, selection, synteny, expression).
# Inputs: GPCR IDs from step 02, expression data, phylogeny and selection results from steps 04 and 05, synteny from step 06
# Outputs: Ranked candidates CSV (${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv), ranking plot
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

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
mkdir -p "${RESULTS_DIR}/ranking" "${LOGS_DIR}"

# Check dependency from step 06
if [ ! -f "${RESULTS_DIR}/step_completed_map_berghia_index.txt" ]; then
    log "Error: Mapping step not completed."
    exit 1
fi

log "Starting candidate ranking."

# Parse expression data if not already done
if [ ! -f "${RESULTS_DIR}/step_completed_expression.txt" ]; then
    for tissue in rhinophore oral-tentacle; do
        for quant_file in "${QUANT_DIR}/${tissue}_"*/quant.sf; do
            awk -v tpm_thresh="${TPM_THRESHOLD}" -v tissue="$tissue" -v w_r="${TPM_WEIGHT_RHINOPHORE}" -v w_o="${TPM_WEIGHT_ORAL}" \
                '$5 > tpm_thresh {print $1 "," (tissue == "rhinophore" ? w_r : w_o)}' "$quant_file" >> "${RESULTS_DIR}/ranking/expressed_ids.txt"
        done
    done
    touch "${RESULTS_DIR}/step_completed_expression.txt"
fi

# Rank candidates using integrated data
run_command "ranking" python3 "${SCRIPTS_DIR}/rank_candidates.py" "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt" "${RESULTS_DIR}/ranking/expressed_ids.txt" \
    "${RESULTS_DIR}/phylogenies/protein" "${RESULTS_DIR}/selective_pressure" "${RESULTS_DIR}/synteny/gpcr_synteny_ids.txt" "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv"

# Plot ranking results
python3 "${SCRIPTS_DIR}/plot_ranking.py" "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" "${RESULTS_DIR}/ranking/ranking_plot"

log "Candidate ranking completed."
