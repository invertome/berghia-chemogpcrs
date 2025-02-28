#!/bin/bash
# 07_candidate_ranking.sh
# Purpose: Rank GPCR candidates using integrated data and expression.
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

mkdir -p "${RESULTS_DIR}/ranking" "${LOGS_DIR}"

if [ ! -f "${RESULTS_DIR}/step_completed_map_berghia_index.txt" ]; then
    log "Error: Mapping step not completed."
    exit 1
fi

log "Starting candidate ranking."

if [ ! -f "${RESULTS_DIR}/step_completed_expression.txt" ]; then
    for tissue in rhinophore oral-tentacle; do
        for quant_file in "${QUANT_DIR}/${tissue}_"*/quant.sf; do
            awk -v tpm_thresh="${TPM_THRESHOLD}" -v tissue="$tissue" -v w_r="${TPM_WEIGHT_RHINOPHORE}" -v w_o="${TPM_WEIGHT_ORAL}" \
                '$5 > tpm_thresh {print $1 "," (tissue == "rhinophore" ? w_r : w_o)}' "$quant_file" >> "${RESULTS_DIR}/ranking/expressed_ids.txt"
        done
    done
    touch "${RESULTS_DIR}/step_completed_expression.txt"
fi

run_command "ranking" python3 "${SCRIPTS_DIR}/rank_candidates.py" "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt" "${RESULTS_DIR}/ranking/expressed_ids.txt" \
    "${RESULTS_DIR}/phylogenies/protein" "${RESULTS_DIR}/selective_pressure" "${RESULTS_DIR}/synteny/gpcr_synteny_ids.txt" "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv"

python3 "${SCRIPTS_DIR}/plot_ranking.py" "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" "${RESULTS_DIR}/ranking/ranking_plot.png"

log "Candidate ranking completed."
