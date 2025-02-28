#!/bin/bash
# 03b_lse_classification.sh
# Purpose: Classify orthogroups into multilevel LSE datasets based on taxonomic criteria from config.sh.
# Inputs: Orthogroups from step 03, species tree from step 03a
# Outputs: LSE FASTA files in ${RESULTS_DIR}/lse_classification/lse_*.fa
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=lse_classification
#SBATCH --output=${LOGS_DIR}/03b_lse_classification_%j.out
#SBATCH --error=${LOGS_DIR}/03b_lse_classification_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directory
mkdir -p "${RESULTS_DIR}/lse_classification" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependencies
check_file "${RESULTS_DIR}/step_completed_orthofinder.txt" "${SPECIES_TREE}"

log "Starting LSE classification."

# --- Parse LSE levels from config ---
declare -A lse_taxids
for level in "${LSE_LEVELS[@]}"; do
    level_name=$(echo "$level" | cut -d':' -f1)
    taxids=$(echo "$level" | cut -d':' -f2 | tr ',' ' ')
    lse_taxids["$level_name"]="$taxids"
done

# --- Process each orthogroup ---
for og in "${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/Orthogroups/OG"*.fa; do
    base=$(basename "$og" .fa)
    # Extract taxids from headers (before first '_')
    taxids=$(grep "^>" "$og" | sed 's/>//' | cut -d'_' -f1 | sort -u | tr '\n' ' ' | sed 's/ $//')
    
    # Check against each LSE level
    for level in "${!lse_taxids[@]}"; do
        level_taxids="${lse_taxids[$level]}"
        level_taxid_array=($level_taxids)
        taxid_array=($taxids)
        match=true
        for t in "${taxid_array[@]}"; do
            if [[ ! " ${level_taxids} " =~ " ${t} " ]]; then
                match=false
                break
            fi
        done
        if [ "$match" = true ] && [ ${#taxid_array[@]} -eq ${#level_taxid_array[@]} ]; then
            echo "$base" >> "${RESULTS_DIR}/lse_classification/lse_${level}.txt"
            break
        fi
    done
done

# --- Combine orthogroups into LSE FASTA files ---
for level in "${!lse_taxids[@]}"; do
    if [ -f "${RESULTS_DIR}/lse_classification/lse_${level}.txt" ]; then
        while read -r base; do
            cat "${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/Orthogroups/${base}.fa" >> "${RESULTS_DIR}/lse_classification/lse_${level}.fa" 2>/dev/null || log "Warning: Failed to append $base to lse_${level}.fa"
        done < "${RESULTS_DIR}/lse_classification/lse_${level}.txt"
    fi
done

touch "${RESULTS_DIR}/step_completed_lse_classification.txt"
log "LSE classification completed."
