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
    # Normalize separators: replace commas with spaces, collapse multiple spaces
    taxids=$(echo "$level" | cut -d':' -f2 | tr ',' ' ' | tr -s ' ')
    lse_taxids["$level_name"]="$taxids"
done

# --- Process each orthogroup ---
for og in "${RESULTS_DIR}/orthogroups/OrthoFinder/Results"*/Orthogroups/OG*.fa; do
    # Skip if glob doesn't match any files
    [ -e "$og" ] || continue

    base=$(basename "$og" .fa)

    # Extract taxids using centralized metadata lookup (excludes references)
    # This uses get_taxids_from_fasta from functions.sh which queries the metadata CSV
    # Fallback to header parsing is built into the function if CSV doesn't exist
    taxids=$(get_taxids_from_fasta "$og" --exclude-refs)

    # Check against each LSE level
    for level in "${!lse_taxids[@]}"; do
        level_taxids="${lse_taxids[$level]}"

        # Use read with IFS to safely split into arrays
        IFS=' ' read -r -a level_taxid_array <<< "$level_taxids"
        IFS=' ' read -r -a taxid_array <<< "$taxids"

        match=true
        for t in "${taxid_array[@]}"; do
            # Skip empty elements
            [ -z "$t" ] && continue
            # Check if taxid is in the level's taxid list
            if [[ ! " ${level_taxids} " =~ " ${t} " ]]; then
                match=false
                break
            fi
        done

        # Verify exact match (same number of taxids)
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
