#!/bin/bash
# 03b_lse_classification.sh
# Purpose: Classify orthogroups into multilevel LSE datasets using NCBI taxonomy.
#          Uses lse_refine.py for robust taxonomy-based classification.
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

# Check dependencies (step 03 creates step_completed_03.txt)
check_file "${RESULTS_DIR}/step_completed_03.txt"

log "Starting LSE classification."

# --- Configuration ---
# USE_PYTHON_LSE: Use the Python-based lse_refine.py (recommended for NCBI taxonomy lookups)
# Set to false to use the legacy bash-based classification from config.sh LSE_LEVELS
USE_PYTHON_LSE=${USE_PYTHON_LSE:-true}

# --- Get orthogroup directory ---
OG_DIR=$(find "${RESULTS_DIR}/orthogroups" -maxdepth 5 -type d -name "Orthogroups" -path "*/Results_*/*" 2>/dev/null | head -1)
if [ -z "$OG_DIR" ] || [ ! -d "$OG_DIR" ]; then
    log "Error: Orthogroups directory not found"
    exit 1
fi

# Count orthogroups
OG_COUNT=$(ls "$OG_DIR"/OG*.fa 2>/dev/null | wc -l)
log "Found ${OG_COUNT} orthogroups to classify"

if [ "$USE_PYTHON_LSE" = true ]; then
    # --- Python-based LSE classification (recommended) ---
    # Uses NCBI taxonomy for accurate lineage-based classification
    log "Using Python-based LSE classification (lse_refine.py)"

    # Clear previous classification files
    rm -f "${RESULTS_DIR}/lse_classification/lse_"*.txt

    # Optional inputs
    GENE_TREE_ARG=""
    SYNTENY_ARG=""
    ID_MAP_ARG=""

    if [ -f "${ID_MAP}" ]; then
        ID_MAP_ARG="--id-map ${ID_MAP}"
    fi

    # Check for gene trees (if phylogenetic step completed)
    if [ -d "${RESULTS_DIR}/phylogenies/protein" ]; then
        GENE_TREE_DIR="${RESULTS_DIR}/phylogenies/protein"
    fi

    # Check for synteny data
    if [ -f "${RESULTS_DIR}/synteny/synteny_ids.txt" ]; then
        SYNTENY_ARG="--synteny-ids ${RESULTS_DIR}/synteny/synteny_ids.txt"
    fi

    # Process each orthogroup
    for og in "$OG_DIR"/OG*.fa; do
        [ -e "$og" ] || continue
        base=$(basename "$og" .fa)

        # Build gene tree argument if available
        if [ -n "$GENE_TREE_DIR" ] && [ -f "${GENE_TREE_DIR}/${base}.treefile" ]; then
            GENE_TREE_ARG="--gene-tree ${GENE_TREE_DIR}/${base}.treefile"
        else
            GENE_TREE_ARG=""
        fi

        # Run lse_refine.py for this orthogroup
        # Use absolute path to ensure it works in SLURM jobs
        python3 "${SCRIPTS_DIR}/lse_refine.py" "$og" \
            ${ID_MAP_ARG} \
            ${GENE_TREE_ARG} \
            ${SYNTENY_ARG} \
            --output-dir "${RESULTS_DIR}/lse_classification" \
            --exclude-refs 2>/dev/null || true
    done

    # Count results
    for level_file in "${RESULTS_DIR}/lse_classification/lse_"*.txt; do
        [ -e "$level_file" ] || continue
        level=$(basename "$level_file" .txt | sed 's/lse_//')
        count=$(wc -l < "$level_file")
        log "LSE level '${level}': ${count} orthogroups"
    done

else
    # --- Legacy bash-based LSE classification ---
    # Uses config.sh LSE_LEVELS for simple taxid matching
    log "Using legacy bash-based LSE classification"

    # Parse LSE levels from config
    declare -A lse_taxids
    for level in "${LSE_LEVELS[@]}"; do
        level_name=$(echo "$level" | cut -d':' -f1)
        # Normalize separators: replace commas with spaces, collapse multiple spaces
        taxids=$(echo "$level" | cut -d':' -f2 | tr ',' ' ' | tr -s ' ')
        lse_taxids["$level_name"]="$taxids"
    done

    # Process each orthogroup
    for og in "$OG_DIR"/OG*.fa; do
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
fi

# --- Combine orthogroups into LSE FASTA files ---
log "Combining classified orthogroups into FASTA files..."
for level_file in "${RESULTS_DIR}/lse_classification/lse_"*.txt; do
    [ -e "$level_file" ] || continue
    level=$(basename "$level_file" .txt | sed 's/lse_//')
    output_fasta="${RESULTS_DIR}/lse_classification/lse_${level}.fa"

    # Clear existing file
    > "$output_fasta"

    while read -r base; do
        og_file="$OG_DIR/${base}.fa"
        if [ -f "$og_file" ]; then
            cat "$og_file" >> "$output_fasta"
        else
            log "Warning: Orthogroup file not found: ${base}.fa"
        fi
    done < "$level_file"

    seq_count=$(grep -c "^>" "$output_fasta" 2>/dev/null || echo 0)
    log "Created ${output_fasta} (${seq_count} sequences)"
done

touch "${RESULTS_DIR}/step_completed_lse_classification.txt"
log "LSE classification completed."
