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
activate_conda_env   # self-activate the project conda env before running tools
# ONE deterministic, chronologically-correct rule for which OrthoFinder
# run is authoritative (mtime of Orthogroups.tsv). Shared by stages
# 03/03b/04/05/06c/07 so they can no longer resolve different runs.
# shellcheck source=scripts/orthofinder_paths.sh
source "${SCRIPTS_DIR:-scripts}/orthofinder_paths.sh"

# Create output directory
mkdir -p "${RESULTS_DIR}/lse_classification" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependencies (step 03 creates step_completed_03.txt)
check_file "${RESULTS_DIR}/step_completed_03.txt"

log "Starting LSE classification."

# --- Configuration ---
# USE_PYTHON_LSE: Use the Python-based lse_refine.py (recommended for NCBI taxonomy lookups)
# Set to false to use the legacy bash-based classification from config.sh LSE_LEVELS
USE_PYTHON_LSE=${USE_PYTHON_LSE:-true}

# --- Get orthogroup FASTA directory ---
# OrthoFinder puts per-OG FASTA files in Orthogroup_Sequences/, not Orthogroups/.
# Resolved off the SAME authoritative run every other stage uses (see
# scripts/orthofinder_paths.sh) — the previous unsorted `find | head -1` could
# pick a different run than stages 04/05/06c/07, so an OG id classified here
# denoted a different gene set downstream while every join still succeeded.
OG_DIR=$(resolve_orthogroup_sequences_dir "${RESULTS_DIR}/orthogroups" || true)
if [ -z "$OG_DIR" ] || [ ! -d "$OG_DIR" ]; then
    log "Error: Orthogroup_Sequences directory not found"
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

    # Process each orthogroup.
    # A failing orthogroup is skipped rather than aborting the stage (existing
    # behaviour), but its stderr is kept -- it used to go to /dev/null, so a
    # total classification wipe-out (missing taxonomy cache, unparseable id_map,
    # import error) was indistinguishable from a clean run.
    lse_refine_failures=0
    lse_refine_attempts=0
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
            --exclude-refs 2>>"${LOGS_DIR}/lse_refine_${base}.log" || lse_refine_failures=$((lse_refine_failures + 1))
        lse_refine_attempts=$((lse_refine_attempts + 1))
    done

    if [ "${lse_refine_failures}" -gt 0 ]; then
        log "Warning: lse_refine.py failed for ${lse_refine_failures}/${lse_refine_attempts} orthogroups (see ${LOGS_DIR}/lse_refine_*.log)"
        if [ "${lse_refine_failures}" -eq "${lse_refine_attempts}" ]; then
            log "Warning: NO orthogroup was classified - LSE results below are empty, not negative"
        fi
    else
        log "lse_refine.py classified ${lse_refine_attempts}/${lse_refine_attempts} orthogroups (0 failures)"
    fi

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

    seq_count=$(grep -c "^>" "$output_fasta" 2>/dev/null || true)
    seq_count=${seq_count:-0}
    log "Created ${output_fasta} (${seq_count} sequences)"
done

touch "${RESULTS_DIR}/step_completed_lse_classification.txt"
log "LSE classification completed."
