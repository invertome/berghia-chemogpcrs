#!/bin/bash
# 03b_lse_classification.sh
# Purpose: Classify orthogroups into multilevel LSE datasets (Aeolids, Nudibranchs, Gastropods) using BUSCO gene trees.
# Inputs: Orthogroups from step 03 (${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/Orthogroups/), BUSCO species tree and gene trees from step 03a
# Outputs: LSE FASTA files (${RESULTS_DIR}/lse_classification/lse_*.fa)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

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
mkdir -p "${RESULTS_DIR}/lse_classification" "${LOGS_DIR}"

# Check dependencies from steps 03 and 03a
if [ ! -f "${RESULTS_DIR}/step_completed_orthofinder.txt" ] || [ ! -f "${RESULTS_DIR}/busco/busco_species_tree.tre" ]; then
    log "Error: Orthology clustering or BUSCO species tree step not completed."
    exit 1
fi

log "Starting LSE classification."

# Verify species tree existence
check_file "${SPECIES_TREE}"

# Process each orthogroup for LSE classification
for og in "${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/Orthogroups/OG"*.fa; do
    base=$(basename "$og" .fa)
    # Locate corresponding BUSCO gene tree (not GPCR-specific)
    gene_tree=$(find "${RESULTS_DIR}/busco/orthofinder/Results*/Gene_Trees/" -name "OG*_tree.txt" | head -n 1)  # Use a representative BUSCO tree
    check_file "$gene_tree"
    python3 "${SCRIPTS_DIR}/lse_refine.py" "$og" "${ID_MAP}" "${SPECIES_TREE}" "$gene_tree" "${RESULTS_DIR}/synteny/gpcr_synteny_ids.txt" "${RESULTS_DIR}/lse_classification" 2>> "${LOGS_DIR}/03b_lse_classification.log"
done

# Combine orthogroups into LSE FASTA files
for level in "aeolids" "nudibranchs" "gastropods"; do
    if [ -f "${RESULTS_DIR}/lse_classification/lse_${level}.txt" ]; then
        while read -r base; do
            cat "${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/Orthogroups/${base}.fa" >> "${RESULTS_DIR}/lse_classification/lse_${level}.fa"
        done < "${RESULTS_DIR}/lse_classification/lse_${level}.txt"
    fi
done

touch "${RESULTS_DIR}/step_completed_lse_classification.txt"
log "LSE classification completed."
