#!/bin/bash
# 03_orthology_clustering.sh
# Purpose: Cluster orthologous groups with OrthoFinder using one FASTA per species.
# Inputs: GPCR FASTA files from step 02, reference sequences from step 01
# Outputs: Orthogroups in ${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/Orthogroups/
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=orthology_clustering
#SBATCH --output=${LOGS_DIR}/03_orthology_clustering_%j.out
#SBATCH --error=${LOGS_DIR}/03_orthology_clustering_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directory
mkdir -p "${RESULTS_DIR}/orthogroups/input" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_extract_berghia.txt"

log "Starting orthology clustering."

# --- Prepare OrthoFinder input: one FASTA per species ---
# Berghia
cp "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" "${RESULTS_DIR}/orthogroups/input/${BERGHIA_TAXID}_berghia.fa" || { log "Error: Failed to copy Berghia GPCR FASTA"; exit 1; }

# Additional transcriptomes
for trans in "${TRANSCRIPTOME_DIR}"/*.aa; do
    sample=$(basename "$trans" .aa)
    taxid_sample="${sample}"
    cp "${RESULTS_DIR}/chemogpcrs/chemogpcrs_${taxid_sample}.fa" "${RESULTS_DIR}/orthogroups/input/${taxid_sample}.fa" 2>/dev/null || log "Warning: No GPCR FASTA for ${taxid_sample}"
done

# Include all reference sequences as a single outgroup file
# References have been renamed to ref_TAXID_N format in step 01 (preserving taxonomy)
# Group references by their taxid prefix
for taxid in "${TAXA[@]}"; do
    # Match references that contain the taxid in their ID (format: ref_TAXID_N)
    grep -A1 "^>ref_${taxid}_" "${RESULTS_DIR}/reference_sequences/all_references.fa" | sed '/^--$/d' > "${RESULTS_DIR}/orthogroups/input/ref_${taxid}.fa" 2>/dev/null

    # If no matches, try alternate format or include all refs
    if [ ! -s "${RESULTS_DIR}/orthogroups/input/ref_${taxid}.fa" ]; then
        rm -f "${RESULTS_DIR}/orthogroups/input/ref_${taxid}.fa"
        log "Note: No references found for taxid ${taxid}"
    fi
done

# If no taxid-specific refs found, include all references as single file
if [ -z "$(ls "${RESULTS_DIR}/orthogroups/input/ref_"*.fa 2>/dev/null)" ]; then
    log "Including all references as single outgroup"
    cp "${RESULTS_DIR}/reference_sequences/all_references.fa" "${RESULTS_DIR}/orthogroups/input/references.fa"
fi

# Run OrthoFinder
# Note: -M msa removed (invalid in OrthoFinder 2.5+), using default MSA method
# -a for number of BLAST threads, -t for tree inference threads
run_command "orthofinder" ${ORTHOFINDER} -f "${RESULTS_DIR}/orthogroups/input" -t "${CPUS}" -a "${CPUS}" -I "${ORTHOFINDER_INFLATION}" -S diamond -A mafft -T fasttree

# Verify output
if [ ! -d "${RESULTS_DIR}/orthogroups/OrthoFinder/Results"* ]; then
    log "Error: OrthoFinder failed to produce results"
    exit 1
fi

log "Orthology clustering completed."
