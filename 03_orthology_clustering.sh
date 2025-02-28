#!/bin/bash
# 03_orthology_clustering.sh
# Purpose: Cluster orthologous groups with OrthoFinder, optimized for GPCR diversity.
# Inputs: GPCR FASTA files from step 02 (${RESULTS_DIR}/chemogpcrs/chemogpcrs_*.fa), reference sequences from step 01
# Outputs: Orthogroups (${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/Orthogroups/)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

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
mkdir -p "${RESULTS_DIR}/orthogroups" "${LOGS_DIR}"

# Check dependency from step 02
if [ ! -f "${RESULTS_DIR}/step_completed_extract_berghia.txt" ]; then
    log "Error: Chemoreceptive GPCR identification step not completed."
    exit 1
fi

log "Starting orthology clustering."

# Combine all sequences for OrthoFinder input
cat "${RESULTS_DIR}/reference_sequences/all_references.fa" "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" "${TRANSCRIPTOME_DIR}"/*.aa > "${RESULTS_DIR}/orthogroups/input_seqs.fa"

# Run OrthoFinder to cluster orthologous groups
run_command "orthofinder" ${ORTHOFINDER} -f "${RESULTS_DIR}/orthogroups/input_seqs.fa" -t "${CPUS}" -I "${ORTHOFINDER_INFLATION}" -S diamond -M msa -A mafft -T fasttree

log "Orthology clustering completed."
