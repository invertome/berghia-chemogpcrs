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

# Check dependency (step 02 creates step_completed_02.txt)
check_file "${RESULTS_DIR}/step_completed_02.txt"

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
    # Use awk to properly handle multi-line FASTA sequences (grep -A1 would truncate them)
    awk -v pattern="^>ref_${taxid}_" '
        /^>/ {
            if (match($0, pattern)) { print_seq = 1 }
            else { print_seq = 0 }
        }
        print_seq { print }
    ' "${RESULTS_DIR}/reference_sequences/all_references.fa" > "${RESULTS_DIR}/orthogroups/input/ref_${taxid}.fa" 2>/dev/null

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

# --- Pre-flight resource check ---
# Detect available resources
detect_resources

# Estimate memory requirements for OrthoFinder based on input size
log "Checking resource requirements for OrthoFinder..."
get_dataset_stats "${RESULTS_DIR}/orthogroups/input"

# Check if we have sufficient memory for the dataset
if ! check_resource_requirements "${RESULTS_DIR}/orthogroups/input" orthofinder; then
    log --level=WARN "Proceeding despite resource warning - monitor for OOM errors"
fi

# Run OrthoFinder
# Note: -M msa removed (invalid in OrthoFinder 2.5+), using default MSA method
# -a for number of BLAST threads, -t for tree inference threads
run_command "orthofinder" ${ORTHOFINDER} -f "${RESULTS_DIR}/orthogroups/input" -t "${CPUS}" -a "${CPUS}" -I "${ORTHOFINDER_INFLATION}" -S diamond -A mafft -T fasttree

# Verify output
ORTHOFINDER_RESULTS=$(find "${RESULTS_DIR}/orthogroups" -maxdepth 4 -type d -name "Results_*" 2>/dev/null | sort -r | head -1)
if [ -z "$ORTHOFINDER_RESULTS" ] || [ ! -d "$ORTHOFINDER_RESULTS" ]; then
    log "Error: OrthoFinder failed to produce results"
    exit 1
fi

# Create orthogroup manifest for downstream array jobs
log "Creating orthogroup manifest..."
create_orthogroup_manifest "${RESULTS_DIR}/orthogroups" "${RESULTS_DIR}/orthogroup_manifest.tsv"

# Validate outputs
validate_outputs --warn-only \
    "${ORTHOFINDER_RESULTS}/Orthogroups/Orthogroups.tsv" \
    "${RESULTS_DIR}/orthogroup_manifest.tsv"

# Create completion flag for downstream steps
touch "${RESULTS_DIR}/step_completed_03.txt"
log "Orthology clustering completed."
