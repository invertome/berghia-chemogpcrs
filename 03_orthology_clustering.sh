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
    [ -f "$trans" ] || continue
    # Skip Berghia transcriptome (already copied above)
    [ "$(realpath "$trans")" = "$(realpath "${TRANSCRIPTOME}")" ] && continue
    sample=$(basename "$trans" .aa)
    taxid_sample="${sample}"
    cp "${RESULTS_DIR}/chemogpcrs/chemogpcrs_${taxid_sample}.fa" "${RESULTS_DIR}/orthogroups/input/${taxid_sample}.fa" 2>/dev/null || log "Warning: No GPCR FASTA for ${taxid_sample}"
done

# --- Include reference sequences: one file per species ---
NATH_ET_AL_DIR="${REFERENCE_DIR}/nath_et_al"

if [ -d "${NATH_ET_AL_DIR}" ]; then
    # Nath et al. structure: include per-species FASTA files directly
    # ORTHOFINDER_REF_GROUPS controls which taxonomic groups to include.
    # Default: all groups. Set to e.g. "gastropoda" for faster preliminary runs.
    REF_GROUPS="${ORTHOFINDER_REF_GROUPS:-all}"

    ref_species_count=0
    for category_dir in "${NATH_ET_AL_DIR}/one_to_one_ortholog" "${NATH_ET_AL_DIR}/lse"; do
        [ -d "$category_dir" ] || continue
        for group_dir in "$category_dir"/*/; do
            [ -d "$group_dir" ] || continue
            group_name=$(basename "$group_dir")

            # Filter by group if specified
            if [ "$REF_GROUPS" != "all" ]; then
                echo "$REF_GROUPS" | tr ',' '\n' | grep -qx "$group_name" || continue
            fi

            for faa in "$group_dir"/*.faa; do
                [ -f "$faa" ] || continue
                species=$(basename "$faa" .faa)
                category=$(basename "$category_dir")
                dest="${RESULTS_DIR}/orthogroups/input/ref_${category}_${species}.fa"

                # Only copy if not already there (avoid LSE overwriting conserved)
                if [ ! -f "$dest" ]; then
                    cp "$faa" "$dest"
                    ref_species_count=$((ref_species_count + 1))
                fi
            done
        done
    done
    log "Included ${ref_species_count} reference species files for OrthoFinder"
else
    # Legacy structure: group references by taxid prefix
    for taxid in "${TAXA[@]}"; do
        awk -v pattern="^>ref_${taxid}_" '
            /^>/ {
                if (match($0, pattern)) { print_seq = 1 }
                else { print_seq = 0 }
            }
            print_seq { print }
        ' "${RESULTS_DIR}/reference_sequences/all_references.fa" > "${RESULTS_DIR}/orthogroups/input/ref_${taxid}.fa" 2>/dev/null

        if [ ! -s "${RESULTS_DIR}/orthogroups/input/ref_${taxid}.fa" ]; then
            rm -f "${RESULTS_DIR}/orthogroups/input/ref_${taxid}.fa"
            log "Note: No references found for taxid ${taxid}"
        fi
    done

    # Fallback: include all references as single file
    if [ -z "$(ls "${RESULTS_DIR}/orthogroups/input/ref_"*.fa 2>/dev/null)" ]; then
        log "Including all references as single outgroup"
        cp "${RESULTS_DIR}/reference_sequences/all_references.fa" "${RESULTS_DIR}/orthogroups/input/references.fa"
    fi
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
# -M msa required when specifying -A (aligner) or -T (tree builder)
# -a for number of BLAST threads, -t for tree inference threads
run_command "orthofinder" ${ORTHOFINDER} -f "${RESULTS_DIR}/orthogroups/input" -t "${CPUS}" -a "${CPUS}" -I "${ORTHOFINDER_INFLATION}" -M msa -S diamond -A mafft -T fasttree

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
