#!/bin/bash
# 03_orthology_clustering.sh
# Purpose: Cluster orthologous groups with OrthoFinder using one FASTA per species.
# Inputs: Berghia chemoreceptor candidates (stage 02), scan-derived candidates from
#         P5/P6 (${SCAN_DERIVED_CANDIDATES_DIR}/*.chemo_candidates.fa), outgroup.
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

# Berghia chemoreceptor candidates (stage 02 output — unchanged)
cp "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" "${RESULTS_DIR}/orthogroups/input/${BERGHIA_TAXID}_berghia.fa" || { log "Error: Failed to copy Berghia GPCR FASTA"; exit 1; }

# --- Scan-derived candidates (P5 / P6) ---
# SCAN_DERIVED_CANDIDATES_DIR must exist and contain at least one
# *.chemo_candidates.fa file produced by scan_proteome_for_chemoreceptors.sh.
# Default: ${RESULTS_DIR}/p5_phase1a_validation/scan  (P5 Phase 1a output).
# Override with SCAN_DERIVED_CANDIDATES_DIR env var to point at the P6
# full-557 scan directory when that run is available.
SCAN_DIR="${SCAN_DERIVED_CANDIDATES_DIR:-${RESULTS_DIR}/p5_phase1a_validation/scan}"

if [ ! -d "${SCAN_DIR}" ]; then
    log "Error: SCAN_DERIVED_CANDIDATES_DIR does not exist: ${SCAN_DIR}"
    log "Run P5 (sbatch sbatch_run_p5_phase1a_scan.sh + sbatch_run_p5_classify_and_pool.sh) or P6 (full 557 scan) first to produce scan-derived candidates."
    exit 1
fi

shopt -s nullglob
scan_files=( "${SCAN_DIR}"/*.chemo_candidates.fa )
shopt -u nullglob

if [ "${#scan_files[@]}" -eq 0 ]; then
    log "Error: No *.chemo_candidates.fa files found in ${SCAN_DIR}"
    log "Run P5 (sbatch sbatch_run_p5_phase1a_scan.sh + sbatch_run_p5_classify_and_pool.sh) or P6 (full 557 scan) first to produce scan-derived candidates."
    exit 1
fi

scan_species_count=0
for fa in "${scan_files[@]}"; do
    # Derive a safe destination name from the filename stem.
    # Files are named <taxid>_<sanitized_binomial>.chemo_candidates.fa.
    # Strip .chemo_candidates.fa suffix; sanitize any remaining non-[A-Za-z0-9_]
    # characters (same logic as sanitize_sample_name in build_braker4_samples_csv.py).
    stem=$(basename "$fa" .chemo_candidates.fa)
    safe_name=$(printf '%s' "$stem" | tr -c 'A-Za-z0-9_' '_' | sed 's/_\+/_/g; s/^_//; s/_$//')
    dest="${RESULTS_DIR}/orthogroups/input/${safe_name}.fa"
    cp "$fa" "$dest" || { log "Warning: Failed to copy scan candidate ${fa}"; continue; }
    scan_species_count=$((scan_species_count + 1))
done
log "Included ${scan_species_count} scan-derived species files for OrthoFinder (from ${SCAN_DIR})"

# --- Outgroup (optional) ---
if [ -n "${OUTGROUP_FASTA:-}" ] && [ -f "${OUTGROUP_FASTA}" ]; then
    cp "${OUTGROUP_FASTA}" "${RESULTS_DIR}/orthogroups/input/outgroup.fa" || log "Warning: Failed to copy outgroup FASTA"
    log "Included outgroup: ${OUTGROUP_FASTA}"
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
#
# Resume logic: if a previous run exists (e.g., from a crash), use -fg to
# skip the expensive DIAMOND all-vs-all phase and resume from orthogroups.
PREV_RESULTS=$(find "${RESULTS_DIR}/orthogroups/input/OrthoFinder" -maxdepth 1 -type d -name "Results_*" 2>/dev/null | sort -r | head -1)

if [ -n "$PREV_RESULTS" ] && [ -f "$PREV_RESULTS/Orthogroups/Orthogroups.tsv" ]; then
    # Previous run has orthogroups — resume from there (skip DIAMOND + MCL)
    log "Resuming OrthoFinder from previous orthogroups: ${PREV_RESULTS}"
    run_command "orthofinder" ${ORTHOFINDER} -fg "${PREV_RESULTS}" -t "${CPUS}" -a "${CPUS}" -M msa -A mafft -T fasttree
else
    # Fresh run
    run_command "orthofinder" ${ORTHOFINDER} -f "${RESULTS_DIR}/orthogroups/input" -t "${CPUS}" -a "${CPUS}" -I "${ORTHOFINDER_INFLATION}" -M msa -S diamond -A mafft -T fasttree
fi

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
