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
# ONE deterministic, chronologically-correct rule for which OrthoFinder
# run is authoritative (mtime of Orthogroups.tsv). Shared by stages
# 03/03b/04/05/06c/07 so they can no longer resolve different runs.
# shellcheck source=scripts/orthofinder_paths.sh
source "${SCRIPTS_DIR:-scripts}/orthofinder_paths.sh"

# Create output directory
mkdir -p "${RESULTS_DIR}/orthogroups/input" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency (step 02 creates step_completed_02.txt)
check_file "${RESULTS_DIR}/step_completed_02.txt"

log "Starting orthology clustering."

# --- Prepare OrthoFinder input: one FASTA per species ---

# Berghia chemoreceptor candidates: reconciled genome-track set when
# RUN_GENOME_TRACK=1 and stage 02c has produced it, else the stage-02
# transcriptome candidate set (byte-identical to legacy). Source is
# toggle-driven via berghia_candidate_fasta (functions.sh).
BERGHIA_CANDS="$(berghia_candidate_fasta)"
# Anti-silent-degradation guard: if the genome track is ENABLED but stage 02c
# never wrote reconciled_candidates.faa, berghia_candidate_fasta fell back to the
# transcriptome-only set — an operator who skipped 02c would otherwise silently
# lose the entire genome-track candidate contribution. Warn loudly; the fallback
# behavior itself is unchanged (we still proceed with the legacy set).
if [ "${RUN_GENOME_TRACK:-1}" != "0" ] && [ ! -f "${RESULTS_DIR}/reconciliation/reconciled_candidates.faa" ]; then
    log --level=WARN "RUN_GENOME_TRACK is enabled but ${RESULTS_DIR}/reconciliation/reconciled_candidates.faa is missing; falling back to the transcriptome-only candidate set. Stage 02c was likely skipped or failed — run it (sbatch 02c_genome_reconcile.sh) first to include genome-track candidates."
fi
cp "${BERGHIA_CANDS}" "${RESULTS_DIR}/orthogroups/input/${BERGHIA_TAXID}_berghia.fa" || { log "Error: Failed to copy Berghia GPCR FASTA"; exit 1; }

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
# `sort -r` on OrthoFinder's Results_<Mon><DD> stamp is neither chronological
# nor numeric ("Results_Sep05" sorts above "Results_Dec01"; "_10" below "_2"),
# and the stamp carries no year at all — so resume could latch onto a stale run.
# resolve_orthofinder_run picks the newest COMPLETED run by Orthogroups.tsv
# mtime (see scripts/orthofinder_paths.sh).
PREV_RESULTS=$(resolve_orthofinder_run "${RESULTS_DIR}/orthogroups/input/OrthoFinder" || true)

if [ -n "$PREV_RESULTS" ] && [ -f "$PREV_RESULTS/Orthogroups/Orthogroups.tsv" ]; then
    # Previous run has orthogroups — resume from there (skip DIAMOND + MCL)
    log "Resuming OrthoFinder from previous orthogroups: ${PREV_RESULTS}"
    run_command "orthofinder" ${ORTHOFINDER} -fg "${PREV_RESULTS}" -t "${CPUS}" -a "${CPUS}" -M msa -A mafft -T fasttree
else
    # Fresh run
    run_command "orthofinder" ${ORTHOFINDER} -f "${RESULTS_DIR}/orthogroups/input" -t "${CPUS}" -a "${CPUS}" -I "${ORTHOFINDER_INFLATION}" -M msa -S diamond -A mafft -T fasttree
fi

# Verify output
# Same rule as the resume probe above, so the run this stage VERIFIES is
# provably the run every downstream stage will read.
ORTHOFINDER_RESULTS=$(resolve_orthofinder_run "${RESULTS_DIR}/orthogroups" || true)
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
