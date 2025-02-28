#!/bin/bash
# 01_reference_processing.sh
# Purpose: Process reference FASTA files and build HMMs if custom HMMs are not provided, preparing data for downstream searches.
# Inputs: Reference FASTA files in ${REFERENCE_DIR}, optional custom HMMs (${CONSERVED_HMM}, ${LSE_HMM})
# Outputs: Combined reference FASTA (${RESULTS_DIR}/reference_sequences/all_references.fa), HMMs (${RESULTS_DIR}/hmms/), HHsuite DB (${RESULTS_DIR}/hhdb/)
# Notes: Supports flexible use of custom HMMs or reference-derived HMMs, always builds HHsuite DB for HHblits.

#SBATCH --job-name=reference_processing
#SBATCH --output=${LOGS_DIR}/01_reference_processing_%j.out
#SBATCH --error=${LOGS_DIR}/01_reference_processing_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/reference_sequences" "${RESULTS_DIR}/hmms" "${RESULTS_DIR}/hhdb" "${LOGS_DIR}" "${SCRIPTS_DIR}"

log "Starting reference processing."

# Process conserved HMMs
if [ -n "${CONSERVED_HMM}" ] && [ -f "${CONSERVED_HMM}" ]; then
    cp "${CONSERVED_HMM}" "${RESULTS_DIR}/hmms/conserved.hmm"
    log "Using provided custom conserved HMM: ${CONSERVED_HMM}"
else
    log "No custom conserved HMM provided, building from reference sequences."
    for taxid in "${TAXA[@]}"; do
        conserved_aa="${REFERENCE_DIR}/${taxid}_conserved_refs.aa"
        check_file "$conserved_aa"
        run_command "ref_hmm_conserved_${taxid}" ${HMMBUILD} "${RESULTS_DIR}/hmms/${taxid}_conserved.hmm" "$conserved_aa"
    done
    cat "${RESULTS_DIR}/hmms/"*_conserved.hmm > "${RESULTS_DIR}/hmms/conserved.hmm"
    check_file "${RESULTS_DIR}/hmms/conserved.hmm"
fi

# Process LSE HMMs
if [ -n "${LSE_HMM}" ] && [ -f "${LSE_HMM}" ]; then
    cp "${LSE_HMM}" "${RESULTS_DIR}/hmms/lse.hmm"
    log "Using provided custom LSE HMM: ${LSE_HMM}"
else
    log "No custom LSE HMM provided, building from reference sequences."
    for taxid in "${TAXA[@]}"; do
        lse_aa="${REFERENCE_DIR}/${taxid}_lse_refs.aa"
        check_file "$lse_aa"
        run_command "ref_hmm_lse_${taxid}" ${HMMBUILD} "${RESULTS_DIR}/hmms/${taxid}_lse.hmm" "$lse_aa"
    done
    cat "${RESULTS_DIR}/hmms/"*_lse.hmm > "${RESULTS_DIR}/hmms/lse.hmm"
    check_file "${RESULTS_DIR}/hmms/lse.hmm"
fi

# Combine all reference sequences for HHsuite database
log "Combining all reference sequences."
cat "${REFERENCE_DIR}"/*_conserved_refs.aa "${REFERENCE_DIR}"/*_lse_refs.aa > "${RESULTS_DIR}/reference_sequences/all_references.fa"
check_file "${RESULTS_DIR}/reference_sequences/all_references.fa"

# Update headers with short IDs and create mapping
log "Updating reference sequence headers with short IDs."
run_command "ref_id_map" python3 "${SCRIPTS_DIR}/update_headers.py" "${RESULTS_DIR}/reference_sequences/all_references.fa" "${ID_MAP}"
mv "${RESULTS_DIR}/reference_sequences/all_references_updated.fa" "${RESULTS_DIR}/reference_sequences/all_references.fa"

# Build HHsuite database for HHblits searches
log "Building HHsuite database from combined references."
run_command "hhdb_creation" ${HHMAKE} -i "${RESULTS_DIR}/reference_sequences/all_references.fa" -o "${RESULTS_DIR}/hhdb/references.hhm" -v 1

log "Reference processing completed successfully."
