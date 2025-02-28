#!/bin/bash
# 01_reference_processing.sh
# Purpose: Process reference FASTA files and build HMMs if custom HMMs are not provided.
#          Generates an ID mapping file for reference sequences.
# Inputs: Reference FASTA files in ${REFERENCE_DIR} (taxid_conserved_refs.aa, taxid_lse_refs.aa), optional custom HMMs
# Outputs: Combined reference FASTA, HMMs, ID map CSV
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

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
mkdir -p "${RESULTS_DIR}/reference_sequences" "${RESULTS_DIR}/hmms" "${LOGS_DIR}" "${SCRIPTS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

log "Starting reference processing."

# --- Process Conserved HMMs ---
if [ -n "${CONSERVED_HMM}" ] && [ -f "${CONSERVED_HMM}" ]; then
    cp "${CONSERVED_HMM}" "${RESULTS_DIR}/hmms/conserved.hmm" || { log "Error: Failed to copy ${CONSERVED_HMM}"; exit 1; }
    log "Using custom conserved HMM: ${CONSERVED_HMM}"
else
    for taxid in "${TAXA[@]}"; do
        conserved_aa="${REFERENCE_DIR}/${taxid}_conserved_refs.aa"
        check_file "$conserved_aa"
        run_command "ref_hmm_conserved_${taxid}" ${HMMBUILD} "${RESULTS_DIR}/hmms/${taxid}_conserved.hmm" "$conserved_aa"
    done
    cat "${RESULTS_DIR}/hmms/"*_conserved.hmm > "${RESULTS_DIR}/hmms/conserved.hmm" || { log "Error: Failed to combine conserved HMMs"; exit 1; }
fi

# --- Process LSE HMMs ---
if [ -n "${LSE_HMM}" ] && [ -f "${LSE_HMM}" ]; then
    cp "${LSE_HMM}" "${RESULTS_DIR}/hmms/lse.hmm" || { log "Error: Failed to copy ${LSE_HMM}"; exit 1; }
    log "Using custom LSE HMM: ${LSE_HMM}"
else
    for taxid in "${TAXA[@]}"; do
        lse_aa="${REFERENCE_DIR}/${taxid}_lse_refs.aa"
        check_file "$lse_aa"
        run_command "ref_hmm_lse_${taxid}" ${HMMBUILD} "${RESULTS_DIR}/hmms/${taxid}_lse.hmm" "$lse_aa"
    done
    cat "${RESULTS_DIR}/hmms/"*_lse.hmm > "${RESULTS_DIR}/hmms/lse.hmm" || { log "Error: Failed to combine LSE HMMs"; exit 1; }
fi

# --- Combine all reference sequences ---
cat "${REFERENCE_DIR}"/*_conserved_refs.aa "${REFERENCE_DIR}"/*_lse_refs.aa > "${RESULTS_DIR}/reference_sequences/all_references.fa" 2>/dev/null || { log "Error: No reference sequences found"; exit 1; }

# --- Update headers with short IDs and generate ID map ---
run_command "ref_id_map" python3 "${SCRIPTS_DIR}/update_headers.py" "${RESULTS_DIR}/reference_sequences/all_references.fa" "${ID_MAP}"
mv "${RESULTS_DIR}/reference_sequences/all_references_updated.fa" "${RESULTS_DIR}/reference_sequences/all_references.fa" || { log "Error: Failed to move updated reference file"; exit 1; }

log "Reference processing completed."
