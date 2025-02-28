#!/bin/bash
# 01_reference_processing.sh
# Purpose: Process reference FASTA files and build HMMs if custom HMMs are not provided.
# Inputs: Reference FASTA files in ${REFERENCE_DIR}, optional custom HMMs in ${CONSERVED_HMM} and ${LSE_HMM}.
# Outputs: Combined reference FASTA (${RESULTS_DIR}/reference_sequences/all_references.fa), HMMs (${RESULTS_DIR}/hmms/).
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

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
mkdir -p "${RESULTS_DIR}/reference_sequences" "${RESULTS_DIR}/hmms" "${LOGS_DIR}" "${SCRIPTS_DIR}"

log "Starting reference processing."

# --- Process Conserved HMMs ---
if [ -n "${CONSERVED_HMM}" ] && [ -f "${CONSERVED_HMM}" ]; then
    # Use provided custom conserved HMM
    cp "${CONSERVED_HMM}" "${RESULTS_DIR}/hmms/conserved.hmm"
    log "Using custom conserved HMM: ${CONSERVED_HMM}"
else
    # Build conserved HMMs from reference sequences for each taxon
    for taxid in "${TAXA[@]}"; do
        conserved_aa="${REFERENCE_DIR}/${taxid}_conserved_refs.aa"
        check_file "$conserved_aa"
        run_command "ref_hmm_conserved_${taxid}" ${HMMBUILD} "${RESULTS_DIR}/hmms/${taxid}_conserved.hmm" "$conserved_aa"
    done
    # Combine individual HMMs into a single conserved HMM
    cat "${RESULTS_DIR}/hmms/"*_conserved.hmm > "${RESULTS_DIR}/hmms/conserved.hmm"
fi

# --- Process LSE HMMs ---
if [ -n "${LSE_HMM}" ] && [ -f "${LSE_HMM}" ]; then
    # Use provided custom LSE HMM
    cp "${LSE_HMM}" "${RESULTS_DIR}/hmms/lse.hmm"
    log "Using custom LSE HMM: ${LSE_HMM}"
else
    # Build LSE HMMs from reference sequences for each taxon
    for taxid in "${TAXA[@]}"; do
        lse_aa="${REFERENCE_DIR}/${taxid}_lse_refs.aa"
        check_file "$lse_aa"
        run_command "ref_hmm_lse_${taxid}" ${HMMBUILD} "${RESULTS_DIR}/hmms/${taxid}_lse.hmm" "$lse_aa"
    done
    # Combine individual HMMs into a single LSE HMM
    cat "${RESULTS_DIR}/hmms/"*_lse.hmm > "${RESULTS_DIR}/hmms/lse.hmm"
fi

# --- Combine all reference sequences ---
cat "${REFERENCE_DIR}"/*_conserved_refs.aa "${REFERENCE_DIR}"/*_lse_refs.aa > "${RESULTS_DIR}/reference_sequences/all_references.fa"

# --- Update headers with short IDs ---
run_command "ref_id_map" python3 "${SCRIPTS_DIR}/update_headers.py" "${RESULTS_DIR}/reference_sequences/all_references.fa" "${ID_MAP}"

# Move updated file to final location
mv "${RESULTS_DIR}/reference_sequences/all_references_updated.fa" "${RESULTS_DIR}/reference_sequences/all_references.fa"

log "Reference processing completed."
