#!/bin/bash
# 01_reference_processing.sh
# Purpose: Process reference FASTA files and build HMMs if custom HMMs are not provided.
#          Generates an ID mapping file for reference sequences.
#          Supports both legacy format and new nath_et_al directory structure.
# Inputs: Reference FASTA files in ${REFERENCE_DIR}
#         Legacy: taxid_conserved_refs.aa, taxid_lse_refs.aa
#         New: nath_et_al/{lse,one_to_one_ortholog}/{taxonomic_group}/Species_name.faa
# Outputs: Combined reference FASTA, HMMs, ID map CSV, optional CDS files
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
mkdir -p "${RESULTS_DIR}/reference_sequences/cds" || { log "Error: Cannot create CDS directory"; exit 1; }

log "Starting reference processing."

# --- Detect reference file structure ---
NATH_ET_AL_DIR="${REFERENCE_DIR}/nath_et_al"
USE_NATH_ET_AL=false

if [ -d "${NATH_ET_AL_DIR}" ]; then
    USE_NATH_ET_AL=true
    log "Detected nath_et_al reference structure"
fi

# --- Step 1: Rename reference files with taxid if needed (nath_et_al structure) ---
if [ "${USE_NATH_ET_AL}" = true ]; then
    # Check if files already have taxid prefix
    has_taxid_prefix=$(find "${NATH_ET_AL_DIR}" -name "*.faa" | head -1 | xargs basename 2>/dev/null | grep -E "^[0-9]+_" || true)

    if [ -z "${has_taxid_prefix}" ]; then
        log "Renaming reference files to include taxid prefix..."

        if [ -f "${SCRIPTS_DIR}/rename_references_with_taxid.py" ]; then
            run_command "ref_rename_taxid" python3 "${SCRIPTS_DIR}/rename_references_with_taxid.py" \
                "${NATH_ET_AL_DIR}" \
                --output-map "${RESULTS_DIR}/reference_sequences/taxid_rename_map.csv" \
                ${LOCAL_DB_DIR:+--db-dir "${LOCAL_DB_DIR}"}
        else
            log "Warning: rename_references_with_taxid.py not found, skipping taxid rename"
        fi
    else
        log "Reference files already have taxid prefix"
    fi
fi

# --- Step 2: Fetch CDS for reference proteins if needed (nath_et_al structure) ---
if [ "${USE_NATH_ET_AL}" = true ] && [ "${FETCH_REFERENCE_CDS:-true}" = true ]; then
    # Check if CDS files already exist
    cds_count=$(find "${RESULTS_DIR}/reference_sequences/cds" -name "*_cds.fna" 2>/dev/null | wc -l)

    if [ "${cds_count}" -eq 0 ]; then
        log "Fetching CDS sequences for reference proteins..."

        if [ -f "${SCRIPTS_DIR}/fetch_reference_cds.py" ]; then
            run_command "ref_fetch_cds" python3 "${SCRIPTS_DIR}/fetch_reference_cds.py" \
                "${NATH_ET_AL_DIR}" \
                -o "${RESULTS_DIR}/reference_sequences/cds" \
                --log-failed "${RESULTS_DIR}/reference_sequences/failed_cds_accessions.csv"

            log "CDS fetch complete. Check failed_cds_accessions.csv for any missing sequences."
        else
            log "Warning: fetch_reference_cds.py not found, skipping CDS fetch"
        fi
    else
        log "CDS files already exist (${cds_count} files found)"
    fi
fi

# --- Step 3: Combine reference sequences by category ---
if [ "${USE_NATH_ET_AL}" = true ]; then
    # New structure: combine by LSE vs conserved (one_to_one_ortholog)
    log "Combining reference sequences from nath_et_al structure..."

    # Combine LSE references
    find "${NATH_ET_AL_DIR}/lse" -name "*.faa" -exec cat {} + > "${RESULTS_DIR}/reference_sequences/lse_references.fa" 2>/dev/null || true
    lse_count=$(grep -c "^>" "${RESULTS_DIR}/reference_sequences/lse_references.fa" 2>/dev/null || echo 0)
    log "Combined ${lse_count} LSE reference sequences"

    # Combine conserved (one_to_one_ortholog) references
    find "${NATH_ET_AL_DIR}/one_to_one_ortholog" -name "*.faa" -exec cat {} + > "${RESULTS_DIR}/reference_sequences/conserved_references.fa" 2>/dev/null || true
    conserved_count=$(grep -c "^>" "${RESULTS_DIR}/reference_sequences/conserved_references.fa" 2>/dev/null || echo 0)
    log "Combined ${conserved_count} conserved reference sequences"

    # Combine all references
    cat "${RESULTS_DIR}/reference_sequences/lse_references.fa" "${RESULTS_DIR}/reference_sequences/conserved_references.fa" \
        > "${RESULTS_DIR}/reference_sequences/all_references.fa" 2>/dev/null

    # Combine CDS files similarly
    if [ -d "${RESULTS_DIR}/reference_sequences/cds/lse" ]; then
        find "${RESULTS_DIR}/reference_sequences/cds/lse" -name "*.fna" -exec cat {} + \
            > "${RESULTS_DIR}/reference_sequences/cds/lse_references_cds.fna" 2>/dev/null || true
    fi
    if [ -d "${RESULTS_DIR}/reference_sequences/cds/one_to_one_ortholog" ]; then
        find "${RESULTS_DIR}/reference_sequences/cds/one_to_one_ortholog" -name "*.fna" -exec cat {} + \
            > "${RESULTS_DIR}/reference_sequences/cds/conserved_references_cds.fna" 2>/dev/null || true
    fi
    cat "${RESULTS_DIR}/reference_sequences/cds/"*.fna > "${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna" 2>/dev/null || true

else
    # Legacy structure: use taxid_*_refs.aa files
    log "Using legacy reference structure..."

    # --- Process Conserved HMMs ---
    if [ -n "${CONSERVED_HMM}" ] && [ -f "${CONSERVED_HMM}" ]; then
        cp "${CONSERVED_HMM}" "${RESULTS_DIR}/hmms/conserved.hmm" || { log "Error: Failed to copy ${CONSERVED_HMM}"; exit 1; }
        log "Using custom conserved HMM: ${CONSERVED_HMM}"
    else
        for taxid in "${TAXA[@]}"; do
            conserved_aa="${REFERENCE_DIR}/${taxid}_conserved_refs.aa"
            if [ -f "$conserved_aa" ]; then
                run_command "ref_hmm_conserved_${taxid}" ${HMMBUILD} "${RESULTS_DIR}/hmms/${taxid}_conserved.hmm" "$conserved_aa"
            fi
        done
        cat "${RESULTS_DIR}/hmms/"*_conserved.hmm > "${RESULTS_DIR}/hmms/conserved.hmm" 2>/dev/null || log "Warning: No conserved HMMs generated"
    fi

    # --- Process LSE HMMs ---
    if [ -n "${LSE_HMM}" ] && [ -f "${LSE_HMM}" ]; then
        cp "${LSE_HMM}" "${RESULTS_DIR}/hmms/lse.hmm" || { log "Error: Failed to copy ${LSE_HMM}"; exit 1; }
        log "Using custom LSE HMM: ${LSE_HMM}"
    else
        for taxid in "${TAXA[@]}"; do
            lse_aa="${REFERENCE_DIR}/${taxid}_lse_refs.aa"
            if [ -f "$lse_aa" ]; then
                run_command "ref_hmm_lse_${taxid}" ${HMMBUILD} "${RESULTS_DIR}/hmms/${taxid}_lse.hmm" "$lse_aa"
            fi
        done
        cat "${RESULTS_DIR}/hmms/"*_lse.hmm > "${RESULTS_DIR}/hmms/lse.hmm" 2>/dev/null || log "Warning: No LSE HMMs generated"
    fi

    # Combine all reference sequences (legacy)
    cat "${REFERENCE_DIR}"/*_conserved_refs.aa "${REFERENCE_DIR}"/*_lse_refs.aa \
        > "${RESULTS_DIR}/reference_sequences/all_references.fa" 2>/dev/null || { log "Error: No reference sequences found"; exit 1; }
fi

# --- Step 4: Build HMMs from combined references (nath_et_al structure) ---
if [ "${USE_NATH_ET_AL}" = true ]; then
    log "Building HMMs from reference sequences..."

    # Build conserved HMM if not provided
    if [ -n "${CONSERVED_HMM}" ] && [ -f "${CONSERVED_HMM}" ]; then
        cp "${CONSERVED_HMM}" "${RESULTS_DIR}/hmms/conserved.hmm"
        log "Using custom conserved HMM: ${CONSERVED_HMM}"
    elif [ -s "${RESULTS_DIR}/reference_sequences/conserved_references.fa" ]; then
        # Align sequences first
        run_command "ref_align_conserved" ${MAFFT} --auto "${RESULTS_DIR}/reference_sequences/conserved_references.fa" \
            > "${RESULTS_DIR}/reference_sequences/conserved_references_aligned.fa"
        run_command "ref_hmm_conserved" ${HMMBUILD} "${RESULTS_DIR}/hmms/conserved.hmm" \
            "${RESULTS_DIR}/reference_sequences/conserved_references_aligned.fa"
        log "Built conserved HMM from ${conserved_count} sequences"
    fi

    # Build LSE HMM if not provided
    if [ -n "${LSE_HMM}" ] && [ -f "${LSE_HMM}" ]; then
        cp "${LSE_HMM}" "${RESULTS_DIR}/hmms/lse.hmm"
        log "Using custom LSE HMM: ${LSE_HMM}"
    elif [ -s "${RESULTS_DIR}/reference_sequences/lse_references.fa" ]; then
        # Align sequences first
        run_command "ref_align_lse" ${MAFFT} --auto "${RESULTS_DIR}/reference_sequences/lse_references.fa" \
            > "${RESULTS_DIR}/reference_sequences/lse_references_aligned.fa"
        run_command "ref_hmm_lse" ${HMMBUILD} "${RESULTS_DIR}/hmms/lse.hmm" \
            "${RESULTS_DIR}/reference_sequences/lse_references_aligned.fa"
        log "Built LSE HMM from ${lse_count} sequences"
    fi
fi

# --- Step 5: Update headers with short IDs and generate ID map ---
# Check if we have any references
ref_count=$(grep -c "^>" "${RESULTS_DIR}/reference_sequences/all_references.fa" 2>/dev/null || echo 0)
if [ "${ref_count}" -eq 0 ]; then
    log "Error: No reference sequences found in all_references.fa"
    exit 1
fi

log "Processing ${ref_count} total reference sequences..."

run_command "ref_id_map" python3 "${SCRIPTS_DIR}/update_headers.py" \
    "${RESULTS_DIR}/reference_sequences/all_references.fa" "${ID_MAP}" \
    --source-type reference
mv "${RESULTS_DIR}/reference_sequences/all_references_updated.fa" "${RESULTS_DIR}/reference_sequences/all_references.fa" || { log "Error: Failed to move updated reference file"; exit 1; }

# --- Step 6: Update CDS headers to match protein headers ---
if [ -f "${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna" ]; then
    log "Updating CDS headers to match protein ID mapping..."
    python3 "${SCRIPTS_DIR}/update_headers.py" \
        "${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna" "${ID_MAP}" \
        --source-type reference --append
    mv "${RESULTS_DIR}/reference_sequences/cds/all_references_cds_updated.fna" \
       "${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna" 2>/dev/null || true
fi

# Create completion flag
touch "${RESULTS_DIR}/step_completed_01.txt"
log "Reference processing completed. Processed ${ref_count} sequences."
