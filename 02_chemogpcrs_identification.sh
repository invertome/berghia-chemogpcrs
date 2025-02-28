#!/bin/bash
# 02_chemogpcrs_identification.sh
# Purpose: Identify chemoreceptive GPCRs using HMMSEARCH (if custom HMMs provided) and HHblits, filtering with DeepTMHMM.
# Inputs: Transcriptome files in ${TRANSCRIPTOME_DIR}/*.aa, reference HMMs and sequences from step 01
# Outputs: GPCR FASTA files in ${RESULTS_DIR}/chemogpcrs/chemogpcrs_*.fa
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=chemogpcrs_id
#SBATCH --output=${LOGS_DIR}/02_chemogpcrs_id_%j.out
#SBATCH --error=${LOGS_DIR}/02_chemogpcrs_id_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/chemogpcrs" "${RESULTS_DIR}/hhdb" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_ref_id_map.txt"

log "Starting chemoreceptive GPCR identification."

# --- Process Berghia transcriptome ---
run_command "deeptmhmm_berghia" ${DEEPTMHMM} -f "${TRANSCRIPTOME}" -o "${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia"
awk -v min_tm="${MIN_TM_REGIONS}" '$5 >= min_tm {print $1}' "${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia/prediction" > "${RESULTS_DIR}/chemogpcrs/complete_ids_berghia.txt" || { log "Error: Failed to parse DeepTMHMM"; exit 1; }
run_command "seqtk_complete_berghia" ${SEQTK} subseq "${TRANSCRIPTOME}" "${RESULTS_DIR}/chemogpcrs/complete_ids_berghia.txt" > "${RESULTS_DIR}/chemogpcrs/complete_berghia.fa"

# --- Build HHsuite database from references ---
run_command "hhdb_creation" ${HHMAKE} -i "${RESULTS_DIR}/reference_sequences/all_references.fa" -o "${RESULTS_DIR}/hhdb/references.hhm" -v 1

# --- HHblits search ---
run_command "hhblits_berghia" ${HHBLITS} -i "${RESULTS_DIR}/chemogpcrs/complete_berghia.fa" -d "${RESULTS_DIR}/hhdb/references.hhm" -o "${RESULTS_DIR}/chemogpcrs/hhblits_berghia.hhr" -e "${HHBLITS_EVALUE}" -cpu "${CPUS}"
awk -v evalue="${HHBLITS_EVALUE}" '$5 < evalue {print $1}' "${RESULTS_DIR}/chemogpcrs/hhblits_berghia.hhr" | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hhblits_berghia.txt" || { log "Error: Failed to extract HHblits hits"; exit 1; }

# --- HMMSEARCH if custom HMMs provided ---
if [ -f "${RESULTS_DIR}/hmms/conserved.hmm" ]; then
    run_command "hmmsearch_conserved_berghia" ${HMMSEARCH} --domtblout "${RESULTS_DIR}/chemogpcrs/hmmsearch_conserved_berghia.domtbl" -E "${HMM_EVALUE}" "${RESULTS_DIR}/hmms/conserved.hmm" "${RESULTS_DIR}/chemogpcrs/complete_berghia.fa"
    awk '$12 < '"${HMM_EVALUE}"' {print $1}' "${RESULTS_DIR}/chemogpcrs/hmmsearch_conserved_berghia.domtbl" | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hmmsearch_conserved_berghia.txt"
fi

if [ -f "${RESULTS_DIR}/hmms/lse.hmm" ]; then
    run_command "hmmsearch_lse_berghia" ${HMMSEARCH} --domtblout "${RESULTS_DIR}/chemogpcrs/hmmsearch_lse_berghia.domtbl" -E "${HMM_EVALUE}" "${RESULTS_DIR}/hmms/lse.hmm" "${RESULTS_DIR}/chemogpcrs/complete_berghia.fa"
    awk '$12 < '"${HMM_EVALUE}"' {print $1}' "${RESULTS_DIR}/chemogpcrs/hmmsearch_lse_berghia.domtbl" | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hmmsearch_lse_berghia.txt"
fi

# --- Combine IDs from all searches ---
cat "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hhblits_berghia.txt" \
    "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hmmsearch_conserved_berghia.txt" \
    "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hmmsearch_lse_berghia.txt" 2>/dev/null | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt"
[ -s "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt" ] || { log "Error: No GPCR IDs identified"; exit 1; }

# --- Extract GPCR sequences ---
run_command "extract_berghia" ${SEQTK} subseq "${RESULTS_DIR}/chemogpcrs/complete_berghia.fa" "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt" > "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"

# --- Process additional transcriptomes ---
for trans in "${TRANSCRIPTOME_DIR}"/*.aa; do
    sample=$(basename "$trans" .aa)
    taxid_sample="${sample}"
    run_command "deeptmhmm_${taxid_sample}" ${DEEPTMHMM} -f "$trans" -o "${RESULTS_DIR}/chemogpcrs/deeptmhmm_${taxid_sample}"
    awk -v min_tm="${MIN_TM_REGIONS}" '$5 >= min_tm {print $1}' "${RESULTS_DIR}/chemogpcrs/deeptmhmm_${taxid_sample}/prediction" > "${RESULTS_DIR}/chemogpcrs/complete_ids_${taxid_sample}.txt"
    run_command "seqtk_complete_${taxid_sample}" ${SEQTK} subseq "$trans" "${RESULTS_DIR}/chemogpcrs/complete_ids_${taxid_sample}.txt" > "${RESULTS_DIR}/chemogpcrs/complete_${taxid_sample}.fa"
    
    run_command "hhblits_${taxid_sample}" ${HHBLITS} -i "${RESULTS_DIR}/chemogpcrs/complete_${taxid_sample}.fa" -d "${RESULTS_DIR}/hhdb/references.hhm" -o "${RESULTS_DIR}/chemogpcrs/hhblits_${taxid_sample}.hhr" -e "${HHBLITS_EVALUE}" -cpu "${CPUS}"
    awk -v evalue="${HHBLITS_EVALUE}" '$5 < evalue {print $1}' "${RESULTS_DIR}/chemogpcrs/hhblits_${taxid_sample}.hhr" | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hhblits_${taxid_sample}.txt"
    
    if [ -f "${RESULTS_DIR}/hmms/conserved.hmm" ]; then
        run_command "hmmsearch_conserved_${taxid_sample}" ${HMMSEARCH} --domtblout "${RESULTS_DIR}/chemogpcrs/hmmsearch_conserved_${taxid_sample}.domtbl" -E "${HMM_EVALUE}" "${RESULTS_DIR}/hmms/conserved.hmm" "${RESULTS_DIR}/chemogpcrs/complete_${taxid_sample}.fa"
        awk '$12 < '"${HMM_EVALUE}"' {print $1}' "${RESULTS_DIR}/chemogpcrs/hmmsearch_conserved_${taxid_sample}.domtbl" | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hmmsearch_conserved_${taxid_sample}.txt"
    fi
    
    if [ -f "${RESULTS_DIR}/hmms/lse.hmm" ]; then
        run_command "hmmsearch_lse_${taxid_sample}" ${HMMSEARCH} --domtblout "${RESULTS_DIR}/chemogpcrs/hmmsearch_lse_${taxid_sample}.domtbl" -E "${HMM_EVALUE}" "${RESULTS_DIR}/hmms/lse.hmm" "${RESULTS_DIR}/chemogpcrs/complete_${taxid_sample}.fa"
        awk '$12 < '"${HMM_EVALUE}"' {print $1}' "${RESULTS_DIR}/chemogpcrs/hmmsearch_lse_${taxid_sample}.domtbl" | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hmmsearch_lse_${taxid_sample}.txt"
    fi
    
    cat "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hhblits_${taxid_sample}.txt" \
        "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hmmsearch_conserved_${taxid_sample}.txt" \
        "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hmmsearch_lse_${taxid_sample}.txt" 2>/dev/null | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_${taxid_sample}.txt"
    run_command "extract_${taxid_sample}" ${SEQTK} subseq "${RESULTS_DIR}/chemogpcrs/complete_${taxid_sample}.fa" "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_${taxid_sample}.txt" > "${RESULTS_DIR}/chemogpcrs/chemogpcrs_${taxid_sample}.fa"
done

log "Chemoreceptive GPCR identification completed."
