#!/bin/bash
# 02_chemogpcrs_identification.sh
# Purpose: Identify chemoreceptive GPCRs using HMMs and HHblits, filtering for near-complete sequences with TM regions.
# Inputs: Transcriptome files (${TRANSCRIPTOME}, ${TRANSCRIPTOME_DIR}/*.aa), reference HMMs from 01_reference_processing.sh
# Outputs: GPCR FASTA files (${RESULTS_DIR}/chemogpcrs/chemogpcrs_*.fa)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

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
mkdir -p "${RESULTS_DIR}/chemogpcrs" "${RESULTS_DIR}/hhdb" "${LOGS_DIR}"

# Check dependency from previous step
if [ ! -f "${RESULTS_DIR}/step_completed_ref_id_map.txt" ]; then
    log "Error: Reference processing step not completed."
    exit 1
fi

log "Starting chemoreceptive GPCR identification."

# --- Filter Berghia transcriptome for near-complete sequences (≥6 TM regions) ---
run_command "deeptmhmm_berghia" ${DEEPTMHMM} -f "${TRANSCRIPTOME}" -o "${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia"
awk -v min_tm="${MIN_TM_REGIONS}" '$5 >= min_tm {print $1}' "${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia/prediction" > "${RESULTS_DIR}/chemogpcrs/complete_ids_berghia.txt"
run_command "seqtk_complete_berghia" ${SEQTK} subseq "${TRANSCRIPTOME}" "${RESULTS_DIR}/chemogpcrs/complete_ids_berghia.txt" > "${RESULTS_DIR}/chemogpcrs/complete_berghia.fa"

# --- Build HHsuite database from reference sequences ---
run_command "hhdb_creation" ${HHMAKE} -i "${RESULTS_DIR}/reference_sequences/all_references.fa" -o "${RESULTS_DIR}/hhdb/references.hhm" -v 1

# --- Search Berghia transcriptome for GPCRs using HHblits ---
run_command "hhblits_berghia" ${HHBLITS} -i "${RESULTS_DIR}/chemogpcrs/complete_berghia.fa" -d "${RESULTS_DIR}/hhdb/references.hhm" -o "${RESULTS_DIR}/chemogpcrs/hhblits_berghia.hhr" -e "${HHBLITS_EVALUE}" -cpu "${CPUS}"

# --- Extract GPCR IDs based on HHblits hits ---
awk -v evalue="${HHBLITS_EVALUE}" '$5 < evalue {print $1}' "${RESULTS_DIR}/chemogpcrs/hhblits_berghia.hhr" | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt"

# --- Extract GPCR sequences ---
run_command "extract_berghia" ${SEQTK} subseq "${RESULTS_DIR}/chemogpcrs/complete_berghia.fa" "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt" > "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"

# --- Process other transcriptomes ---
for trans in "${TRANSCRIPTOME_DIR}"/*.aa; do
    taxid_sample=$(basename "$trans" .aa)
    # Filter for near-complete sequences
    run_command "deeptmhmm_${taxid_sample}" ${DEEPTMHMM} -f "$trans" -o "${RESULTS_DIR}/chemogpcrs/deeptmhmm_${taxid_sample}"
    awk -v min_tm="${MIN_TM_REGIONS}" '$5 >= min_tm {print $1}' "${RESULTS_DIR}/chemogpcrs/deeptmhmm_${taxid_sample}/prediction" > "${RESULTS_DIR}/chemogpcrs/complete_ids_${taxid_sample}.txt"
    run_command "seqtk_complete_${taxid_sample}" ${SEQTK} subseq "$trans" "${RESULTS_DIR}/chemogpcrs/complete_ids_${taxid_sample}.txt" > "${RESULTS_DIR}/chemogpcrs/complete_${taxid_sample}.fa"
    # Search for GPCRs
    run_command "hhblits_${taxid_sample}" ${HHBLITS} -i "${RESULTS_DIR}/chemogpcrs/complete_${taxid_sample}.fa" -d "${RESULTS_DIR}/hhdb/references.hhm" -o "${RESULTS_DIR}/chemogpcrs/hhblits_${taxid_sample}.hhr" -e "${HHBLITS_EVALUE}" -cpu "${CPUS}"
    awk -v evalue="${HHBLITS_EVALUE}" '$5 < evalue {print $1}' "${RESULTS_DIR}/chemogpcrs/hhblits_${taxid_sample}.hhr" | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_${taxid_sample}.txt"
    run_command "extract_${taxid_sample}" ${SEQTK} subseq "${RESULTS_DIR}/chemogpcrs/complete_${taxid_sample}.fa" "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_${taxid_sample}.txt" > "${RESULTS_DIR}/chemogpcrs/chemogpcrs_${taxid_sample}.fa"
done

log "Chemoreceptive GPCR identification completed."
