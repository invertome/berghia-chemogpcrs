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

# Check dependency (step 01 creates step_completed_01.txt)
check_file "${RESULTS_DIR}/step_completed_01.txt"

log "Starting chemoreceptive GPCR identification."

# --- Process Berghia transcriptome ---
run_command "deeptmhmm_berghia" ${DEEPTMHMM} -f "${TRANSCRIPTOME}" -o "${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia"
# Filter by TM region count AND confidence score
# DeepTMHMM output format: ID, prediction_type, confidence, ..., n_tm_regions
# Column 3 is confidence, column 5 is TM region count
awk -v min_tm="${MIN_TM_REGIONS}" -v min_conf="${DEEPTMHMM_MIN_CONFIDENCE:-0.5}" \
    'NF >= 5 && $5+0 >= min_tm && $3+0 >= min_conf {print $1}' \
    "${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia/prediction" > "${RESULTS_DIR}/chemogpcrs/complete_ids_berghia.txt" || { log "Error: Failed to parse DeepTMHMM"; exit 1; }

# Log filtering stats
total_pred=$(wc -l < "${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia/prediction" 2>/dev/null || echo 0)
passed_pred=$(wc -l < "${RESULTS_DIR}/chemogpcrs/complete_ids_berghia.txt" 2>/dev/null || echo 0)
log "DeepTMHMM: ${passed_pred}/${total_pred} sequences passed filters (>=${MIN_TM_REGIONS} TM regions, >=${DEEPTMHMM_MIN_CONFIDENCE:-0.5} confidence)"
run_command "seqtk_complete_berghia" ${SEQTK} subseq "${TRANSCRIPTOME}" "${RESULTS_DIR}/chemogpcrs/complete_ids_berghia.txt" > "${RESULTS_DIR}/chemogpcrs/complete_berghia.fa"

# --- HHblits search (optional - requires proper HH-suite database) ---
# hhmake requires an MSA input, not a multi-FASTA. Building a proper HH-suite
# database requires ffindex. Skip HHblits if tools are unavailable or skipped.
HHBLITS_AVAILABLE=false
if command -v "${HHMAKE%% *}" &>/dev/null && command -v "${HHBLITS%% *}" &>/dev/null \
   && [[ "${HHMAKE}" != *SKIPPED* ]] && [[ "${HHBLITS}" != *SKIPPED* ]]; then
    # Check if we have an aligned reference to build HHM from
    if [ -f "${RESULTS_DIR}/reference_sequences/conserved_references_aligned.fa" ]; then
        run_command "hhdb_creation" ${HHMAKE} -i "${RESULTS_DIR}/reference_sequences/conserved_references_aligned.fa" -o "${RESULTS_DIR}/hhdb/references.hhm" -v 1
        HHBLITS_AVAILABLE=true
    else
        log "Warning: No aligned references for HHM building, skipping HHblits"
    fi
fi

touch "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hhblits_berghia.txt"
if [ "$HHBLITS_AVAILABLE" = true ] && [ -f "${RESULTS_DIR}/hhdb/references.hhm" ]; then
    run_command "hhblits_berghia" ${HHBLITS} -i "${RESULTS_DIR}/chemogpcrs/complete_berghia.fa" -d "${RESULTS_DIR}/hhdb/references.hhm" -o "${RESULTS_DIR}/chemogpcrs/hhblits_berghia.hhr" -e "${HHBLITS_EVALUE}" -cpu "${CPUS}"
    # Parse HHR format: Extract query IDs from hits with E-value below threshold
    python3 -c "
import sys, re
evalue_thresh = float('${HHBLITS_EVALUE}')
current_query = None
with open('${RESULTS_DIR}/chemogpcrs/hhblits_berghia.hhr') as f:
    for line in f:
        if line.startswith('Query'):
            current_query = line.split()[1]
        elif re.match(r'^\s*\d+\s+\S+', line) and current_query:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    evalue = float(parts[4])
                    if evalue < evalue_thresh:
                        print(current_query)
                except (ValueError, IndexError):
                    pass
" | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hhblits_berghia.txt" || log "Warning: HHblits parsing failed"
else
    log "HHblits skipped — using HMMSEARCH as primary identification method"
fi

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
    [ -f "$trans" ] || continue
    # Skip Berghia transcriptome (already processed above)
    [ "$(realpath "$trans")" = "$(realpath "${TRANSCRIPTOME}")" ] && continue
    sample=$(basename "$trans" .aa)
    taxid_sample="${sample}"
    run_command "deeptmhmm_${taxid_sample}" ${DEEPTMHMM} -f "$trans" -o "${RESULTS_DIR}/chemogpcrs/deeptmhmm_${taxid_sample}"
    # Filter by TM region count AND confidence score
    awk -v min_tm="${MIN_TM_REGIONS}" -v min_conf="${DEEPTMHMM_MIN_CONFIDENCE:-0.5}" \
        'NF >= 5 && $5+0 >= min_tm && $3+0 >= min_conf {print $1}' \
        "${RESULTS_DIR}/chemogpcrs/deeptmhmm_${taxid_sample}/prediction" > "${RESULTS_DIR}/chemogpcrs/complete_ids_${taxid_sample}.txt"
    log "DeepTMHMM ${taxid_sample}: $(wc -l < "${RESULTS_DIR}/chemogpcrs/complete_ids_${taxid_sample}.txt") sequences passed filters"
    run_command "seqtk_complete_${taxid_sample}" ${SEQTK} subseq "$trans" "${RESULTS_DIR}/chemogpcrs/complete_ids_${taxid_sample}.txt" > "${RESULTS_DIR}/chemogpcrs/complete_${taxid_sample}.fa"
    
    touch "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hhblits_${taxid_sample}.txt"
    if [ "$HHBLITS_AVAILABLE" = true ] && [ -f "${RESULTS_DIR}/hhdb/references.hhm" ]; then
        run_command "hhblits_${taxid_sample}" ${HHBLITS} -i "${RESULTS_DIR}/chemogpcrs/complete_${taxid_sample}.fa" -d "${RESULTS_DIR}/hhdb/references.hhm" -o "${RESULTS_DIR}/chemogpcrs/hhblits_${taxid_sample}.hhr" -e "${HHBLITS_EVALUE}" -cpu "${CPUS}"
        python3 -c "
import sys, re
evalue_thresh = float('${HHBLITS_EVALUE}')
current_query = None
with open('${RESULTS_DIR}/chemogpcrs/hhblits_${taxid_sample}.hhr') as f:
    for line in f:
        if line.startswith('Query'):
            current_query = line.split()[1]
        elif re.match(r'^\s*\d+\s+\S+', line) and current_query:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    evalue = float(parts[4])
                    if evalue < evalue_thresh:
                        print(current_query)
                except (ValueError, IndexError):
                    pass
" | sort -u > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_hhblits_${taxid_sample}.txt" || true
    fi
    
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
    touch "${RESULTS_DIR}/step_completed_extract_${taxid_sample}.txt"
done

# --- Register GPCR sequences in the ID map for downstream taxid lookups ---
# References are registered in step 01; target species sequences need registration too
export ID_MAP BERGHIA_TAXID
python3 << 'PYEOF'
import os, sys

id_map = os.environ.get('ID_MAP', '')
if not id_map or not os.path.exists(id_map):
    print("Warning: ID map not found, skipping target registration", file=sys.stderr)
    sys.exit(0)

# Collect all GPCR FASTA files and their taxids
entries = []

# Berghia
berghia_fasta = os.path.join(os.environ.get('RESULTS_DIR', ''), 'chemogpcrs', 'chemogpcrs_berghia.fa')
berghia_taxid = os.environ.get('BERGHIA_TAXID', '')
if os.path.exists(berghia_fasta) and berghia_taxid:
    with open(berghia_fasta) as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]
                entries.append((seq_id, berghia_taxid, 'Berghia stephanieae'))

# Load existing IDs to avoid duplicates
existing_ids = set()
with open(id_map) as f:
    for line in f:
        if line.strip():
            existing_ids.add(line.split(',')[1] if ',' in line else '')

# Append new entries
added = 0
with open(id_map, 'a') as out:
    for seq_id, taxid, species in entries:
        if seq_id not in existing_ids:
            out.write(f"{seq_id},{seq_id},{taxid},target,{species},,0\n")
            added += 1

print(f"Registered {added} target sequences in ID map", file=sys.stderr)
PYEOF

touch "${RESULTS_DIR}/step_completed_extract_berghia.txt"

# Create main step completion flag for downstream dependencies
touch "${RESULTS_DIR}/step_completed_02.txt"
log "Chemoreceptive GPCR identification completed."
