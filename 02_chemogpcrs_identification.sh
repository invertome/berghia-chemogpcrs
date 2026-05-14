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

# --- Process Berghia transcriptome (HMM-first restructure, bead -m1f) ---
#
# Architecture: HMM-first GPCR filter -> TMbed only on GPCR-positive subset
# -> >=6 TM filter for chemoreceptors. The HMM-positive set IS the GPCR
# census for the follow-up Berghia GPCR/brain-expression paper.
#
# Why this order: previously we ran TMbed on the full 86k-protein
# transcriptome, then post-filtered for 6+ TM + HMM hits. 99%+ of TMbed
# compute was wasted on non-GPCRs and the few thousand multi-kDa transcript
# assembly outliers timed out TMbed via ProtT5's quadratic-memory slow-tail.
# HMM-first is faster (HMM scan is linear in length) and gives us the
# all-GPCRs census output for free.

mkdir -p "${RESULTS_DIR}/all_gpcrs"

# Step 1: HMM-first filter — full transcriptome -> GPCR-positive subset.
# Scans classification HMMs (curated bioamine/peptide/opsin/lipid/nucleotide/
# class-B/C/F) + Pfam fallback (7tm_1/2/3, Frizzled) + lse.hmm (mollusc
# chemoreceptor LSE OG HMMs). Census TSV preserved for the follow-up paper.
identify_gpcr_candidates "${TRANSCRIPTOME}" \
    "${RESULTS_DIR}/chemogpcrs/gpcrs_berghia.fa" \
    "${RESULTS_DIR}/all_gpcrs/berghia_gpcr_census.tsv"

[ -s "${RESULTS_DIR}/chemogpcrs/gpcrs_berghia.fa" ] || {
    log "Error: HMM-first filter found no GPCRs in Berghia transcriptome"; exit 1; }

# Step 2: Length pre-filter (defense-in-depth) — drops the rare LGR/ADGR/
# mGluR >MAX_AA_LENGTH from the TMbed input. Real chemoreceptors are
# 300-500 aa. The full GPCR-positive set (including long LGRs/ADGRs/mGluRs)
# remains in gpcrs_berghia.fa and the census TSV for the follow-up paper.
filter_fasta_by_length "${RESULTS_DIR}/chemogpcrs/gpcrs_berghia.fa" \
    "${RESULTS_DIR}/chemogpcrs/_tmbed_input_berghia.aa"

# Step 3: TMbed on GPCR-positive + length-OK subset (small, fast on GPU)
run_command "deeptmhmm_berghia" ${DEEPTMHMM} \
    -f "${RESULTS_DIR}/chemogpcrs/_tmbed_input_berghia.aa" \
    -o "${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia"

# Step 4: Apply >=6 TM region filter (chemoreceptor signature)
awk -v min_tm="${MIN_TM_REGIONS}" -v min_conf="${DEEPTMHMM_MIN_CONFIDENCE:-0.5}" \
    'NF >= 5 && $5+0 >= min_tm && $3+0 >= min_conf {print $1}' \
    "${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia/prediction" \
    > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt" || {
        log "Error: Failed to parse TMbed prediction"; exit 1; }

total_pred=$(wc -l < "${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia/prediction" 2>/dev/null || echo 0)
passed_pred=$(wc -l < "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt" 2>/dev/null || echo 0)
log "TMbed Berghia: ${passed_pred}/${total_pred} GPCR-positive sequences passed >=${MIN_TM_REGIONS} TM filter"

[ -s "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt" ] || {
    log "Error: No chemoreceptor candidates after TM filter"; exit 1; }

# Step 5: Extract chemoreceptor candidate sequences
run_command "extract_berghia" \
    --stdout="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" \
    ${SEQTK} subseq "${RESULTS_DIR}/chemogpcrs/gpcrs_berghia.fa" \
    "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_berghia.txt"

# --- Process additional reference transcriptomes (same HMM-first flow) ---
for trans in "${TRANSCRIPTOME_DIR}"/*.aa; do
    [ -f "$trans" ] || continue
    # Skip Berghia transcriptome (already processed above)
    [ "$(realpath "$trans")" = "$(realpath "${TRANSCRIPTOME}")" ] && continue
    sample=$(basename "$trans" .aa)
    taxid_sample="${sample}"

    # Step 1: HMM-first GPCR filter (no census for reference species — they
    # feed orthology/phylogeny stages, not the Berghia census paper)
    identify_gpcr_candidates "$trans" \
        "${RESULTS_DIR}/chemogpcrs/gpcrs_${taxid_sample}.fa"

    if [ ! -s "${RESULTS_DIR}/chemogpcrs/gpcrs_${taxid_sample}.fa" ]; then
        log --level=WARN "${taxid_sample}: no GPCR-positive sequences; skipping TMbed"
        : > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_${taxid_sample}.txt"
        : > "${RESULTS_DIR}/chemogpcrs/chemogpcrs_${taxid_sample}.fa"
        touch "${RESULTS_DIR}/step_completed_extract_${taxid_sample}.txt"
        continue
    fi

    # Step 2: Length pre-filter (defense-in-depth)
    filter_fasta_by_length "${RESULTS_DIR}/chemogpcrs/gpcrs_${taxid_sample}.fa" \
        "${RESULTS_DIR}/chemogpcrs/_tmbed_input_${taxid_sample}.aa"

    # Step 3: TMbed on GPCR-positive + length-OK subset
    run_command "deeptmhmm_${taxid_sample}" ${DEEPTMHMM} \
        -f "${RESULTS_DIR}/chemogpcrs/_tmbed_input_${taxid_sample}.aa" \
        -o "${RESULTS_DIR}/chemogpcrs/deeptmhmm_${taxid_sample}"

    # Step 4: Apply >=6 TM filter
    awk -v min_tm="${MIN_TM_REGIONS}" -v min_conf="${DEEPTMHMM_MIN_CONFIDENCE:-0.5}" \
        'NF >= 5 && $5+0 >= min_tm && $3+0 >= min_conf {print $1}' \
        "${RESULTS_DIR}/chemogpcrs/deeptmhmm_${taxid_sample}/prediction" \
        > "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_${taxid_sample}.txt"
    log "TMbed ${taxid_sample}: $(wc -l < "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_${taxid_sample}.txt") chemoreceptor candidates"

    # Step 5: Extract chemoreceptor candidate sequences
    run_command "extract_${taxid_sample}" \
        --stdout="${RESULTS_DIR}/chemogpcrs/chemogpcrs_${taxid_sample}.fa" \
        ${SEQTK} subseq "${RESULTS_DIR}/chemogpcrs/gpcrs_${taxid_sample}.fa" \
        "${RESULTS_DIR}/chemogpcrs/chemogpcr_ids_${taxid_sample}.txt"
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
