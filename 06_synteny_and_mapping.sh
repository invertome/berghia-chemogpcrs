#!/bin/bash
# 06_synteny_and_mapping.sh
# Purpose: Analyze synteny across available genomes using MCScanX; mapping is optional if transcriptomes exist.
# Inputs: Genomes in ${GENOME_DIR}/*.fasta, transcriptomes for mapping if available
# Outputs: Synteny results in ${RESULTS_DIR}/synteny/, optional mapping BAMs in ${RESULTS_DIR}/mapping/
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=synteny_mapping
#SBATCH --output=${LOGS_DIR}/06_synteny_and_mapping_%j.out
#SBATCH --error=${LOGS_DIR}/06_synteny_and_mapping_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/mapping" "${RESULTS_DIR}/synteny" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_extract_berghia.txt"

log "Starting synteny and mapping analysis."

# --- Optional mapping of transcriptomes to genomes ---
for genome in "${GENOME_DIR}"/*.fasta; do
    taxid_sample=$(basename "$genome" .fasta)
    trans="${TRANSCRIPTOME_DIR}/${taxid_sample}.aa"
    if [ -f "$trans" ]; then
        run_command "minimap2_${taxid_sample}" ${MINIMAP2} -ax splice -uf -k14 "$genome" "$trans" > "${RESULTS_DIR}/mapping/${taxid_sample}.sam"
        run_command "samtools_${taxid_sample}" ${SAMTOOLS} view -bS "${RESULTS_DIR}/mapping/${taxid_sample}.sam" | ${SAMTOOLS} sort -o "${RESULTS_DIR}/mapping/${taxid_sample}.bam"
        run_command "samtools_index_${taxid_sample}" ${SAMTOOLS} index "${RESULTS_DIR}/mapping/${taxid_sample}.bam"
    else
        log "Note: No transcriptome for $taxid_sample, skipping mapping"
    fi
done

# --- Synteny analysis with MCScanX for all available genomes ---
for genome1 in "${GENOME_DIR}"/*.fasta; do
    for genome2 in "${GENOME_DIR}"/*.fasta; do
        if [ "$genome1" != "$genome2" ]; then
            taxid1=$(basename "$genome1" .fasta)
            taxid2=$(basename "$genome2" .fasta)
            run_command "synteny_${taxid1}_vs_${taxid2}" ${MCSCANX} -a -b 2 "$genome1" "$genome2" "${RESULTS_DIR}/synteny/${taxid1}_vs_${taxid2}"
        fi
    done
done

# --- Plot synteny at different taxonomic levels ---
python3 "${SCRIPTS_DIR}/plot_synteny.py" "${RESULTS_DIR}/synteny" "${RESULTS_DIR}/synteny/synteny_plot" || log "Warning: Synteny plotting failed"

log "Synteny and mapping analysis completed."
