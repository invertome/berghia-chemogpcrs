#!/bin/bash
# 06_synteny_and_mapping.sh
# Purpose: Map transcripts to genomes and perform synteny analyses.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

#SBATCH --job-name=synteny_mapping
#SBATCH --output=${LOGS_DIR}/06_synteny_mapping_%j.out
#SBATCH --error=${LOGS_DIR}/06_synteny_mapping_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

mkdir -p "${RESULTS_DIR}/mapping" "${RESULTS_DIR}/synteny" "${LOGS_DIR}"

if [ ! -f "${RESULTS_DIR}/step_completed_extract_berghia.txt" ]; then
    log "Error: Chemoreceptive GPCR identification step not completed."
    exit 1
fi

log "Starting synteny and mapping analysis."

run_command "map_berghia_minimap2" ${MINIMAP2} -ax splice -uf --secondary=no "${GENOME}" "${TRANSCRIPTOME}" > "${RESULTS_DIR}/mapping/berghia_transcripts.sam"
run_command "map_berghia_samtools" ${SAMTOOLS} view -bS "${RESULTS_DIR}/mapping/berghia_transcripts.sam" | ${SAMTOOLS} sort -o "${RESULTS_DIR}/mapping/berghia_transcripts_sorted.bam"
run_command "map_berghia_index" ${SAMTOOLS} index "${RESULTS_DIR}/mapping/berghia_transcripts_sorted.bam"

for genome in "${GENOME_DIR}"/*/*.fa; do
    if [ "$genome" != "${GENOME}" ]; then
        taxid_sample=$(basename "$(dirname "$genome")")
        run_command "synteny_mcscanx_${taxid_sample}" ${MCSCANX} -a -e "${BLAST_EVALUE}" -m 25 -s 5 "${GENOME}" "$genome" "${RESULTS_DIR}/synteny/berghia_vs_${taxid_sample}_mcscanx"
        run_command "synteny_iadhore_${taxid_sample}" ${IADHORE} -i "${GENOME}" -j "$genome" -o "${RESULTS_DIR}/synteny/berghia_vs_${taxid_sample}_iadhore" -c "${DATA_DIR}/iadhore_config.ini"
    fi
done

for col_file in "${RESULTS_DIR}/synteny/"*/*_mcscanx.collinearity; do
    grep "GPCR" "$col_file" | awk '{print $1}' | sort -u >> "${RESULTS_DIR}/synteny/gpcr_synteny_ids.txt"
done
python3 "${SCRIPTS_DIR}/plot_synteny.py" "${RESULTS_DIR}/synteny/berghia_vs_*_mcscanx.collinearity" "${RESULTS_DIR}/synteny/synteny_plot.png"

log "Synteny and mapping analysis completed."
