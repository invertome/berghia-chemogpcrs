#!/bin/bash
# 03a_busco_species_tree.sh
# Purpose: Generate a species tree using BUSCO genes with ASTRAL and OrthoFinder for gene tree reconciliation.
# Inputs: Transcriptome files (${TRANSCRIPTOME}, ${TRANSCRIPTOME_DIR}/*.aa, GPCR output from step 02)
# Outputs: ASTRAL species tree (${SPECIES_TREE}), BUSCO gene trees (${RESULTS_DIR}/busco/orthofinder/Results*/Gene_Trees/)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

#SBATCH --job-name=busco_species_tree
#SBATCH --output=${LOGS_DIR}/03a_busco_species_tree_%j.out
#SBATCH --error=${LOGS_DIR}/03a_busco_species_tree_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH --cpus-per-task=16  # Cap for IQ-TREE/ASTRAL
#SBATCH --mem=${DEFAULT_MEM}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/busco" "${RESULTS_DIR}/busco/single_copy" "${RESULTS_DIR}/busco/alignments" "${RESULTS_DIR}/busco/gene_trees" "${RESULTS_DIR}/busco/orthofinder" "${LOGS_DIR}"

# Check dependency from step 02
if [ ! -f "${RESULTS_DIR}/step_completed_extract_berghia.txt" ]; then
    log "Error: Chemoreceptive GPCR identification step not completed."
    exit 1
fi

log "Starting BUSCO species tree generation."

# --- Run BUSCO on all transcriptomes to identify single-copy orthologs ---
for trans in "${TRANSCRIPTOME_DIR}"/*.aa "${TRANSCRIPTOME}" "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"; do
    taxid_sample=$(basename "$trans" .fa | sed 's/chemogpcrs_//')
    run_command "busco_${taxid_sample}" ${BUSCO} -i "$trans" -o "${RESULTS_DIR}/busco/busco_${taxid_sample}" -m transcriptome -l mollusca_odb10 -c "${CPUS}"
done

# --- Extract complete single-copy BUSCOs ---
for busco_dir in "${RESULTS_DIR}/busco/busco_"*; do
    taxid_sample=$(basename "$busco_dir" | sed 's/busco_//')
    cp "${busco_dir}/run_mollusca_odb10/single_copy_busco_sequences/"*.faa "${RESULTS_DIR}/busco/single_copy/${taxid_sample}_" 2>/dev/null || log "No single-copy BUSCOs found for ${taxid_sample}, skipping."
done

# --- Align and trim BUSCOs, infer gene trees ---
for busco in "${RESULTS_DIR}/busco/single_copy/"*.faa; do
    busco_id=$(basename "$busco" .faa | cut -d'_' -f2-)
    run_command "align_busco_${busco_id}" ${MAFFT} --auto --thread "${CPUS}" "$busco" > "${RESULTS_DIR}/busco/alignments/${busco_id}.afa"
    run_command "trim_busco_${busco_id}" ${TRIMAL} -in "${RESULTS_DIR}/busco/alignments/${busco_id}.afa" -out "${RESULTS_DIR}/busco/alignments/${busco_id}_trimmed.afa" -automated1
    run_command "tree_busco_${busco_id}" ${IQTREE} -s "${RESULTS_DIR}/busco/alignments/${busco_id}_trimmed.afa" -m "${IQTREE_MODEL}" -B 100 -nt 16 -pre "${RESULTS_DIR}/busco/gene_trees/${busco_id}"
done

# --- Run OrthoFinder on BUSCO genes for gene tree reconciliation ---
cat "${RESULTS_DIR}/busco/single_copy/"*.faa > "${RESULTS_DIR}/busco/orthofinder/busco_input.fa"
run_command "orthofinder_busco" ${ORTHOFINDER} -f "${RESULTS_DIR}/busco/orthofinder/busco_input.fa" -t "${CPUS}" -S diamond -M msa -A mafft -T fasttree

# --- Combine gene trees into species tree with ASTRAL ---
find "${RESULTS_DIR}/busco/gene_trees/" -name "*.treefile" > "${RESULTS_DIR}/busco/gene_trees_list.txt"
run_command "astral_species_tree" ${ASTRAL} -i "${RESULTS_DIR}/busco/gene_trees_list.txt" -o "${RESULTS_DIR}/busco/busco_species_tree.tre"

log "BUSCO species tree generation completed."
