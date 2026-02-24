#!/bin/bash
# 03a_busco_species_tree.sh
# Purpose: Generate a species tree using BUSCO single-copy orthologs with IQ-TREE for gene trees.
# Inputs: Transcriptome files in ${TRANSCRIPTOME_DIR}/*.aa, GPCR output from step 02
# Outputs: Species tree in ${SPECIES_TREE}, gene trees in ${RESULTS_DIR}/busco/gene_trees/
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=busco_species_tree
#SBATCH --output=${LOGS_DIR}/03a_busco_species_tree_%j.out
#SBATCH --error=${LOGS_DIR}/03a_busco_species_tree_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/busco" "${RESULTS_DIR}/busco/single_copy" "${RESULTS_DIR}/busco/alignments" "${RESULTS_DIR}/busco/gene_trees" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_extract_berghia.txt"

log "Starting BUSCO species tree generation."

# --- Run BUSCO on all transcriptomes ---
for trans in "${TRANSCRIPTOME_DIR}"/*.aa "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"; do
    sample=$(basename "$trans" .fa | sed 's/chemogpcrs_//')
    taxid_sample="${sample}"
    run_command "busco_${taxid_sample}" ${BUSCO} -i "$trans" -o "${RESULTS_DIR}/busco/busco_${taxid_sample}" -m transcriptome -l mollusca_odb10 -c "${CPUS}"
done

# --- Extract single-copy BUSCOs ---
for busco_dir in "${RESULTS_DIR}/busco/busco_"*; do
    [ -d "$busco_dir" ] || continue
    taxid_sample=$(basename "$busco_dir" | sed 's/busco_//')

    # Create per-sample directory for BUSCO sequences
    mkdir -p "${RESULTS_DIR}/busco/single_copy/${taxid_sample}"

    # Copy single-copy BUSCO sequences, prefixing with sample name
    for busco_faa in "${busco_dir}/run_mollusca_odb10/single_copy_busco_sequences/"*.faa; do
        [ -f "$busco_faa" ] || continue
        busco_name=$(basename "$busco_faa")
        cp "$busco_faa" "${RESULTS_DIR}/busco/single_copy/${taxid_sample}/${busco_name}"
    done

    [ -z "$(ls -A "${RESULTS_DIR}/busco/single_copy/${taxid_sample}" 2>/dev/null)" ] && log "Warning: No single-copy BUSCOs for ${taxid_sample}"
done

# --- Combine orthologous BUSCOs across species and align ---
# First, identify all unique BUSCO IDs present in multiple species
declare -A busco_files
for sample_dir in "${RESULTS_DIR}/busco/single_copy/"*/; do
    [ -d "$sample_dir" ] || continue
    sample=$(basename "$sample_dir")
    for busco_faa in "${sample_dir}"*.faa; do
        [ -f "$busco_faa" ] || continue
        busco_id=$(basename "$busco_faa" .faa)
        busco_files["$busco_id"]+="${busco_faa} "
    done
done

# Concatenate same-ID BUSCOs from different species and align
for busco_id in "${!busco_files[@]}"; do
    files=(${busco_files[$busco_id]})
    # Only process if present in at least 2 species
    if [ ${#files[@]} -ge 2 ]; then
        # Combine sequences, adding sample name to headers
        combined_file="${RESULTS_DIR}/busco/alignments/${busco_id}_combined.fa"
        > "$combined_file"
        for faa in "${files[@]}"; do
            sample=$(basename "$(dirname "$faa")")
            # Add sample prefix to sequence headers
            sed "s/^>/>${sample}_/" "$faa" >> "$combined_file"
        done

        run_command "align_busco_${busco_id}" --stdout="${RESULTS_DIR}/busco/alignments/${busco_id}.afa" ${MAFFT} --auto --thread "${CPUS}" "$combined_file"
        run_command "trim_busco_${busco_id}" ${TRIMAL} -in "${RESULTS_DIR}/busco/alignments/${busco_id}.afa" -out "${RESULTS_DIR}/busco/alignments/${busco_id}_trimmed.afa" -automated1
    fi
done

# --- Infer gene trees with IQ-TREE ---
for aln in "${RESULTS_DIR}/busco/alignments/"*_trimmed.afa; do
    base=$(basename "$aln" _trimmed.afa)
    run_command "tree_busco_${base}" ${IQTREE} -s "$aln" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -nt "${CPUS}" -pre "${RESULTS_DIR}/busco/gene_trees/${base}"
done

# --- Generate species tree with ASTRAL ---
find "${RESULTS_DIR}/busco/gene_trees/" -name "*.treefile" > "${RESULTS_DIR}/busco/gene_trees_list.txt" || { log "Error: No gene trees found"; exit 1; }
run_command "astral_species_tree" ${ASTRAL} -i "${RESULTS_DIR}/busco/gene_trees_list.txt" -o "${SPECIES_TREE}"

touch "${RESULTS_DIR}/step_completed_busco_species_tree.txt"
log "BUSCO species tree generation completed."
