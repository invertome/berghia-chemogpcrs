#!/bin/bash
# 05_selective_pressure_and_asr.sh
# Purpose: Analyze selective pressure with HyPhy's aBSREL and reconstruct ancestral sequences with FastML.
# Inputs: Orthogroups from step 03, alignments from step 04
# Outputs: Selective pressure results in ${RESULTS_DIR}/selective_pressure/, ASR sequences in ${RESULTS_DIR}/asr/
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=selective_pressure_asr
#SBATCH --output=${LOGS_DIR}/05_selective_pressure_asr_%j_%a.out
#SBATCH --error=${LOGS_DIR}/05_selective_pressure_asr_%j_%a.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=${DEFAULT_MEM}
#SBATCH --array=0-999%50  # Adjusted for orthogroup processing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/selective_pressure" "${RESULTS_DIR}/asr" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_lse_classification.txt"

log "Starting selective pressure and ASR analysis."

# Get orthogroups
ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/Orthogroups/OG"*.fa)
[ ${#ORTHOGROUPS[@]} -eq 0 ] && { log "Error: No orthogroups found"; exit 1; }

og="${ORTHOGROUPS[$SLURM_ARRAY_TASK_ID]}"
[ -z "$og" ] || [ ! -f "$og" ] && { log "Skipping missing orthogroup at index ${SLURM_ARRAY_TASK_ID}"; exit 0; }

base=$(basename "$og" .fa)
taxids=$(grep "^>" "$og" | sed 's/>//' | cut -d'_' -f1 | sort -u)
taxa_count=$(echo "$taxids" | wc -l)

# --- Selective Pressure with aBSREL ---
if [ "$taxa_count" -gt 1 ]; then
    protein_align="${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"
    nuc_align="${TRANSCRIPTOME_DIR}/$(grep "^>" "$protein_align" | head -1 | sed 's/>//' | cut -d'_' -f1,2).mrna"
    tree="${RESULTS_DIR}/phylogenies/protein/${base}.treefile"
    check_file "$protein_align" "$nuc_align" "$tree"
    run_command "${base}_codon" pal2nal.pl "$protein_align" "$nuc_align" -output paml > "${RESULTS_DIR}/selective_pressure/${base}_codon.phy"
    run_command "${base}_absrel" hyphy aBSREL --alignment "${RESULTS_DIR}/selective_pressure/${base}_codon.phy" --tree "$tree" --output "${RESULTS_DIR}/selective_pressure/${base}_absrel.json"
    python3 "${SCRIPTS_DIR}/parse_absrel.py" "${RESULTS_DIR}/selective_pressure/${base}_absrel.json" "${RESULTS_DIR}/selective_pressure/absrel_results.csv" || log "Warning: Failed to parse aBSREL for $base"
fi

# --- Berghia-specific LSEs ---
if [ "$(echo "$taxids" | grep -c "${BERGHIA_TAXID}")" -eq "$taxa_count" ] && [ "$taxa_count" -eq 1 ]; then
    protein_align="${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"
    nuc_align="${TRANSCRIPTOME_DIR}/taxid_berghia_berghia.mrna"
    tree="${RESULTS_DIR}/phylogenies/protein/${base}.treefile"
    check_file "$protein_align" "$nuc_align" "$tree"
    run_command "${base}_codon_lse" pal2nal.pl "$protein_align" "$nuc_align" -output paml > "${RESULTS_DIR}/selective_pressure/${base}_codon_lse.phy"
    run_command "${base}_absrel_lse" hyphy aBSREL --alignment "${RESULTS_DIR}/selective_pressure/${base}_codon_lse.phy" --tree "$tree" --output "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.json"
    python3 "${SCRIPTS_DIR}/parse_absrel.py" "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.json" "${RESULTS_DIR}/selective_pressure/absrel_results_lse.csv"
    
    # ASR for deep nodes
    deep_nodes=$(python3 "${SCRIPTS_DIR}/select_deep_nodes.py" "$tree" "${BERGHIA_TAXID}" "${MIN_ASR_DISTANCE}")
    for node in $deep_nodes; do
        run_command "${base}_${node}_asr" ${FASTML} --seq "$nuc_align" --tree "$tree" --out_seq "${RESULTS_DIR}/asr/${base}_${node}_asr.fa" --out_tree "${RESULTS_DIR}/asr/${base}_${node}_asr.tree" --node "$node" -t 8 --verbose
    done
    [ -f "${RESULTS_DIR}/asr/${base}_${deep_nodes%% *}_asr.fa" ] && python3 "${SCRIPTS_DIR}/plot_asr.py" "$tree" "${RESULTS_DIR}/asr/${base}_${deep_nodes%% *}_asr.fa" "${RESULTS_DIR}/asr/${base}_asr_plot"
fi

log "Selective pressure and ASR completed for ${base}."
