#!/bin/bash
# 08_structural_analysis.sh
# Purpose: Predict structures with AlphaFold, build structural phylogenies including references, and compare with sequence trees.
# Inputs: Ranked candidates from step 07, reference structures from GPCRdb
# Outputs: Predicted structures in ${RESULTS_DIR}/structural_analysis/alphafold/, structural trees
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=structural_analysis
#SBATCH --output=${LOGS_DIR}/08_structural_analysis_%j.out
#SBATCH --error=${LOGS_DIR}/08_structural_analysis_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
[ "$GPU_ENABLED" = true ] && echo "#SBATCH --gres=gpu:1"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/structural_analysis/alphafold" "${RESULTS_DIR}/structural_analysis/references" "${RESULTS_DIR}/structural_analysis/all_pdb" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_rank_candidates.txt"

log "Starting structural analysis."

# --- Select diverse candidates ---
run_command "select_diverse" python3 "${SCRIPTS_DIR}/select_diverse_candidates.py" \
    "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" \
    "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" \
    "${NUM_STRUCTURAL_CANDIDATES}" \
    "${RESULTS_DIR}/structural_analysis/top_ids.txt"

# --- Include ASR sequences ---
cat "${RESULTS_DIR}/asr/"*_asr.fa > "${RESULTS_DIR}/structural_analysis/asr_seqs.fa" 2>/dev/null
grep "^>" "${RESULTS_DIR}/structural_analysis/asr_seqs.fa" | sed 's/>//' >> "${RESULTS_DIR}/structural_analysis/top_ids.txt" 2>/dev/null

# --- Deduplicate IDs ---
sort -u "${RESULTS_DIR}/structural_analysis/top_ids.txt" > "${RESULTS_DIR}/structural_analysis/top_ids_unique.txt"
mv "${RESULTS_DIR}/structural_analysis/top_ids_unique.txt" "${RESULTS_DIR}/structural_analysis/top_ids.txt"

# --- Extract sequences for AlphaFold ---
run_command "seqtk_top" ${SEQTK} subseq "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" "${RESULTS_DIR}/structural_analysis/top_ids.txt" > "${RESULTS_DIR}/structural_analysis/top_seqs.fa"

# --- Predict structures with AlphaFold ---
run_command "alphafold" ${ALPHAFOLD} --fasta_paths="${RESULTS_DIR}/structural_analysis/top_seqs.fa" --output_dir="${RESULTS_DIR}/structural_analysis/alphafold" --max_template_date=2023-01-01 --use_gpu=${GPU_ENABLE} --model_preset=monomer

# --- Fetch reference structures and metadata ---
run_command "fetch_references" python3 "${SCRIPTS_DIR}/fetch_ligands.py" \
    "${RESULTS_DIR}/structural_analysis/references" \
    "${RESULTS_DIR}/structural_analysis/reference_ligands.csv" \
    "${GPCRDB_SEARCH_TERMS}" \
    "${GPCRDB_SPECIES}"

# --- Combine predicted and reference PDBs for structural phylogeny ---
cp "${RESULTS_DIR}/structural_analysis/alphafold/"*.pdb "${RESULTS_DIR}/structural_analysis/all_pdb/" || log "Warning: No AlphaFold PDBs found"
cp "${RESULTS_DIR}/structural_analysis/references/"*.pdb "${RESULTS_DIR}/structural_analysis/all_pdb/" || log "Warning: No reference PDBs found"

# --- Build structural phylogeny with FoldTree ---
run_command "foldtree" ${FOLDTREE} --input_dir "${RESULTS_DIR}/structural_analysis/all_pdb" --output "${RESULTS_DIR}/structural_analysis/foldtree.tre" --method "${FOLDTREE_METHOD}"

# --- Compare structural and sequence trees ---
python3 "${SCRIPTS_DIR}/plot_struct_vs_seq.py" "${RESULTS_DIR}/structural_analysis/foldtree.tre" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" "${RESULTS_DIR}/structural_analysis/struct_vs_seq_plot"

log "Structural analysis completed."
