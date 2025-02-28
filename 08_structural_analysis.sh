#!/bin/bash
# 08_structural_analysis.sh
# Purpose: Predict protein structures with AlphaFold, perform structural phylogeny with FoldTree, and compare with sequence trees.
# Inputs: Top-ranked candidates from step 07, reference structures from GPCRdb
# Outputs: Predicted structures (${RESULTS_DIR}/structural_analysis/alphafold/), structural trees (${RESULTS_DIR}/structural_analysis/*.tre), comparisons
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

#SBATCH --job-name=structural_analysis
#SBATCH --output=${LOGS_DIR}/08_structural_analysis_%j.out
#SBATCH --error=${LOGS_DIR}/08_structural_analysis_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
if ${GPU_ENABLED}; then
    #SBATCH --gres=gpu:1
fi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/structural_analysis/alphafold" "${RESULTS_DIR}/structural_analysis/references" "${LOGS_DIR}"

# Check dependency from step 07
if [ ! -f "${RESULTS_DIR}/step_completed_ranking.txt" ]; then
    log "Error: Candidate ranking step not completed."
    exit 1
fi

log "Starting structural analysis."

# Select top candidates for structural prediction
head -n 6 "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" | tail -n 5 > "${RESULTS_DIR}/structural_analysis/top_conserved.csv"
tail -n +2 "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" | awk -F, '$6 > 0.5' | head -n 5 > "${RESULTS_DIR}/structural_analysis/top_lse.csv"
cat "${RESULTS_DIR}/structural_analysis/top_conserved.csv" "${RESULTS_DIR}/structural_analysis/top_lse.csv" | cut -d, -f1 > "${RESULTS_DIR}/structural_analysis/top_ids.txt"
run_command "seqtk_top" ${SEQTK} subseq "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" "${RESULTS_DIR}/structural_analysis/top_ids.txt" > "${RESULTS_DIR}/structural_analysis/top_seqs.fa"

# Predict structures with AlphaFold
run_command "alphafold" ${ALPHAFOLD} --fasta_paths="${RESULTS_DIR}/structural_analysis/top_seqs.fa" --output_dir="${RESULTS_DIR}/structural_analysis/alphafold" --max_template_date=2023-01-01 --use_gpu=${GPU_ENABLED} --model_preset=monomer

# Fetch reference structures from GPCRdb
run_command "fetch_references" python3 "${SCRIPTS_DIR}/fetch_ligands.py" "${RESULTS_DIR}/structural_analysis/references" "${RESULTS_DIR}/structural_analysis/reference_ligands.csv"

# Perform TM-align comparisons
for berghia_pdb in "${RESULTS_DIR}/structural_analysis/alphafold/"*.pdb; do
    berghia_id=$(basename "$berghia_pdb" .pdb)
    for ref_pdb in "${RESULTS_DIR}/structural_analysis/references/"*.pdb; do
        ref_id=$(basename "$ref_pdb" .pdb)
        run_command "tmalign_${berghia_id}_${ref_id}" ${TMALIGN} "$berghia_pdb" "$ref_pdb" > "${RESULTS_DIR}/structural_analysis/tmalign_${berghia_id}_${ref_id}.txt"
    done
done

# Parse TM-align scores
run_command "tmalign_parse" awk '/TM-score=/ {print FILENAME "," $2}' "${RESULTS_DIR}/structural_analysis/tmalign_"*.txt > "${RESULTS_DIR}/structural_analysis/tmalign_scores.csv"

# Cluster structures with UPGMA
run_command "structural_clustering" python3 "${SCRIPTS_DIR}/cluster_structures.py" "${RESULTS_DIR}/structural_analysis/tmalign_scores.csv" "${RESULTS_DIR}/structural_analysis/clustered_tree.tre"

# Run FoldTree for structural phylogeny
run_command "foldtree" ${FOLDTREE} --input_dir "${RESULTS_DIR}/structural_analysis/alphafold" --output "${RESULTS_DIR}/structural_analysis/foldtree.tre" --method "${FOLDTREE_METHOD}"

# Plot structural analysis results
python3 "${SCRIPTS_DIR}/plot_heatmap.py" "${RESULTS_DIR}/structural_analysis/tmalign_scores.csv" "${RESULTS_DIR}/structural_analysis/tmalign_heatmap"
python3 "${SCRIPTS_DIR}/plot_pca.py" "${RESULTS_DIR}/structural_analysis/tmalign_scores.csv" "${RESULTS_DIR}/structural_analysis/tmalign_pca"
python3 "${SCRIPTS_DIR}/plot_struct_vs_seq.py" "${RESULTS_DIR}/structural_analysis/foldtree.tre" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" "${RESULTS_DIR}/structural_analysis/foldtree_vs_seq_plot"
python3 "${SCRIPTS_DIR}/plot_struct_vs_seq.py" "${RESULTS_DIR}/structural_analysis/clustered_tree.tre" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" "${RESULTS_DIR}/structural_analysis/clustered_vs_seq_plot"

# Estimate putative ligands
run_command "fetch_ligands" python3 "${SCRIPTS_DIR}/fetch_ligands.py" "${RESULTS_DIR}/structural_analysis/references" "${RESULTS_DIR}/structural_analysis/reference_ligands.csv" "${RESULTS_DIR}/structural_analysis/tmalign_scores.csv" "${RESULTS_DIR}/structural_analysis/clustered_tree.tre" "${RESULTS_DIR}/structural_analysis/foldtree.tre" "${RESULTS_DIR}/structural_analysis/putative_ligands.csv"

log "Structural analysis completed."
