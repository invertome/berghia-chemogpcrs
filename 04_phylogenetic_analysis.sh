#!/bin/bash
# 04_phylogenetic_analysis.sh
# Purpose: Build and refine phylogenetic trees at multiple levels using optimized methods.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

#SBATCH --job-name=phylogenetic_analysis
#SBATCH --output=${LOGS_DIR}/04_phylogenetic_analysis_%A_%a.out
#SBATCH --error=${LOGS_DIR}/04_phylogenetic_analysis_%A_%a.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH --cpus-per-task=16  # Cap for IQ-TREE/MrBayes
#SBATCH --mem=${DEFAULT_MEM}
#SBATCH --array=0-999%100   # Increased concurrency
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

mkdir -p "${RESULTS_DIR}/phylogenies/protein" "${RESULTS_DIR}/phylogenies/nucleotide" "${RESULTS_DIR}/phylogenies/visualizations" "${LOGS_DIR}"

if [ ! -f "${RESULTS_DIR}/step_completed_lse_classification.txt" ]; then
    log "Error: LSE classification step not completed."
    exit 1
fi

log "Starting phylogenetic analysis."

# Quality control function
check_alignment() {
    local file="$1"
    local min_len="${MIN_SEQ_LENGTH}"
    local max_gap="${MAX_GAP_PERCENT}"
    awk -v min_len="$min_len" -v max_gap="$max_gap" '
    BEGIN { seq_len = 0; gaps = 0 }
    /^>/ { if (seq_len > 0) {
        if (seq_len < min_len || (gaps/seq_len)*100 > max_gap) exit 1
        seq_len = 0; gaps = 0
    } }
    !/^>/ { seq_len += length($0); gaps += gsub("-", "-", $0) }
    END { if (seq_len < min_len || (gaps/seq_len)*100 > max_gap) exit 1 }
    ' "$file"
}

# All Berghia candidates + references
if [ ! -f "${RESULTS_DIR}/step_completed_all_berghia_refs_iqtree.txt" ]; then
    cat "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" "${RESULTS_DIR}/reference_sequences/all_references.fa" > "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.fa"
    run_command "all_berghia_refs_mafft" ${MAFFT} --auto --thread "${CPUS}" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.fa" > "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa"
    check_alignment "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa" || { log "Error: Initial alignment quality check failed for all_berghia_refs"; exit 1; }
    run_command "all_berghia_refs_trimal" ${TRIMAL} -in "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa" -out "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa" -automated1
    run_command "all_berghia_refs_bmge" ${BMGE} -i "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa" -t AA -o "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_bmge.fa"
    python3 "${SCRIPTS_DIR}/prune_alignment.py" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_bmge.fa" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_pruned.fa" "${MIN_SEQ_LENGTH}" "${MAX_GAP_PERCENT}"
    run_command "all_berghia_refs_fasttree" ${FASTTREE} "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_pruned.fa" > "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_fasttree.tre"
    run_command "all_berghia_refs_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_pruned.fa" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -nt 16 -pre "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs" -t "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_fasttree.tre"
    if [ "$USE_MRBAYES" = true ]; then
        echo "begin mrbayes;" > "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
        echo "set autoclose=yes nowarn=yes;" >> "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
        echo "execute ${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_pruned.fa;" >> "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
        echo "lset rates=gamma;" >> "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
        echo "mcmc ngen=1000000 samplefreq=100;" >> "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
        echo "sump;" >> "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
        echo "sumt;" >> "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
        echo "end;" >> "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
        run_command "all_berghia_refs_mrbayes" ${MRBAYES} -i "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
    fi
    run_command "all_berghia_refs_phyloformer" ${PHYLOFORMER} -i "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_pruned.fa" -o "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_phyloformer.tre" --model LG+G --threads "${CPUS}"
fi

python3 "${SCRIPTS_DIR}/visualize_tree.py" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" "${RESULTS_DIR}/phylogenies/visualizations/all_berghia_refs_iqtree"
python3 "${SCRIPTS_DIR}/visualize_tree.py" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_phyloformer.tre" "${RESULTS_DIR}/phylogenies/visualizations/all_berghia_refs_phyloformer"
python3 "${SCRIPTS_DIR}/plot_phyloformer_iqtree.py" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_phyloformer.tre" "${RESULTS_DIR}/phylogenies/visualizations/all_berghia_refs_comparison.png"

# Multi-level LSE trees
for level in "aeolids" "nudibranchs" "gastropods"; do
    if [ -f "${RESULTS_DIR}/lse_classification/lse_${level}.fa" ] && [ ! -f "${RESULTS_DIR}/step_completed_lse_${level}_iqtree.txt" ]; then
        mkdir -p "${RESULTS_DIR}/phylogenies/protein/lse_${level}"
        run_command "lse_${level}_mafft" ${MAFFT} --auto --thread "${CPUS}" "${RESULTS_DIR}/lse_classification/lse_${level}.fa" > "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa"
        check_alignment "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa" || { log "Error: Alignment quality check failed for lse_${level}"; exit 1; }
        run_command "lse_${level}_trimal" ${TRIMAL} -in "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa" -out "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" -automated1
        run_command "lse_${level}_fasttree" ${FASTTREE} "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" > "${RESULTS_DIR}/phylogenies/protein/lse_${level}/fasttree.tre"
        run_command "lse_${level}_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -nt 16 -pre "${RESULTS_DIR}/phylogenies/protein/lse_${level}" -t "${RESULTS_DIR}/phylogenies/protein/lse_${level}/fasttree.tre"
        python3 "${SCRIPTS_DIR}/visualize_tree.py" "${RESULTS_DIR}/phylogenies/protein/lse_${level}.treefile" "${RESULTS_DIR}/phylogenies/visualizations/lse_${level}"
    fi
done

# Orthogroup-specific trees
ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/Orthogroups/OG"*.fa)
if [ ${#ORTHOGROUPS[@]} -eq 0 ]; then
    log "Error: No orthogroups found."
    exit 1
fi

og="${ORTHOGROUPS[$SLURM_ARRAY_TASK_ID]}"
if [ -z "$og" ] || [ ! -f "$og" ]; then
    log "Skipping empty or missing orthogroup at index ${SLURM_ARRAY_TASK_ID}."
    exit 0
fi

base=$(basename "$og" .fa")
if [ ! -f "${RESULTS_DIR}/step_completed_${base}_iqtree.txt" ]; then
    run_command "${base}_mafft" ${MAFFT} --auto --thread "${CPUS}" "$og" > "${RESULTS_DIR}/phylogenies/protein/${base}_aligned.fa"
    check_alignment "${RESULTS_DIR}/phylogenies/protein/${base}_aligned.fa" || { log "Error: Alignment quality check failed for ${base}"; exit 1; }
    run_command "${base}_trimal" ${TRIMAL} -in "${RESULTS_DIR}/phylogenies/protein/${base}_aligned.fa" -out "${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa" -automated1
    run_command "${base}_fasttree" ${FASTTREE} "${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa" > "${RESULTS_DIR}/phylogenies/protein/${base}_fasttree.tre"
    run_command "${base}_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -nt 16 -pre "${RESULTS_DIR}/phylogenies/protein/${base}" -t "${RESULTS_DIR}/phylogenies/protein/${base}_fasttree.tre"
fi

if [ ! -f "${RESULTS_DIR}/step_completed_${base}_nuc_iqtree.txt" ]; then
    nuc_og="${RESULTS_DIR}/phylogenies/nucleotide/${base}_nuc.fa"
    run_command "${base}_nuc_seqtk" ${SEQTK} subseq "${TRANSCRIPTOME}" <(awk '{print $1}' "$og") > "$nuc_og"
    run_command "${base}_nuc_mafft" ${MAFFT} --auto --thread "${CPUS}" "$nuc_og" > "${RESULTS_DIR}/phylogenies/nucleotide/${base}_aligned.fa"
    check_alignment "${RESULTS_DIR}/phylogenies/nucleotide/${base}_aligned.fa" || { log "Error: Alignment quality check failed for ${base}_nuc"; exit 1; }
    run_command "${base}_nuc_trimal" ${TRIMAL} -in "${RESULTS_DIR}/phylogenies/nucleotide/${base}_aligned.fa" -out "${RESULTS_DIR}/phylogenies/nucleotide/${base}_trimmed.fa" -automated1
    run_command "${base}_nuc_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/nucleotide/${base}_trimmed.fa" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -nt 16 -pre "${RESULTS_DIR}/phylogenies/nucleotide/${base}"
fi

python3 "${SCRIPTS_DIR}/visualize_tree.py" "${RESULTS_DIR}/phylogenies/protein/${base}.treefile" "${RESULTS_DIR}/phylogenies/visualizations/${base}"

log "Phylogeny for ${base} completed."
