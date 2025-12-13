#!/bin/bash
# 04_phylogenetic_analysis.sh
# Purpose: Construct phylogenetic trees using IQ-TREE, Phyloformer, and optionally MrBayes, with alignment quality checks.
# Inputs: GPCR FASTA files from step 02, LSE FASTAs from step 03b, reference sequences from step 01
# Outputs: Phylogenetic trees in ${RESULTS_DIR}/phylogenies/protein/*.treefile, visualizations
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=phylogenetic_analysis
#SBATCH --output=${LOGS_DIR}/04_phylogenetic_analysis_%A_%a.out
#SBATCH --error=${LOGS_DIR}/04_phylogenetic_analysis_%A_%a.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --array=0-999%50  # Adjusted for orthogroup processing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/phylogenies/protein" "${RESULTS_DIR}/phylogenies/nucleotide" "${RESULTS_DIR}/phylogenies/visualizations" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_lse_classification.txt"

log "Starting phylogenetic analysis."

# --- Alignment Quality Check Function ---
# Validates alignment has sufficient sequence length and not too many gaps
check_alignment() {
    local file="$1"
    local min_len="${MIN_SEQ_LENGTH:-100}"
    local max_gap="${MAX_GAP_PERCENT:-50}"

    awk -v min_len="$min_len" -v max_gap="$max_gap" '
    BEGIN {
        seq_len = 0
        gaps = 0
        seq_count = 0
        total_len = 0
        total_gaps = 0
        failed = 0
    }
    /^>/ {
        # Process previous sequence if exists
        if (seq_len > 0) {
            seq_count++
            total_len += seq_len
            total_gaps += gaps
            gap_pct = (gaps / seq_len) * 100
            if (seq_len < min_len) {
                failed = 1
            }
        }
        # Reset for new sequence
        seq_len = 0
        gaps = 0
        next
    }
    !/^>/ {
        # Accumulate sequence stats (handles multi-line FASTA)
        line = $0
        seq_len += length(line)
        # Count gap characters
        n = gsub(/-/, "-", line)
        gaps += n
    }
    END {
        # Process last sequence
        if (seq_len > 0) {
            seq_count++
            total_len += seq_len
            total_gaps += gaps
            if (seq_len < min_len) {
                failed = 1
            }
        }

        # Check overall alignment quality
        if (seq_count == 0) {
            print "ERROR: No sequences found" > "/dev/stderr"
            exit 1
        }

        avg_gap_pct = (total_gaps / total_len) * 100
        if (avg_gap_pct > max_gap) {
            print "ERROR: Average gap percentage " avg_gap_pct "% exceeds " max_gap "%" > "/dev/stderr"
            exit 1
        }

        if (failed) {
            print "ERROR: Some sequences shorter than " min_len " residues" > "/dev/stderr"
            exit 1
        }

        # Success
        exit 0
    }
    ' "$file"
}

# --- All Berghia candidates + references ---
if [ ! -f "${RESULTS_DIR}/step_completed_all_berghia_refs_iqtree.txt" ]; then
    cat "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" "${RESULTS_DIR}/reference_sequences/all_references.fa" > "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.fa"
    run_command "all_berghia_refs_mafft" ${MAFFT} --auto --thread "${CPUS}" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.fa" > "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa"
    check_alignment "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa" || { log "Error: Alignment quality check failed"; exit 1; }
    run_command "all_berghia_refs_clipkit" ${CLIPKIT} smart-gap -i "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa" -o "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa"

    # FastTree seed strategy: Generate approximate ML tree first to avoid local optima
    # This is important for large, divergent gene families like GPCRs
    run_command "all_berghia_refs_fasttree" ${FASTTREE} -lg -gamma "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa" > "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_fasttree.tre"

    # Use FastTree result as starting tree for IQ-TREE (-t option)
    run_command "all_berghia_refs_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -nt "${CPUS}" -t "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_fasttree.tre" -pre "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs"
    run_command "phyloformer_all_berghia" python3 "${SCRIPTS_DIR}/test_phyloformer_models.py" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_phyloformer" "${CPUS}"
    if [ "$USE_MRBAYES" = true ]; then
        cat <<EOF > "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
begin mrbayes;
set autoclose=yes nowarn=yes;
execute ${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa;
lset rates=gamma;
mcmc ngen=1000000 samplefreq=100;
sump;
sumt;
end;
EOF
        run_command "all_berghia_refs_mrbayes" ${MRBAYES} -i "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.nex"
    fi
fi

# --- Visualizations for all_berghia_refs ---
python3 "${SCRIPTS_DIR}/visualize_tree.py" "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" "${RESULTS_DIR}/phylogenies/visualizations/all_berghia_refs"

# --- LSE Trees ---
for level in "${!lse_taxids[@]}"; do
    if [ -f "${RESULTS_DIR}/lse_classification/lse_${level}.fa" ] && [ ! -f "${RESULTS_DIR}/step_completed_lse_${level}_iqtree.txt" ]; then
        mkdir -p "${RESULTS_DIR}/phylogenies/protein/lse_${level}"
        run_command "lse_${level}_mafft" ${MAFFT} --auto --thread "${CPUS}" "${RESULTS_DIR}/lse_classification/lse_${level}.fa" > "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa"
        check_alignment "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa" || { log "Error: Alignment quality check failed for lse_${level}"; exit 1; }
        run_command "lse_${level}_trimal" ${TRIMAL} -in "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa" -out "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" -automated1

        # FastTree seed strategy for LSE trees
        run_command "lse_${level}_fasttree" ${FASTTREE} -lg -gamma "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" > "${RESULTS_DIR}/phylogenies/protein/lse_${level}/fasttree.tre"
        run_command "lse_${level}_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -nt "${CPUS}" -t "${RESULTS_DIR}/phylogenies/protein/lse_${level}/fasttree.tre" -pre "${RESULTS_DIR}/phylogenies/protein/lse_${level}"

        python3 "${SCRIPTS_DIR}/visualize_tree.py" "${RESULTS_DIR}/phylogenies/protein/lse_${level}.treefile" "${RESULTS_DIR}/phylogenies/visualizations/lse_${level}"
    fi
done

# --- Orthogroup Trees ---
ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results"*/Orthogroups/OG*.fa)

# Handle case where no orthogroups exist
if [ ${#ORTHOGROUPS[@]} -eq 0 ] || [ ! -f "${ORTHOGROUPS[0]}" ]; then
    log "Warning: No orthogroups found for phylogenetic analysis"
    exit 0
fi

# Handle SLURM array indexing - skip if index exceeds available orthogroups
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    if [ "$SLURM_ARRAY_TASK_ID" -ge ${#ORTHOGROUPS[@]} ]; then
        log "Skipping: Array index ${SLURM_ARRAY_TASK_ID} exceeds available orthogroups (${#ORTHOGROUPS[@]})"
        exit 0
    fi
    og="${ORTHOGROUPS[$SLURM_ARRAY_TASK_ID]}"
else
    # Non-array mode: process all orthogroups sequentially (for testing)
    log "Running in non-array mode, processing first orthogroup only"
    og="${ORTHOGROUPS[0]}"
fi

[ -z "$og" ] || [ ! -f "$og" ] && { log "Skipping missing orthogroup"; exit 0; }

base=$(basename "$og" .fa)
if [ ! -f "${RESULTS_DIR}/step_completed_${base}_iqtree.txt" ]; then
    run_command "${base}_mafft" ${MAFFT} --auto --thread "${CPUS}" "$og" > "${RESULTS_DIR}/phylogenies/protein/${base}_aligned.fa"
    check_alignment "${RESULTS_DIR}/phylogenies/protein/${base}_aligned.fa" || { log "Error: Alignment quality check failed for ${base}"; exit 1; }
    run_command "${base}_trimal" ${TRIMAL} -in "${RESULTS_DIR}/phylogenies/protein/${base}_aligned.fa" -out "${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa" -automated1

    # FastTree seed strategy for orthogroup trees
    run_command "${base}_fasttree" ${FASTTREE} -lg -gamma "${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa" > "${RESULTS_DIR}/phylogenies/protein/${base}_fasttree.tre"
    run_command "${base}_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -nt "${CPUS}" -t "${RESULTS_DIR}/phylogenies/protein/${base}_fasttree.tre" -pre "${RESULTS_DIR}/phylogenies/protein/${base}"
fi

python3 "${SCRIPTS_DIR}/visualize_tree.py" "${RESULTS_DIR}/phylogenies/protein/${base}.treefile" "${RESULTS_DIR}/phylogenies/visualizations/${base}"

log "Phylogenetic analysis completed for ${base}."
