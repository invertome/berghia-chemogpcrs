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

# --- Pre-flight resource check ---
detect_resources

# --- All Berghia candidates + references ---
# Reference subsampling: when raw references exceed MAX_PHYLO_REFS, use CD-HIT clustering
# + taxonomy-weighted proportional allocation to select a diverse representative subset.
# This ensures the tree always includes reference sequences for meaningful phylo_scores.
MAX_PHYLO_REFS=${MAX_PHYLO_REFS:-2000}
REF_CLUSTER_IDENTITY=${REF_CLUSTER_IDENTITY:-0.7}
REF_LSE_WEIGHT=${REF_LSE_WEIGHT:-1.5}
REF_TAXONOMY_WEIGHTS=${REF_TAXONOMY_WEIGHTS:-"gastropoda:3.0,cephalopoda:1.5,bivalvia:1.5,other_molluscan_classes:1.2,annelida:1.0,platyhelminthes:1.0,other_lophotrochozoan_phyla:1.0"}

if [ ! -f "${RESULTS_DIR}/step_completed_all_berghia_refs_iqtree.txt" ]; then
    ALL_REFS="${RESULTS_DIR}/reference_sequences/all_references.fa"
    BERGHIA_FA="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"
    COMBINED="${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.fa"
    SUBSAMPLED_REFS="${RESULTS_DIR}/reference_sequences/subsampled_references.fa"
    SUBSAMPLED_RAW="${RESULTS_DIR}/reference_sequences/subsampled_refs_raw.fa"
    SUBSAMPLED_REPORT="${RESULTS_DIR}/reference_sequences/subsampling_report.json"
    SUBSAMPLED_ID_MAP="${RESULTS_DIR}/reference_sequences/subsampled_id_map.csv"
    NATH_REF_DIR="${REFERENCE_DIR}/nath_et_al"

    BERGHIA_COUNT=$(grep -c "^>" "$BERGHIA_FA")
    REF_COUNT=$(grep -c "^>" "$ALL_REFS" 2>/dev/null || echo 0)

    if [ "$REF_COUNT" -le "$MAX_PHYLO_REFS" ]; then
        # References fit within limit — use all
        cat "$BERGHIA_FA" "$ALL_REFS" > "$COMBINED"
        log "Combined Berghia ($BERGHIA_COUNT) + references ($REF_COUNT) = $((BERGHIA_COUNT + REF_COUNT)) sequences"
    elif [ -f "$SUBSAMPLED_REFS" ]; then
        # Reuse existing subsampled references (idempotent)
        SUB_COUNT=$(grep -c "^>" "$SUBSAMPLED_REFS")
        cat "$BERGHIA_FA" "$SUBSAMPLED_REFS" > "$COMBINED"
        log "Reusing subsampled references ($SUB_COUNT) + Berghia ($BERGHIA_COUNT) = $((BERGHIA_COUNT + SUB_COUNT)) sequences"
    elif [ -d "$NATH_REF_DIR" ]; then
        # Subsample: run CD-HIT + taxonomy-weighted selection on raw .faa files
        log "References ($REF_COUNT) exceed MAX_PHYLO_REFS ($MAX_PHYLO_REFS) — subsampling"
        run_command "subsample_references" python3 "${SCRIPTS_DIR}/subsample_references.py" \
            --ref-dir "$NATH_REF_DIR" \
            --target-size "$MAX_PHYLO_REFS" \
            --cluster-identity "$REF_CLUSTER_IDENTITY" \
            --taxonomy-weights "$REF_TAXONOMY_WEIGHTS" \
            --lse-weight "$REF_LSE_WEIGHT" \
            --cdhit-memory "${CDHIT_MEMORY:-8000}" \
            --threads "${CPUS}" \
            --output "$SUBSAMPLED_RAW" \
            --report "$SUBSAMPLED_REPORT"

        # Standardize headers (same pipeline as step 01)
        run_command "subsample_update_headers" python3 "${SCRIPTS_DIR}/update_headers.py" \
            "$SUBSAMPLED_RAW" "$SUBSAMPLED_ID_MAP" --source-type reference
        mv "${SUBSAMPLED_RAW}_updated.fa" "$SUBSAMPLED_REFS"

        SUB_COUNT=$(grep -c "^>" "$SUBSAMPLED_REFS")
        cat "$BERGHIA_FA" "$SUBSAMPLED_REFS" > "$COMBINED"
        log "Subsampled references ($SUB_COUNT) + Berghia ($BERGHIA_COUNT) = $((BERGHIA_COUNT + SUB_COUNT)) sequences"
    else
        # Fallback: no nath_et_al directory and refs too large — use Berghia only
        log "WARNING: References ($REF_COUNT) exceed limit and no nath_et_al/ directory found for subsampling"
        log "Using Berghia candidates only for overview tree"
        cp "$BERGHIA_FA" "$COMBINED"
    fi

    SEQ_COUNT=$(grep -c "^>" "$COMBINED")
    log "Building tree from $SEQ_COUNT sequences"

    # Check resource requirements for alignment and tree building
    log "Checking resource requirements for phylogenetic analysis..."
    get_dataset_stats "$COMBINED"
    check_resource_requirements "$COMBINED" alignment || \
        log --level=WARN "Proceeding despite alignment resource warning"

    run_command "all_berghia_refs_mafft" --stdout="${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa" ${MAFFT} --auto --thread "${CPUS}" "$COMBINED"
    check_alignment "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa" || { log "Error: Alignment quality check failed"; exit 1; }
    run_command "all_berghia_refs_clipkit" ${CLIPKIT} "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa" -m smart-gap -o "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa"

    # FastTree seed strategy: Generate approximate ML tree first to avoid local optima
    # This is important for large, divergent gene families like GPCRs
    run_command "all_berghia_refs_fasttree" --stdout="${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_fasttree.tre" ${FASTTREE} -lg -gamma "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa"

    # Check resource requirements for IQ-TREE (more memory-intensive than FastTree)
    check_resource_requirements "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa" iqtree || \
        log --level=WARN "Proceeding despite IQ-TREE resource warning"

    # Use FastTree result as starting tree for IQ-TREE (-t option)
    # IQ-TREE 3.x uses -T (threads) and --prefix instead of -nt and -pre
    run_command "all_berghia_refs_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -T "${CPUS}" -t "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_fasttree.tre" --prefix "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs"
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
# Parse LSE levels from config (same format as 03b_lse_classification.sh)
declare -A lse_taxids
for level in "${LSE_LEVELS[@]}"; do
    level_name=$(echo "$level" | cut -d':' -f1)
    taxids=$(echo "$level" | cut -d':' -f2 | tr ',' ' ')
    lse_taxids["$level_name"]="$taxids"
done

for level in "${!lse_taxids[@]}"; do
    if [ -f "${RESULTS_DIR}/lse_classification/lse_${level}.fa" ] && [ ! -f "${RESULTS_DIR}/step_completed_lse_${level}_iqtree.txt" ]; then
        mkdir -p "${RESULTS_DIR}/phylogenies/protein/lse_${level}"
        run_command "lse_${level}_mafft" --stdout="${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa" ${MAFFT} --auto --thread "${CPUS}" "${RESULTS_DIR}/lse_classification/lse_${level}.fa"
        check_alignment "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa" || { log "Error: Alignment quality check failed for lse_${level}"; exit 1; }
        run_command "lse_${level}_trimal" ${TRIMAL} -in "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa" -out "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" -automated1

        # FastTree seed strategy for LSE trees
        run_command "lse_${level}_fasttree" --stdout="${RESULTS_DIR}/phylogenies/protein/lse_${level}/fasttree.tre" ${FASTTREE} -lg -gamma "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa"
        run_command "lse_${level}_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -T "${CPUS}" -t "${RESULTS_DIR}/phylogenies/protein/lse_${level}/fasttree.tre" --prefix "${RESULTS_DIR}/phylogenies/protein/lse_${level}"

        python3 "${SCRIPTS_DIR}/visualize_tree.py" "${RESULTS_DIR}/phylogenies/protein/lse_${level}.treefile" "${RESULTS_DIR}/phylogenies/visualizations/lse_${level}"
    fi
done

# --- Orthogroup Trees ---
# Use manifest if available, otherwise fall back to globbing
MANIFEST_FILE="${RESULTS_DIR}/orthogroup_manifest.tsv"
OG_COUNT=$(get_orthogroup_count "$MANIFEST_FILE")

if [ "$OG_COUNT" -eq 0 ]; then
    log "Warning: No orthogroups found for phylogenetic analysis"
    # Still create completion flag since main tree may have been built
    if [ -f "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" ]; then
        touch "${RESULTS_DIR}/step_completed_04.txt"
    fi
    exit 0
fi

# Handle SLURM array indexing using new helper
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    validate_array_index "$SLURM_ARRAY_TASK_ID" "$OG_COUNT"

    # Get orthogroup name from manifest
    if [ -f "$MANIFEST_FILE" ]; then
        base=$(get_orthogroup_by_index "$SLURM_ARRAY_TASK_ID" "$MANIFEST_FILE")
    else
        # Fallback to globbing
        ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results"*/Orthogroups/OG*.fa)
        og="${ORTHOGROUPS[$SLURM_ARRAY_TASK_ID]}"
        base=$(basename "$og" .fa)
    fi
else
    # Non-array mode: process first orthogroup only (for testing)
    log "Running in non-array mode, processing first orthogroup only"
    if [ -f "$MANIFEST_FILE" ]; then
        base=$(get_orthogroup_by_index 0 "$MANIFEST_FILE")
    else
        ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results"*/Orthogroups/OG*.fa)
        base=$(basename "${ORTHOGROUPS[0]}" .fa)
    fi
fi

# Find the orthogroup FASTA file
og=$(find "${RESULTS_DIR}/orthogroups" -name "${base}.fa" -type f 2>/dev/null | head -1)
[ -z "$og" ] || [ ! -f "$og" ] && { log "Skipping missing orthogroup: ${base}"; exit 0; }

if [ ! -f "${RESULTS_DIR}/step_completed_${base}_iqtree.txt" ]; then
    run_command "${base}_mafft" --stdout="${RESULTS_DIR}/phylogenies/protein/${base}_aligned.fa" ${MAFFT} --auto --thread "${CPUS}" "$og"
    check_alignment "${RESULTS_DIR}/phylogenies/protein/${base}_aligned.fa" || { log "Error: Alignment quality check failed for ${base}"; exit 1; }
    run_command "${base}_trimal" ${TRIMAL} -in "${RESULTS_DIR}/phylogenies/protein/${base}_aligned.fa" -out "${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa" -automated1

    # FastTree seed strategy for orthogroup trees
    run_command "${base}_fasttree" --stdout="${RESULTS_DIR}/phylogenies/protein/${base}_fasttree.tre" ${FASTTREE} -lg -gamma "${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"
    run_command "${base}_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa" -m "${IQTREE_MODEL}" -B "${IQTREE_BOOTSTRAP}" -T "${CPUS}" -t "${RESULTS_DIR}/phylogenies/protein/${base}_fasttree.tre" --prefix "${RESULTS_DIR}/phylogenies/protein/${base}"

    # Create per-orthogroup completion flag
    touch "${RESULTS_DIR}/step_completed_${base}_iqtree.txt"

    # Also create array checkpoint for resume capability
    [ -n "$SLURM_ARRAY_TASK_ID" ] && create_array_checkpoint "04_phylo" "$SLURM_ARRAY_TASK_ID"
fi

python3 "${SCRIPTS_DIR}/visualize_tree.py" "${RESULTS_DIR}/phylogenies/protein/${base}.treefile" "${RESULTS_DIR}/phylogenies/visualizations/${base}"

log "Phylogenetic analysis completed for ${base}."

# Create overall completion flag when main tree is built (checked by downstream steps)
# Note: This flag is created even if not all array tasks are complete, since downstream
# steps only require the main all_berghia_refs tree
if [ -f "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" ]; then
    touch "${RESULTS_DIR}/step_completed_04.txt"
fi
