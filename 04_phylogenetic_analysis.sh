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
# Bead hto: SLURM does not expand command/variable substitution in #SBATCH
# directives, so '$(scale_resources)' was silently ignored and a direct
# `sbatch 04_...sh` got cluster-default resources (≈1 CPU) for a large IQ-TREE
# job. Static defaults below (match config CPUS=16 / DEFAULT_MEM=64G); the
# production wrapper sbatch_run_04.sh sets higher values and takes precedence.
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
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

# --- Phase gating (bead vo8.1) ---
# Stage 04 builds two kinds of products:
#   1. "globals": the per-class A/B/C/F trees + LSE-level trees — built ONCE.
#   2. "ogs":     per-orthogroup trees — one per SLURM array task.
# The global products are NOT array-safe: their completion markers are checked
# then written non-atomically, so 50 array tasks starting together would all
# see "no marker" and race to build the same class_*/ files concurrently.
# STAGE04_PHASE selects which half runs:
#   globals -> build globals, then exit (run as a single non-array job)
#   ogs     -> skip globals, run this array task's OG tree only
#   all     -> both (default; for non-array / local runs)
# The production wrapper submits a non-array globals job, then the OG array
# with --dependency=afterok on it.
STAGE04_PHASE="${STAGE04_PHASE:-all}"
case "${STAGE04_PHASE}" in
    globals|ogs|all) ;;
    *) log "ERROR: invalid STAGE04_PHASE='${STAGE04_PHASE}' (want globals|ogs|all)"; exit 2 ;;
esac
log "STAGE04_PHASE=${STAGE04_PHASE}"

# Bead -m6k: Build IQ-TREE bootstrap-flag string. Use UFBoot + SH-aLRT
# always; add --tbe (Transfer Bootstrap Expectation, Lemoine 2018) when
# IQTREE_TBE=1 in config.sh — robust to rogue taxa, important for paralog-
# rich GPCR trees.
IQTREE_BOOT_FLAGS="-B ${IQTREE_BOOTSTRAP} -alrt 1000"
if [ "${IQTREE_TBE:-1}" = "1" ]; then
    IQTREE_BOOT_FLAGS+=" --tbe"
fi

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

# ===================== GLOBALS PHASE (globals|all) =====================
# Per-class global trees + LSE-level trees. Built once; skipped when
# STAGE04_PHASE=ogs (the OG array tasks rely on a prior globals job).
if [ "${STAGE04_PHASE}" != "ogs" ]; then

# --- Per-class global trees (A, B, C, F) ---
# Replaces the former single "all_berghia_refs" tree. Each class gets its own
# combined FASTA (P2 pool refs + Berghia subset + outgroup), filter stack,
# ClipKit trim, FastTree seed, and IQ-TREE run. Completion is gated per-class
# so individual classes can be resumed independently.
#
# Per-class pool dir and Berghia class TSV are set in config.sh.
# Outgroup swap-map: A<-C, B<-A, C<-A, F<-A (locked decision 2026-05-28).
#
# Back-compat: after the loop the Class A outputs are symlinked to the legacy
# all_berghia_refs.* paths so stage 04b / 07 / 08 consumers keep working.
BERGHIA_FA="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"

for CLASS in ${GPCR_CLASSES:-A B C F}; do
    CLASS_MARKER="${RESULTS_DIR}/step_completed_class_${CLASS}_iqtree.txt"
    if [ -f "${CLASS_MARKER}" ]; then
        log "Class ${CLASS}: skipping (marker present)"
        continue
    fi

    log "=== Class ${CLASS}: building global tree ==="

    CLASS_DIR="${RESULTS_DIR}/phylogenies/protein/class_${CLASS}"
    mkdir -p "${CLASS_DIR}"

    # --- 1. Assemble class input FASTA: refs + Berghia(class) + outgroup ---
    REFS_FA="${PER_CLASS_POOL_DIR}/refs_class_${CLASS}.fa"
    if [ ! -f "${REFS_FA}" ]; then
        log "ERROR: missing P2 pool ${REFS_FA} — run P2 (per-class pool build) first"
        exit 1
    fi

    # Berghia subset for this class — extract from the full Berghia FASTA
    # using the class column from the P1 classifier output.
    BERGHIA_CLASS_FA="${CLASS_DIR}/berghia_class_${CLASS}.fa"
    if [ -f "${BERGHIA_CLASS_TSV}" ] && [ -f "${BERGHIA_FA}" ]; then
        # Extract sequence IDs for this class from TSV (col 1=seq_id, col 2=class; skip header)
        awk -F'\t' -v cls="${CLASS}" 'NR>1 && $2==cls {print $1}' \
            "${BERGHIA_CLASS_TSV}" > "${CLASS_DIR}/berghia_class_${CLASS}_ids.txt"
        N_IDS=$(wc -l < "${CLASS_DIR}/berghia_class_${CLASS}_ids.txt")
        if [ "${N_IDS}" -gt 0 ]; then
            # Use seqtk (already on PATH per config.sh SEQTK) for ID-based extraction
            "${SEQTK:-seqtk}" subseq "${BERGHIA_FA}" \
                "${CLASS_DIR}/berghia_class_${CLASS}_ids.txt" > "${BERGHIA_CLASS_FA}"
            log "Class ${CLASS}: extracted ${N_IDS} Berghia sequences"
        else
            log "WARN: no Berghia sequences assigned to class ${CLASS} in ${BERGHIA_CLASS_TSV}"
            : > "${BERGHIA_CLASS_FA}"
        fi
    else
        log "WARN: ${BERGHIA_CLASS_TSV} not found — class ${CLASS} tree skips Berghia inclusion"
        : > "${BERGHIA_CLASS_FA}"
    fi

    # Outgroup: first ${OUTGROUP_BUDGET_PER_CLASS} sequences from the
    # swap-map source class's ref pool (locked decision 2026-05-28).
    OUTGROUP_SOURCE_VAR="OUTGROUP_SOURCE_CLASS_${CLASS}"
    OUTGROUP_SOURCE_CLASS="${!OUTGROUP_SOURCE_VAR}"
    OUTGROUP_SOURCE_FA="${PER_CLASS_POOL_DIR}/refs_class_${OUTGROUP_SOURCE_CLASS}.fa"
    OUTGROUP_FA="${CLASS_DIR}/outgroup_class_${CLASS}.fa"
    if [ -f "${OUTGROUP_SOURCE_FA}" ]; then
        # Take the first N complete FASTA records (header + sequence block)
        awk -v n="${OUTGROUP_BUDGET_PER_CLASS:-10}" \
            '/^>/{c++} c>n{exit} {print}' \
            "${OUTGROUP_SOURCE_FA}" > "${OUTGROUP_FA}"
        OG_COUNT=$(grep -c '^>' "${OUTGROUP_FA}" 2>/dev/null || echo 0)
        log "Class ${CLASS}: outgroup from class ${OUTGROUP_SOURCE_CLASS} = ${OG_COUNT} seqs"
    else
        log "WARN: no outgroup source ${OUTGROUP_SOURCE_FA} for class ${CLASS}; tree will be unrooted"
        : > "${OUTGROUP_FA}"
    fi

    # Combine into one FASTA for the filter stack
    COMBINED="${CLASS_DIR}/class_${CLASS}_combined.fa"
    cat "${REFS_FA}" "${BERGHIA_CLASS_FA}" "${OUTGROUP_FA}" > "${COMBINED}"
    N_INPUT=$(grep -c '^>' "${COMBINED}")
    log "Class ${CLASS}: combined input = ${N_INPUT} sequences (refs + Berghia + outgroup, cap 3000)"

    # Length filter: remove extreme outliers before alignment
    COMBINED_FILTERED="${CLASS_DIR}/class_${CLASS}_combined_filtered.fa"
    python3 "${SCRIPTS_DIR}/filter_sequences_by_length.py" \
        "${COMBINED}" "${COMBINED_FILTERED}" \
        --min-length "${SEQ_LENGTH_FILTER_MIN:-250}" \
        --max-length-method "${SEQ_LENGTH_FILTER_MAX_METHOD:-tukey}" \
        --max-length-floor "${SEQ_LENGTH_FILTER_MAX_FLOOR:-800}" \
        --report "${CLASS_DIR}/class_${CLASS}_length_filter_report.tsv"
    COMBINED="${COMBINED_FILTERED}"
    N_INPUT=$(grep -c '^>' "${COMBINED}")
    log "Class ${CLASS}: after length filter: ${N_INPUT} sequences"

    # Resource pre-flight
    get_dataset_stats "${COMBINED}"
    check_resource_requirements "${COMBINED}" alignment || \
        log --level=WARN "Class ${CLASS}: proceeding despite alignment resource warning"

    # --- 2. Filter stack (PREQUAL → ensemble → CLOAK → TAPER) ---
    # MAFFT_DASH=1 is the default (config.sh); run_alignment_ensemble.sh and
    # run_aligner.sh both respect MAFFT_DASH env, so --dash is applied to all
    # MAFFT calls within the ensemble automatically.
    run_alignment_filter_stack \
        "${COMBINED}" \
        "${CLASS_DIR}/class_${CLASS}_aligned.fa" \
        "${CLASS_DIR}/_filter_stack_class_${CLASS}" \
        "class_${CLASS}" \
        "${CPUS}" || { log "Error: filter stack failed for class ${CLASS}"; exit 1; }
    check_alignment "${CLASS_DIR}/class_${CLASS}_aligned.fa" || \
        { log "Error: alignment quality check failed for class ${CLASS}"; exit 1; }

    # --- 3. ClipKit kpic-smart-gap ---
    # kpic-smart-gap: drops both heavily-gappy columns (smart-gap) AND
    # singleton variants (kpic), correctly cleaning up CLOAK's super-
    # alignment artifact while preserving parsimony-informative + constant
    # sites for IQ-TREE branch-length estimation.
    # shellcheck disable=SC2086
    run_command "class_${CLASS}_clipkit" \
        ${CLIPKIT} "${CLASS_DIR}/class_${CLASS}_aligned.fa" \
        -m kpic-smart-gap \
        -o "${CLASS_DIR}/class_${CLASS}_trimmed.fa"

    # Drop near-all-gap rows defensively (see per-OG comment below)
    python3 "${SCRIPTS_DIR}/drop_near_all_gap_rows.py" \
        --input  "${CLASS_DIR}/class_${CLASS}_trimmed.fa" \
        --output "${CLASS_DIR}/class_${CLASS}_trimmed.fa" \
        2>>"${LOGS_DIR}/class_${CLASS}_drop_gappy.log" || true

    # --- 4. FastTree seed → IQ-TREE ---
    run_command "class_${CLASS}_fasttree" \
        --stdout="${CLASS_DIR}/class_${CLASS}_fasttree.tre" \
        ${FASTTREE} -seed "${FASTTREE_SEED}" -lg -gamma \
        "${CLASS_DIR}/class_${CLASS}_trimmed.fa"

    check_resource_requirements "${CLASS_DIR}/class_${CLASS}_trimmed.fa" iqtree || \
        log --level=WARN "Class ${CLASS}: proceeding despite IQ-TREE resource warning"

    class_iqtree_seed_arg=""
    if [ -s "${CLASS_DIR}/class_${CLASS}_fasttree.tre" ]; then
        class_iqtree_seed_arg="-t ${CLASS_DIR}/class_${CLASS}_fasttree.tre"
    fi
    # shellcheck disable=SC2086
    run_command "class_${CLASS}_iqtree" \
        ${IQTREE} \
        -s "${CLASS_DIR}/class_${CLASS}_trimmed.fa" \
        -st AA \
        -m "${IQTREE_MODEL_FIND}" -mset "${IQTREE_MODEL_SET}" \
        ${IQTREE_BOOT_FLAGS} \
        -seed "${IQTREE_SEED}" \
        -T "${CPUS}" \
        ${class_iqtree_seed_arg} \
        --prefix "${CLASS_DIR}/class_${CLASS}"

    if [ "${USE_MRBAYES:-false}" = true ]; then
        cat <<MBEOF > "${CLASS_DIR}/class_${CLASS}.nex"
begin mrbayes;
set autoclose=yes nowarn=yes;
execute ${CLASS_DIR}/class_${CLASS}_trimmed.fa;
lset rates=gamma;
mcmc ngen=1000000 samplefreq=100;
sump;
sumt;
end;
MBEOF
        run_command "class_${CLASS}_mrbayes" ${MRBAYES} -i "${CLASS_DIR}/class_${CLASS}.nex"
    fi

    touch "${CLASS_MARKER}"
    log "=== Class ${CLASS}: done ==="
done

# --- Class A back-compat aliases ---
# Stage 04b, 07, and 08 reference all_berghia_refs.* paths. Symlink Class A
# outputs to those canonical names so existing consumers work without
# modification until P4 updates them.
if [ -f "${RESULTS_DIR}/phylogenies/protein/class_A/class_A.treefile" ]; then
    ln -sf "class_A/class_A.treefile" \
        "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile"
    ln -sf "class_A/class_A_aligned.fa" \
        "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned.fa"
    ln -sf "class_A/class_A_trimmed.fa" \
        "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_trimmed.fa"
    # Also alias canonical sibling (used by stage 04b ECL divergence)
    if [ -f "${RESULTS_DIR}/phylogenies/protein/class_A/class_A_aligned_canonical.fa" ]; then
        ln -sf "class_A/class_A_aligned_canonical.fa" \
            "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs_aligned_canonical.fa"
    fi
    log "Class A back-compat aliases created (all_berghia_refs.* -> class_A/*)"
fi

# --- Visualizations for Class A (back-compat visualization path) ---
python3 "${SCRIPTS_DIR}/visualize_tree.py" \
    "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" \
    "${RESULTS_DIR}/phylogenies/visualizations/all_berghia_refs"

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
        # Bead -i61: full filter stack (PREQUAL -> ensemble -> CLOAK -> TAPER).
        run_alignment_filter_stack \
            "${RESULTS_DIR}/lse_classification/lse_${level}.fa" \
            "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa" \
            "${RESULTS_DIR}/phylogenies/protein/lse_${level}/_filter_stack" \
            "lse_${level}" \
            "${CPUS}" || { log "Error: filter stack failed for lse_${level}"; exit 1; }
        check_alignment "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa" || { log "Error: Alignment quality check failed for lse_${level}"; exit 1; }
        # ClipKit kpic-smart-gap (consistent with global path; correctly handles
        # CLOAK super-alignment by dropping singleton-variant + heavily-gappy cols).
        run_command "lse_${level}_clipkit" ${CLIPKIT} "${RESULTS_DIR}/phylogenies/protein/lse_${level}/aligned.fa" -m kpic-smart-gap -o "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa"

        # Drop near-all-gap rows defensively (see comment in the global path).
        python3 "${SCRIPTS_DIR}/drop_near_all_gap_rows.py" \
            --input "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" \
            --output "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" \
            2>>"${LOGS_DIR}/lse_${level}_drop_gappy.log" || true

        # FastTree seed strategy for LSE trees
        run_command "lse_${level}_fasttree" --stdout="${RESULTS_DIR}/phylogenies/protein/lse_${level}/fasttree.tre" ${FASTTREE} -seed "${FASTTREE_SEED}" -lg -gamma "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa"
        # Guard the FastTree seed file (see global path comment) and force -st AA.
        lse_iqtree_seed_arg=""
        if [ -s "${RESULTS_DIR}/phylogenies/protein/lse_${level}/fasttree.tre" ]; then
            lse_iqtree_seed_arg="-t ${RESULTS_DIR}/phylogenies/protein/lse_${level}/fasttree.tre"
        fi
        # shellcheck disable=SC2086
        run_command "lse_${level}_iqtree" ${IQTREE} -s "${RESULTS_DIR}/phylogenies/protein/lse_${level}/trimmed.fa" -st AA -m "${IQTREE_MODEL_FIND}" -mset "${IQTREE_MODEL_SET}" ${IQTREE_BOOT_FLAGS} -seed "${IQTREE_SEED}" -T "${CPUS}" ${lse_iqtree_seed_arg} --prefix "${RESULTS_DIR}/phylogenies/protein/lse_${level}"

        python3 "${SCRIPTS_DIR}/visualize_tree.py" "${RESULTS_DIR}/phylogenies/protein/lse_${level}.treefile" "${RESULTS_DIR}/phylogenies/visualizations/lse_${level}"
    fi
done

# --- End globals phase ---
if [ "${STAGE04_PHASE}" = "globals" ]; then
    # Signal globals done so 04b/07/08 can proceed on the Class A tree even
    # while the OG array is still running.
    touch "${RESULTS_DIR}/step_completed_04_globals.txt"
    if [ -f "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" ]; then
        touch "${RESULTS_DIR}/step_completed_04.txt"
    fi
    log "STAGE04_PHASE=globals complete — globals built, exiting before OG array"
    exit 0
fi

fi  # ===================== END GLOBALS PHASE =====================

# ===================== OG ARRAY PHASE (ogs|all) =====================
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

# --- Per-OG class tagging ---
# Route each OG into its class subdirectory so outputs are co-located with
# the per-class global tree. Reads og_class_majority.tsv when available
# (produced by a future P3b/stage-03 step); defaults to class A with a
# warning when the file is absent (back-compat for runs without per-OG
# class assignment).
OG_CLASS_TSV="${RESULTS_DIR}/classification/og_class_majority.tsv"
if [ -f "${OG_CLASS_TSV}" ]; then
    OG_CLASS=$(awk -F'\t' -v og="${base}" 'NR>1 && $1==og {print $2; exit}' "${OG_CLASS_TSV}")
    if [ -z "${OG_CLASS}" ]; then
        log --level=WARN "OG '${base}' not found in ${OG_CLASS_TSV}; routing to class_unclassified"
        OG_CLASS="unclassified"
    fi
else
    log "WARN: ${OG_CLASS_TSV} not found; routing ${base} to class_A (back-compat)"
    OG_CLASS="A"
fi
OG_OUT_DIR="${RESULTS_DIR}/phylogenies/protein/class_${OG_CLASS}/${base}"
mkdir -p "${OG_OUT_DIR}"

if [ ! -f "${RESULTS_DIR}/step_completed_${base}_iqtree.txt" ]; then
    # Bead -i61 (May 2026 v2): full filter stack per OG.
    # MAFFT_DASH=1 is inherited from the environment (config.sh default);
    # run_alignment_ensemble.sh and run_aligner.sh apply --dash automatically.
    run_alignment_filter_stack \
        "$og" \
        "${OG_OUT_DIR}/${base}_aligned.fa" \
        "${OG_OUT_DIR}/_filter_stack_${base}" \
        "${base}" \
        "${CPUS}" || { log "Error: filter stack failed for ${base}"; exit 1; }
    check_alignment "${OG_OUT_DIR}/${base}_aligned.fa" || { log "Error: Alignment quality check failed for ${base}"; exit 1; }
    # ClipKit kpic-smart-gap (consistent with global+LSE paths; correctly
    # handles CLOAK super-alignment).
    # shellcheck disable=SC2086
    run_command "${base}_clipkit" \
        ${CLIPKIT} "${OG_OUT_DIR}/${base}_aligned.fa" \
        -m kpic-smart-gap \
        -o "${OG_OUT_DIR}/${base}_trimmed.fa"

    # Drop near-all-gap rows defensively (see comment in the global path).
    python3 "${SCRIPTS_DIR}/drop_near_all_gap_rows.py" \
        --input  "${OG_OUT_DIR}/${base}_trimmed.fa" \
        --output "${OG_OUT_DIR}/${base}_trimmed.fa" \
        2>>"${LOGS_DIR}/${base}_drop_gappy.log" || true

    # FastTree seed strategy for orthogroup trees
    run_command "${base}_fasttree" \
        --stdout="${OG_OUT_DIR}/${base}_fasttree.tre" \
        ${FASTTREE} -seed "${FASTTREE_SEED}" -lg -gamma \
        "${OG_OUT_DIR}/${base}_trimmed.fa"
    # Guard the FastTree seed file (see global path comment) and force -st AA.
    og_iqtree_seed_arg=""
    if [ -s "${OG_OUT_DIR}/${base}_fasttree.tre" ]; then
        og_iqtree_seed_arg="-t ${OG_OUT_DIR}/${base}_fasttree.tre"
    fi
    # shellcheck disable=SC2086
    run_command "${base}_iqtree" \
        ${IQTREE} \
        -s "${OG_OUT_DIR}/${base}_trimmed.fa" \
        -st AA \
        -m "${IQTREE_MODEL_FIND}" -mset "${IQTREE_MODEL_SET}" \
        ${IQTREE_BOOT_FLAGS} \
        -seed "${IQTREE_SEED}" \
        -T "${CPUS}" \
        ${og_iqtree_seed_arg} \
        --prefix "${OG_OUT_DIR}/${base}"

    # Bead -iof: TreeShrink rogue-taxon cleaning. Replaces the tree in place
    # and leaves a rollback copy at <tree>.original.treefile.
    if [ "${RUN_TREESHRINK:-1}" = "1" ] && [ -f "${OG_OUT_DIR}/${base}.treefile" ]; then
        bash "${SCRIPTS_DIR}/run_treeshrink.sh" \
            --single-tree="${OG_OUT_DIR}/${base}.treefile" \
            --output-dir="${OG_OUT_DIR}/treeshrink" \
            --quantile="${TREESHRINK_QUANTILE:-0.05}" \
            2>> "${LOGS_DIR}/treeshrink_${base}.err" \
            || log --level=WARN "TreeShrink failed for ${base} (kept original tree)"
    fi

    # Create per-orthogroup completion flag
    touch "${RESULTS_DIR}/step_completed_${base}_iqtree.txt"

    # Also create array checkpoint for resume capability
    [ -n "$SLURM_ARRAY_TASK_ID" ] && create_array_checkpoint "04_phylo" "$SLURM_ARRAY_TASK_ID"
fi

python3 "${SCRIPTS_DIR}/visualize_tree.py" \
    "${OG_OUT_DIR}/${base}.treefile" \
    "${RESULTS_DIR}/phylogenies/visualizations/${base}"

log "Phylogenetic analysis completed for ${base} (class ${OG_CLASS})."

# Create overall completion flag when main Class A tree is built (checked by
# downstream steps via the back-compat all_berghia_refs.treefile symlink).
if [ -f "${RESULTS_DIR}/phylogenies/protein/all_berghia_refs.treefile" ]; then
    touch "${RESULTS_DIR}/step_completed_04.txt"
fi
