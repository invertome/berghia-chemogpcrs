#!/bin/bash
# 07_candidate_ranking.sh
# Purpose: Rank GPCR candidates based on phylogeny, selection, expression, and synteny data.
# Inputs: GPCR IDs from step 02, expression data, results from steps 04-06
# Outputs: Ranked candidates CSV in ${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv, plots
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=candidate_ranking
#SBATCH --output=${LOGS_DIR}/07_candidate_ranking_%j.out
#SBATCH --error=${LOGS_DIR}/07_candidate_ranking_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directory
mkdir -p "${RESULTS_DIR}/ranking" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check required dependencies
check_file "${RESULTS_DIR}/step_completed_extract_berghia.txt" "${PHYLO_DIR}/${PHYLO_TREE_FILENAME}"

# Check optional dependencies (warn but don't fail)
check_file --warn-only "${RESULTS_DIR}/selective_pressure/absrel_results.csv"
if [ ! -f "${RESULTS_DIR}/selective_pressure/absrel_results.csv" ]; then
    log "Note: aBSREL results not found - dN/dS scoring will be skipped"
fi

log "Starting candidate ranking."

# Create expression output directory
mkdir -p "${RESULTS_DIR}/expression" "${RESULTS_DIR}/gproteins" "${RESULTS_DIR}/ecl_analysis"

# --- Phase 1: Process expression data if available ---
if [ -d "${SALMON_QUANT_DIR}" ]; then
    log "Processing expression data from ${SALMON_QUANT_DIR}..."
    python3 "${SCRIPTS_DIR}/process_expression.py" \
        --quant-dir "${SALMON_QUANT_DIR}" \
        --chemosensory-tissues "${CHEMOSENSORY_TISSUES}" \
        --other-tissues "${NON_CHEMOSENSORY_TISSUES}" \
        --min-tpm "${MIN_TPM_THRESHOLD}" \
        --tau-threshold "${TAU_THRESHOLD}" \
        --output "${RESULTS_DIR}/expression/expression_summary.csv" \
        || log "Warning: Expression processing failed or no data available"
else
    log "Note: No expression data directory found at ${SALMON_QUANT_DIR}"
fi

# --- Phase 2: Classify G-proteins if references available ---
if [ -f "${GPROTEIN_REF_FASTA}" ]; then
    log "Classifying G-proteins in transcriptomes..."
    python3 "${SCRIPTS_DIR}/classify_gproteins.py" \
        --transcriptomes "${TRANSCRIPTOME_DIR}/*.aa" \
        --reference "${GPROTEIN_REF_FASTA}" \
        --classes "${GPROTEIN_REF_CLASSES}" \
        --mollusc-reference "${GPROTEIN_MOLLUSC_FASTA}" \
        --mollusc-classes "${GPROTEIN_MOLLUSC_CLASSES}" \
        --evalue "${GPROTEIN_EVALUE}" \
        --output "${RESULTS_DIR}/gproteins/classified_gproteins.csv" \
        || log "Warning: G-protein classification failed or no references available"

    # Calculate co-expression if expression data available
    if [ -f "${RESULTS_DIR}/expression/expression_summary.csv" ]; then
        log "Analyzing G-protein co-expression..."
        python3 "${SCRIPTS_DIR}/coexpression_analysis.py" \
            --expression "${RESULTS_DIR}/expression/expression_summary.csv" \
            --gproteins "${RESULTS_DIR}/gproteins/classified_gproteins.csv" \
            --candidates "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" \
            --chemosensory-tissues "${CHEMOSENSORY_TISSUES}" \
            --output "${RESULTS_DIR}/gproteins/gprotein_coexpression.csv" \
            || log "Warning: G-protein co-expression analysis failed"
    fi
else
    log "Note: No G-protein reference found at ${GPROTEIN_REF_FASTA}"
fi

# --- Phase 3: Analyze ECL divergence ---
if [ -d "${RESULTS_DIR}/deeptmhmm" ] && [ -d "${RESULTS_DIR}/phylogenies" ]; then
    log "Analyzing extracellular loop divergence..."
    python3 "${SCRIPTS_DIR}/analyze_ecl.py" \
        --alignments "${RESULTS_DIR}/phylogenies/protein/aligned_*.fasta" \
        --deeptmhmm "${RESULTS_DIR}/deeptmhmm/" \
        --min-ecl-length "${MIN_ECL_LENGTH}" \
        --output "${RESULTS_DIR}/ecl_analysis/ecl_divergence.csv" \
        || log "Warning: ECL divergence analysis failed"
else
    log "Note: DeepTMHMM or phylogeny results not found for ECL analysis"
fi

# Extract all GPCR IDs
awk '/^>/ {print substr($1,2)}' "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" > "${RESULTS_DIR}/ranking/candidate_ids.txt"

# --- RANK_METHOD=rankagg: audit signal independence before ranking ---
# The label-free rank-aggregation reranker must not let signals that share a
# confound (e.g. phylo + og_confidence, both from the same OrthoFinder tree)
# count as several independent votes. The audit groups correlated signals and
# rank_candidates.py fuses each group into one vote. It reads the PRIOR run's
# ranked CSV (the signal-correlation structure is stable across runs); on the
# first run none exists, so rankagg simply treats every signal independently.
# The default RANK_METHOD=weighted path is unchanged (this block is skipped).
RANKED_CSV="${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv"
if [ "${RANK_METHOD:-weighted}" = "rankagg" ]; then
    if [ -f "${RANKED_CSV}" ]; then
        log "RANK_METHOD=rankagg: auditing signal-ranking independence..."
        python3 "${SCRIPTS_DIR}/audit_signal_ranking_independence.py" \
            --ranked-csv "${RANKED_CSV}" \
            --out-prefix "${RESULTS_DIR}/ranking/signal_independence" \
            2>> "${LOGS_DIR}/signal_independence_audit.err" \
            || log --level=WARN "Signal-independence audit failed (rankagg falls back to ungrouped signals)"
    else
        log "RANK_METHOD=rankagg: no prior ranked CSV; rankagg will treat signals as independent (first run)"
    fi
fi

# Run ranking script with improved algorithm
# Note: Now takes synteny directory (not file) for quantitative scoring
run_command "rank_candidates" python3 "${SCRIPTS_DIR}/rank_candidates.py" \
    "${RESULTS_DIR}/ranking/candidate_ids.txt" \
    "${EXPRESSION_DATA}" \
    "${PHYLO_DIR}" \
    "${RESULTS_DIR}/selective_pressure" \
    "${RESULTS_DIR}/synteny" \
    "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv"

# Bead -xqz: HCR-friendliness diagnostic columns (cds_length_bp,
# paralog_min_identity, hcr_probe_friendly). Augments the ranked CSV in place.
HCR_AUG_INPUT="${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv"
if [ -f "$HCR_AUG_INPUT" ]; then
    HCR_ALIGNMENT="${RESULTS_DIR}/phylogenies/protein/class_A/class_A_trimmed.fa"
    [ -f "$HCR_ALIGNMENT" ] || HCR_ALIGNMENT=""
    HCR_CDS="${GENOME_CDS:-}"
    [ -f "$HCR_CDS" ] || HCR_CDS=""
    python3 "${SCRIPTS_DIR}/add_hcr_columns.py" \
        --ranked-csv "$HCR_AUG_INPUT" \
        --cds-fasta "$HCR_CDS" \
        --alignment "$HCR_ALIGNMENT" \
        --out "$HCR_AUG_INPUT" \
        2>> "${LOGS_DIR}/hcr_columns.err" \
        || log --level=WARN "HCR-friendliness column augmentation failed (kept original CSV)"
fi

# Bead -ogc: per-OG ref-CDS coverage transparency columns
# (og_n_ref_cds, og_n_total, og_dnds_reliability). Annotates each row
# with the data-quality context for its dN/dS contribution. Doesn't
# change ranking; lets reviewers see which candidates' rankings depend
# on dN/dS estimates from sparse reference CDS.
COV_REF_CDS="${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna"
COV_OG_TSV=$(find "${RESULTS_DIR}/orthogroups" -name "Orthogroups.tsv" -path "*/Orthogroups/*" 2>/dev/null | head -1)
if [ -f "$HCR_AUG_INPUT" ] && [ -f "$COV_REF_CDS" ] && [ -n "$COV_OG_TSV" ]; then
    python3 "${SCRIPTS_DIR}/add_og_coverage_columns.py" \
        --ranked-csv "$HCR_AUG_INPUT" \
        --cds-fasta "$COV_REF_CDS" \
        --orthogroups-tsv "$COV_OG_TSV" \
        --out "$HCR_AUG_INPUT" \
        2>> "${LOGS_DIR}/og_coverage_columns.err" \
        || log --level=WARN "OG-coverage column augmentation failed (kept original CSV)"
elif [ -f "$HCR_AUG_INPUT" ]; then
    log --level=WARN "Skipping OG-coverage columns: missing $COV_REF_CDS or Orthogroups.tsv"
fi

# Phase 4 / Task 5.1: non-chemoreceptor classification columns
# (classification, classification_confidence, classification_family,
#  classification_subfamily, classification_evidence). Adds columns from
# the 06c consensus TSV; defaults to 'chemoreceptor-candidate' when 06c
# hasn't run or the candidate has no consensus row. Doesn't change rank.
CLASS_TSV="${RESULTS_DIR}/classification/candidate_classifications.tsv"
if [ -f "$HCR_AUG_INPUT" ] && [ -f "$CLASS_TSV" ]; then
    python3 "${SCRIPTS_DIR}/add_classification_columns.py" \
        --ranked-csv "$HCR_AUG_INPUT" \
        --consensus-tsv "$CLASS_TSV" \
        --out "$HCR_AUG_INPUT" \
        2>> "${LOGS_DIR}/classification_columns.err" \
        || log --level=WARN "Classification column augmentation failed (kept original CSV)"
elif [ -f "$HCR_AUG_INPUT" ]; then
    log --level=WARN "Skipping classification columns: $CLASS_TSV not found (run 06c first)"
fi

# --- Weighted-vs-rank-aggregation comparison (always emitted; non-fatal) ---
# Descriptive audit only: reports how far the label-free rank-aggregation order
# would differ from the production weighted order (Spearman, top-k overlap,
# biggest movers) plus an honest, permutation-null positive-control readout.
# The production shortlist is NOT changed here. Uses the signal-independence
# groups.json when the rankagg audit above produced one.
if [ -f "${RANKED_CSV}" ]; then
    COMPARE_ARGS=(--ranked-csv "${RANKED_CSV}" \
        --out "${RESULTS_DIR}/ranking/ranking_method_comparison.md")
    COMPARE_GROUPS="${RESULTS_DIR}/ranking/signal_independence_groups.json"
    [ -f "${COMPARE_GROUPS}" ] && COMPARE_ARGS+=(--groups-json "${COMPARE_GROUPS}")
    [ -f "${HCR_CONTROLS_CSV:-${REFERENCE_DIR}/hcr_positive_controls.csv}" ] && \
        COMPARE_ARGS+=(--controls-csv "${REFERENCE_DIR}/hcr_positive_controls.csv")
    python3 "${SCRIPTS_DIR}/compare_ranking_methods.py" "${COMPARE_ARGS[@]}" \
        2>> "${LOGS_DIR}/ranking_method_comparison.err" \
        || log --level=WARN "Ranking-method comparison failed (non-fatal)"
fi

# Generate plots
python3 "${SCRIPTS_DIR}/plot_ranking.py" "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" "${RESULTS_DIR}/ranking/ranking_plot" || log "Warning: Ranking plot failed"

# Bead -edx: positive-control HCR-validated genes sanity check.
# Non-fatal — pipeline continues regardless of result, but the alert is
# logged + emitted to results/ranking/positive_controls_check.tsv.
HCR_CONTROLS_CSV="${REFERENCE_DIR}/hcr_positive_controls.csv"
if [ -f "$HCR_CONTROLS_CSV" ]; then
    log "Running positive-control HCR sanity check..."
    python3 "${SCRIPTS_DIR}/check_positive_controls.py" \
        --ranked-csv "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" \
        --controls-csv "$HCR_CONTROLS_CSV" \
        --out "${RESULTS_DIR}/ranking/positive_controls_check.tsv" \
        --alert-percentile 50 \
        || log --level=WARN "Positive-control check returned a non-zero exit"
fi

# Create completion flag
touch "${RESULTS_DIR}/step_completed_07.txt"
log "Candidate ranking completed."
