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

# --- RANK_METHOD=rankagg: build evidence channels + audit signal
# independence before ranking (Glue G3) ---
# The label-free rank-aggregation reranker must not let signals that share a
# confound (e.g. phylo + og_confidence, both from the same OrthoFinder tree)
# count as several independent votes. The audit groups correlated signals and
# rank_candidates.py fuses each group into one vote. It reads the PRIOR run's
# ranked CSV (the signal-correlation structure is stable across runs); on the
# first run none exists, so rankagg simply treats every signal independently.
#
# The structural + microswitch evidence-channel PRODUCERS (Glue G1/G5:
# build_structural_channel.py / build_microswitch_channel.py) run here, each
# independently gated on its own raw Unity inputs -- a producer whose inputs
# aren't present is skipped (its channel simply stays dormant downstream, never
# an error). The embedding producer (Glue G2) now runs unconditionally ABOVE
# this block (cw3.6) so its emb_novelty surfaces as an always-present column.
# Outputs land in CHANNELS_DIR, which rank_candidates.py's
# _merge_channels_if_present() reads (default location: alongside the output
# CSV, in a channels/ subdirectory).
#
# The default RANK_METHOD=weighted path is unchanged (this whole block is
# skipped).
RANKED_CSV="${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv"
CHANNELS_DIR="${RESULTS_DIR}/ranking/channels"

# --- Embedding evidence channel PRODUCER (cw3.6: always-on, method-independent) ---
# Writes ${CHANNELS_DIR}/embedding_channel.tsv, consumed by BOTH (a)
# add_embedding_columns.py -> the always-present emb_novelty column in the ranked
# CSV + both views, AND (b) rank_candidates.py's rankagg voter
# (_merge_channels_if_present). It therefore runs REGARDLESS of RANK_METHOD.
# Self-gated on its npz inputs: a scorer whose embeddings aren't present is
# skipped and the channel simply stays dormant downstream (never an error).
mkdir -p "${CHANNELS_DIR}"
EMBEDDINGS_DIR="${EMBEDDINGS_DIR:-${RESULTS_DIR}/ranking/embeddings}"
if [ "${EMB_SCORER:-consensus}" = "consensus" ]; then
    # Locked cw3.6 decision: length-deconfounded proteinclip3b + protrek RRA
    # consensus. Build one tag:cand_npz:ref_npz:identity_tsv spec per model,
    # requiring both npz to exist; run only when >=2 models resolve.
    _emb_specs=(); _emb_missing=0
    for _tag in ${EMB_CONSENSUS_MODELS:-proteinclip3b protrek}; do
        _cand="${EMBEDDINGS_DIR}/candidates_${_tag}_classA.npz"
        _ref="${EMBEDDINGS_DIR}/reference_${_tag}_PROD.npz"
        if [ -f "${_cand}" ] && [ -f "${_ref}" ]; then
            _emb_specs+=("${_tag}:${_cand}:${_ref}:${EMB_CANDIDATE_IDENTITY_TSV}")
        else
            _emb_missing=1
            log "Note: consensus model '${_tag}' npz missing (${_cand} / ${_ref}) -- excluded"
        fi
    done
    if [ "${#_emb_specs[@]}" -ge 2 ] && [ -f "${EMB_CANDIDATE_FASTA}" ]; then
        log "Building consensus embedding channel (RRA, deconfound seq_len; models: ${EMB_CONSENSUS_MODELS:-proteinclip3b protrek})..."
        python3 "${SCRIPTS_DIR}/fusion_consensus.py" \
            --models "${_emb_specs[@]}" \
            --ref-labels "${EMB_REF_LABELS}" \
            --candidate-fasta "${EMB_CANDIDATE_FASTA}" \
            --combiner rra --deconfound seq_len \
            --out "${CHANNELS_DIR}/embedding_channel.tsv" \
            2>> "${LOGS_DIR}/embedding_channel.err" \
            || log --level=WARN "Consensus embedding channel producer failed (channel stays dormant)"
    else
        log "Note: consensus embedding channel skipped (need >=2 model npz + candidate FASTA ${EMB_CANDIDATE_FASTA}) -- stays dormant"
    fi
else
    # Legacy single-model scorer: EMB_SCORER=maha -> ESM-C 600M tied-covariance
    # Mahalanobis; else cosine on ESM-C 300M. Kept for comparison (dominated).
    if [ "${EMB_SCORER}" = "maha" ]; then
        _emb_tag="esmc600m"; _emb_scorer_args="--scorer maha"
    else
        _emb_tag="esmc300m"; _emb_scorer_args=""
    fi
    ESMC_CANDIDATE_NPZ="${ESMC_CANDIDATE_NPZ:-${EMBEDDINGS_DIR}/candidates_${_emb_tag}.npz}"
    ESMC_REFERENCE_NPZ="${ESMC_REFERENCE_NPZ:-${EMBEDDINGS_DIR}/reference_${_emb_tag}.npz}"
    if [ -f "${ESMC_CANDIDATE_NPZ}" ] && [ -f "${ESMC_REFERENCE_NPZ}" ]; then
        log "Building ESM-C embedding channel (scorer=${EMB_SCORER})..."
        python3 "${SCRIPTS_DIR}/build_embedding_channel.py" \
            --candidate-npz "${ESMC_CANDIDATE_NPZ}" \
            --ref-npz "${ESMC_REFERENCE_NPZ}" \
            --ref-labels "${REFERENCE_DIR}/anchors/anchor_set.tsv" \
            ${_emb_scorer_args} \
            --out "${CHANNELS_DIR}/embedding_channel.tsv" \
            2>> "${LOGS_DIR}/build_embedding_channel.err" \
            || log --level=WARN "Embedding channel producer failed (channel stays dormant)"
    else
        log "Note: ESM-C embeddings not found (${ESMC_CANDIDATE_NPZ}, ${ESMC_REFERENCE_NPZ}) -- embedding channel skipped, stays dormant"
    fi
fi

if [ "${RANK_METHOD:-weighted}" = "rankagg" ]; then
    mkdir -p "${CHANNELS_DIR}"

    # Structural channel (Foldseek vs PDB/AFDB50/GPCRdb; Glue G1). Gated on
    # at least one of the 3 Foldseek DB tabs scripts/unity/run_foldseek_candidates.sh
    # produces -- that script itself tolerates a subset of its 3 DBs being
    # unavailable, so we mirror that tolerance here rather than requiring all 3.
    FOLDSEEK_CANDIDATES_DIR="${FOLDSEEK_CANDIDATES_DIR:-${RESULTS_DIR}/foldseek/candidates}"
    FOLDSEEK_TSVS=()
    for db in PDB AFDB50 GPCRdb; do
        [ -f "${FOLDSEEK_CANDIDATES_DIR}/${db}.tsv" ] && FOLDSEEK_TSVS+=("${FOLDSEEK_CANDIDATES_DIR}/${db}.tsv")
    done
    if [ "${#FOLDSEEK_TSVS[@]}" -gt 0 ]; then
        log "RANK_METHOD=rankagg: building structural evidence channel (${#FOLDSEEK_TSVS[@]} Foldseek DB tab(s))..."
        python3 "${SCRIPTS_DIR}/build_structural_channel.py" \
            --foldseek-tsvs "${FOLDSEEK_TSVS[@]}" \
            --anchor-set "${REFERENCE_DIR}/anchors/anchor_set.tsv" \
            --out "${CHANNELS_DIR}/structural_channel.tsv" \
            2>> "${LOGS_DIR}/build_structural_channel.err" \
            || log --level=WARN "Structural channel producer failed (channel stays dormant)"
    else
        log "Note: no Foldseek candidate hit TSVs found under ${FOLDSEEK_CANDIDATES_DIR} (run scripts/unity/run_foldseek_candidates.sh) -- structural channel skipped, stays dormant"
    fi

    # (The embedding channel producer now runs unconditionally above, so its
    # emb_novelty can surface as an always-present column on the weighted path
    # too — cw3.6. It still writes ${CHANNELS_DIR}/embedding_channel.tsv, which
    # this block's rank_candidates.py rankagg voter reads as before.)

    # OR-microswitch channel (structural BW number-transfer; Glue G5). Gated
    # on the per-candidate alignment-report directory existing and non-empty.
    MICROSWITCH_ALIGNMENTS_DIR="${MICROSWITCH_ALIGNMENTS_DIR:-${RESULTS_DIR}/ranking/microswitch/alignments}"
    if [ -d "${MICROSWITCH_ALIGNMENTS_DIR}" ] && [ -n "$(ls -A "${MICROSWITCH_ALIGNMENTS_DIR}" 2>/dev/null)" ]; then
        log "RANK_METHOD=rankagg: building OR-microswitch evidence channel..."
        python3 "${SCRIPTS_DIR}/build_microswitch_channel.py" \
            --alignments "${MICROSWITCH_ALIGNMENTS_DIR}" \
            --out "${CHANNELS_DIR}/microswitch_channel.tsv" \
            2>> "${LOGS_DIR}/build_microswitch_channel.err" \
            || log --level=WARN "Microswitch channel producer failed (channel stays dormant)"
    else
        log "Note: no microswitch alignment reports found under ${MICROSWITCH_ALIGNMENTS_DIR} (run scripts/unity/run_microswitch_bw_transfer.sh) -- microswitch channel skipped, stays dormant"
    fi

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

# cw3.6: always-present emb_novelty column. Left-joins the embedding channel's
# emb_novelty (+ emb_nonchemo_family, has_emb_data) into the ranked CSV as a
# sortable descriptive column in BOTH views, independent of RANK_METHOD (the
# channel otherwise only feeds the rankagg voter). Dormant channel -> the column
# is present but empty. Must run BEFORE emit_ranked_views (its discovery score
# reads emb_novelty). The channel is authoritative (overwrites any rankagg-merged
# copy), so running on either RANK_METHOD path yields exactly one emb_novelty.
if [ -f "$HCR_AUG_INPUT" ]; then
    python3 "${SCRIPTS_DIR}/add_embedding_columns.py" \
        --ranked-csv "$HCR_AUG_INPUT" \
        --channel-tsv "${CHANNELS_DIR}/embedding_channel.tsv" \
        --out "$HCR_AUG_INPUT" \
        2>> "${LOGS_DIR}/embedding_columns.err" \
        || log --level=WARN "Embedding column augmentation failed (kept original CSV)"
fi

# Bead 1nr: two ranked views. Re-project the composite-sorted, fully-augmented
# ranked CSV into a CONFIDENCE shortlist (safe-bet chemoreceptor candidates with
# complete evidence) and a DISCOVERY view (high-novelty divergent-LSE candidates
# — reference-poor / manual-review / high tandem signal — that a single
# composite sort would bury). Emits two sibling CSVs; doesn't touch the ranked
# CSV. Must run AFTER add_classification_columns.py (needs classification +
# needs_manual_review) and add_og_coverage_columns.py (needs og_dnds_reliability).
if [ -f "$HCR_AUG_INPUT" ]; then
    python3 "${SCRIPTS_DIR}/emit_ranked_views.py" \
        --ranked-csv "$HCR_AUG_INPUT" \
        --confidence-out "${RESULTS_DIR}/ranking/ranked_candidates_confidence.csv" \
        --discovery-out "${RESULTS_DIR}/ranking/ranked_candidates_discovery.csv" \
        2>> "${LOGS_DIR}/ranked_views.err" \
        || log --level=WARN "Two-ranked-views emission failed (kept ranked CSV)"
fi

# --- Weighted-vs-rank-aggregation comparison (always emitted; non-fatal) ---
# Descriptive audit only: reports how far the hand-weighted-sum order and the
# label-free rank-aggregation order differ (Spearman, top-k overlap, biggest
# movers) plus an honest, permutation-null positive-control readout. Both orders
# are recomputed from the ranked CSV's signal columns, so this stays a valid
# audit whichever RANK_METHOD is the production default. It changes nothing.
# Uses the signal-independence groups.json when the rankagg audit produced one.
if [ -f "${RANKED_CSV}" ]; then
    COMPARE_ARGS=(--ranked-csv "${RANKED_CSV}" \
        --out "${RESULTS_DIR}/ranking/ranking_method_comparison.md")
    COMPARE_GROUPS="${RESULTS_DIR}/ranking/signal_independence_groups.json"
    [ -f "${COMPARE_GROUPS}" ] && COMPARE_ARGS+=(--groups-json "${COMPARE_GROUPS}")
    [ -f "${HCR_CONTROLS_CSV:-${REFERENCE_DIR}/hcr_positive_controls.csv}" ] && \
        COMPARE_ARGS+=(--controls-csv "${HCR_CONTROLS_CSV:-${REFERENCE_DIR}/hcr_positive_controls.csv}")
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
