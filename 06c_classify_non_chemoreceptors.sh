#!/bin/bash
#SBATCH --job-name=berghia_06c_classify
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# 06c_classify_non_chemoreceptors.sh — Stage 06c of the Berghia
# chemoreceptor pipeline. Phase 4 Task 4.5 of the non-chemoreceptor
# classification feature.
#
# Wraps three independent classifiers + their consensus into a single
# stage:
#   1. HMM scan          (scripts/classify_via_hmm.py)
#   2. Orthogroup vote   (scripts/classify_via_og_vote.py)
#   3. Phylogenetic placement (scripts/classify_via_placement.py)
#   4. Consensus         (scripts/classify_consensus.py)
#
# Inputs:
#   - Berghia candidate FASTA (stage 02b output)
#   - Custom HMM library + LOO thresholds (Phase 2: Tasks 2.1, 2.2)
#   - Pfam fallback HMMs (Phase 2 Task 2.3)
#   - OrthoFinder Orthogroups.tsv (stage 03)
#   - Reference annotations (Phase 1: references/non_chemo_gpcr/all_references.tsv)
#   - Reference trees + alignments + leaf TSVs (Phase 3)
#
# Output:
#   ${RESULTS_DIR}/classification/candidate_classifications.tsv
# (consumed by stage 07's add_classification_columns.py)
#
# Idempotent: skips already-completed sub-steps via step_completed_*.txt
# markers consistent with the rest of the pipeline.

set -eo pipefail
mkdir -p logs

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# shellcheck disable=SC1091
source config.sh
# shellcheck disable=SC1091
source functions.sh
# ONE deterministic, chronologically-correct rule for which OrthoFinder
# run is authoritative (mtime of Orthogroups.tsv). Shared by stages
# 03/03b/04/05/06c/07 so they can no longer resolve different runs.
# shellcheck source=scripts/orthofinder_paths.sh
source "${SCRIPTS_DIR:-scripts}/orthofinder_paths.sh"

THREADS="${SLURM_CPUS_PER_TASK:-${CPUS:-4}}"

CANDIDATE_FASTA="${CANDIDATE_FASTA:-${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa}"
HMM_DIR="${RESULTS_DIR}/classification/hmms"
PFAM_DIR="${RESULTS_DIR}/classification/hmms/pfam_fallback"
LOO_METRICS="${RESULTS_DIR}/classification/loo/loo_metrics.tsv"
TREE_DIR="${RESULTS_DIR}/classification/trees"
ANNOTATIONS_TSV="${REFERENCE_DIR:-references}/non_chemo_gpcr/all_references.tsv"
ORTHOGROUPS_TSV=$(resolve_orthogroups_tsv "${RESULTS_DIR}/orthogroups" || true)

OUT_DIR="${RESULTS_DIR}/classification"
mkdir -p "$OUT_DIR"

log "=== Stage 06c: non-chemoreceptor classification ==="

# Prerequisites
check_file "$CANDIDATE_FASTA"        || { log "Error: missing candidate FASTA: $CANDIDATE_FASTA"; exit 1; }
check_file "$ANNOTATIONS_TSV"        || { log "Error: missing reference annotations: $ANNOTATIONS_TSV"; exit 1; }
check_dir "$HMM_DIR"                 || { log "Error: HMM library missing — build with build_classification_hmms.sh"; exit 1; }
[ -n "$ORTHOGROUPS_TSV" ] && [ -f "$ORTHOGROUPS_TSV" ] || {
    log --level=WARN "Orthogroups.tsv not found; OG-vote will be skipped (all candidates 'unclassified-og')"
    ORTHOGROUPS_TSV=""
}
# Orthology quarantine (ORTHOLOGY_SOURCE_TRUSTED, see config.sh). The OG-vote
# classifier is entirely Nath-fed: it votes within Nath-built orthogroups over
# Nath references annotated by step 1b below. Under the quarantine it must not
# vote at all, so reuse the existing, tested "no Orthogroups.tsv" branch rather
# than adding a second way to be absent.
#
# CONSEQUENCE, stated plainly because it is significant: classify_consensus.py
# reaches 'non-chemoreceptor' only on 3-of-3 agreement and
# 'likely-non-chemoreceptor' only when HMM *and* OG agree — placement is
# explicitly barred from substituting for either. With the OG source dark,
# NEITHER demotion tier is reachable, so every candidate stays
# 'chemoreceptor-candidate'. The classifier is conservative-by-default, which is
# the safe direction (it cannot demote a real chemoreceptor), but it means the
# non-chemoreceptor filter provides NO suppression until orthology is rebuilt.
if [ "${ORTHOLOGY_SOURCE_TRUSTED:-0}" != "1" ] && [ -n "$ORTHOGROUPS_TSV" ]; then
    log --level=WARN "Quarantine: OG-vote classifier disabled (orthogroups are Nath-derived and not trusted). Consensus falls to HMM + placement, which cannot reach either demotion tier -- every candidate will stay 'chemoreceptor-candidate'."
    ORTHOGROUPS_TSV=""
fi

N_CANDIDATES=$(grep -c '^>' "$CANDIDATE_FASTA")
log "Candidates: $N_CANDIDATES sequences in $CANDIDATE_FASTA"

# --- 1a. HMM scan classifier on Berghia candidates ---
HMM_OUT="${OUT_DIR}/candidate_hmm_assignments.tsv"
# Audit (-wkb): PFAM_ARGS/LOO_ARGS feed BOTH the candidate HMM scan and the
# Nath-ref scan below; compute them unconditionally so an idempotent resume
# (06c_classify_hmm already done but 06c_classify_nath not) still passes
# --pfam-fallback-dir / --loo-metrics instead of the weaker no-fallback scan.
PFAM_ARGS=""
if [ -d "$PFAM_DIR" ]; then
    PFAM_ARGS="--pfam-fallback-dir $PFAM_DIR"
fi
LOO_ARGS=""
if [ -f "$LOO_METRICS" ]; then
    LOO_ARGS="--loo-metrics $LOO_METRICS"
fi

if ! is_step_completed "06c_classify_hmm"; then
    log "[06c] Running HMM-scan classifier on candidates..."
    # shellcheck disable=SC2086
    python3 "${SCRIPTS_DIR}/classify_via_hmm.py" \
        --candidate-fasta "$CANDIDATE_FASTA" \
        --hmm-dir "$HMM_DIR" \
        $PFAM_ARGS $LOO_ARGS \
        --output-tsv "$HMM_OUT" \
        --threads "$THREADS"
    create_checkpoint "06c_classify_hmm"
else
    log "[06c] HMM-scan classifier (candidates): already completed"
fi

# --- 1b. HMM scan on the pipeline's rewritten reference set ---
# This produces a `ref_TAXID_N -> family` annotation TSV that matches
# OrthoFinder member IDs (post-update_headers.py rewrite). Without it,
# OG-vote can't find any annotated members because Swiss-Prot accessions
# don't match the rewritten OrthoFinder IDs (H6 from 2026-05-07 review).
NATH_REWRITTEN="${RESULTS_DIR}/reference_sequences/all_references.fa"
NATH_CLASSIFICATIONS="${OUT_DIR}/nath_classifications.tsv"
if ! is_step_completed "06c_classify_nath"; then
    if [ -f "$NATH_REWRITTEN" ]; then
        log "[06c] Running HMM-scan classifier on Nath et al rewritten refs..."
        # shellcheck disable=SC2086
        python3 "${SCRIPTS_DIR}/classify_via_hmm.py" \
            --candidate-fasta "$NATH_REWRITTEN" \
            --hmm-dir "$HMM_DIR" \
            $PFAM_ARGS $LOO_ARGS \
            --output-tsv "$NATH_CLASSIFICATIONS" \
            --threads "$THREADS"
        create_checkpoint "06c_classify_nath"
    else
        log --level=WARN "Skipping Nath classification: $NATH_REWRITTEN not found. OG-vote will rely on Swiss-Prot annotations only (likely zero matches)."
        NATH_CLASSIFICATIONS=""
    fi
else
    log "[06c] HMM-scan classifier (Nath refs): already completed"
fi

# --- 2. OG-vote classifier ---
OG_OUT="${OUT_DIR}/candidate_og_assignments.tsv"
if ! is_step_completed "06c_classify_og"; then
    if [ -n "$ORTHOGROUPS_TSV" ]; then
        log "[06c] Running OG-vote classifier..."
        NATH_ARG=""
        if [ -n "$NATH_CLASSIFICATIONS" ] && [ -f "$NATH_CLASSIFICATIONS" ]; then
            NATH_ARG="--hmm-classifications-tsv $NATH_CLASSIFICATIONS"
        fi
        # shellcheck disable=SC2086
        python3 "${SCRIPTS_DIR}/classify_via_og_vote.py" \
            --candidate-fasta "$CANDIDATE_FASTA" \
            --orthogroups-tsv "$ORTHOGROUPS_TSV" \
            --annotations-tsv "$ANNOTATIONS_TSV" \
            $NATH_ARG \
            --output-tsv "$OG_OUT"
    else
        log --level=WARN "Skipping OG-vote (no Orthogroups.tsv); writing empty output"
        echo -e "candidate_id\tog_id\tog_vote_family\tog_vote_subfamily\tn_members\tn_annotated\tconsensus_fraction" > "$OG_OUT"
    fi
    create_checkpoint "06c_classify_og"
else
    log "[06c] OG-vote classifier: already completed"
fi

# --- 3. Placement classifier (optional — needs reference trees) ---
PLACEMENT_OUT="${OUT_DIR}/candidate_placement.tsv"
if ! is_step_completed "06c_classify_placement"; then
    if [ -d "$TREE_DIR" ] && [ -f "$TREE_DIR/backbone.treefile" ]; then
        log "[06c] Running phylogenetic placement classifier..."
        python3 "${SCRIPTS_DIR}/classify_via_placement.py" \
            --candidate-fasta "$CANDIDATE_FASTA" \
            --tree-dir "$TREE_DIR" \
            --output-tsv "$PLACEMENT_OUT" \
            --threads "$THREADS" \
            || log --level=WARN "Placement classifier failed; consensus will fall back to HMM+OG (medium confidence cap)"
    else
        log --level=WARN "Reference trees not found at $TREE_DIR — skipping placement; consensus will use HMM+OG only (medium confidence cap)"
        PLACEMENT_OUT=""
    fi
    create_checkpoint "06c_classify_placement"
else
    log "[06c] Placement classifier: already completed"
fi

# --- 4. Consensus ---
CONSENSUS_OUT="${OUT_DIR}/candidate_classifications.tsv"
if ! is_step_completed "06c_classify_consensus"; then
    log "[06c] Running 3-source consensus..."
    PLACE_ARG=""
    if [ -n "$PLACEMENT_OUT" ] && [ -f "$PLACEMENT_OUT" ]; then
        PLACE_ARG="--placement-tsv $PLACEMENT_OUT"
    fi
    # shellcheck disable=SC2086
    python3 "${SCRIPTS_DIR}/classify_consensus.py" \
        --hmm-tsv "$HMM_OUT" \
        --og-tsv "$OG_OUT" \
        $PLACE_ARG \
        --out "$CONSENSUS_OUT"
    create_checkpoint "06c_classify_consensus"
else
    log "[06c] Consensus: already completed"
fi

# --- Summary ---
if [ -f "$CONSENSUS_OUT" ]; then
    n_total=$(awk -F'\t' 'NR>1' "$CONSENSUS_OUT" | wc -l)
    n_high=$(awk -F'\t' 'NR>1 && $3=="high"' "$CONSENSUS_OUT" | wc -l)
    n_med=$(awk -F'\t' 'NR>1 && $3=="medium"' "$CONSENSUS_OUT" | wc -l)
    n_chem=$(awk -F'\t' 'NR>1 && $2=="chemoreceptor-candidate"' "$CONSENSUS_OUT" | wc -l)
    log "Stage 06c summary:"
    log "  candidates total:                  $n_total"
    log "  non-chemoreceptor (high conf):     $n_high"
    log "  likely-non-chemo (medium conf):    $n_med"
    log "  chemoreceptor-candidate (default): $n_chem"
fi

touch "${RESULTS_DIR}/step_completed_06c.txt"
log "=== Stage 06c complete ==="
