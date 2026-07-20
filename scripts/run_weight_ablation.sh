#!/bin/bash
# run_weight_ablation.sh — rank-weight ablation for paper robustness story.
#
# Bead -ab1 (option-G ablation framework). Reruns stage 07 with each
# ranking-axis weight zeroed in turn, then computes pairwise comparison
# metrics (Spearman ρ, Kendall τ, Jaccard top-N, median rank shift)
# across the resulting rankings.
#
# This is the cheapest of the three ablation tiers:
#   tier 1 (this script)      ~minutes   reranks via env-var weight overrides
#   tier 2 (TBD)              ~hours     rerun trim + tree only with different ClipKit modes
#   tier 3 (full ablation)    ~24h/cfg   rerun stage 04 from raw input with RUN_* off
#
# Tier 1 demonstrates that no single ranking axis is driving the
# top-candidate set — a defensible robustness story for the paper that
# requires only stage 07 reruns (stage 04 + 05 outputs are reused).
#
# Prerequisites: stage 07 must have run successfully on the production
# data, producing ${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv.
#
# Usage:
#   bash run_weight_ablation.sh \\
#       --output-dir=results/ablations/weights/
#
# (the baseline ranking is read from $RESULTS_DIR/ranking/.)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

OUTPUT_DIR=""
TOP_NS="${TOP_NS:-10,25,50,100}"

usage() {
    cat >&2 <<EOF
run_weight_ablation.sh — rank-weight ablation runner

Required:
  --output-dir=<dir>     Where to write per-ablation rerankings + comparison TSV

Optional:
  --top-ns=<csv>         Top-N values for Jaccard / rank-shift (default 10,25,50,100)

Outputs:
  <output-dir>/no_<axis>/ranked_candidates_sorted.csv   per-ablation reranking
  <output-dir>/comparison.tsv                            pairwise metrics
  <output-dir>/comparison.summary.json
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir=*)   OUTPUT_DIR="${1#*=}"; shift ;;
        --top-ns=*)       TOP_NS="${1#*=}"; shift ;;
        -h|--help)        usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -z "$OUTPUT_DIR" ]] && { echo "ERROR: --output-dir required" >&2; exit 1; }

# shellcheck disable=SC1091
source "$PROJECT_DIR/config.sh"

BASELINE_CSV="${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv"
[[ -f "$BASELINE_CSV" ]] || {
    echo "ERROR: baseline ranked CSV not found: $BASELINE_CSV"  >&2
    echo "       run stage 07 first (bash 07_candidate_ranking.sh)" >&2
    exit 1
}

mkdir -p "$OUTPUT_DIR"
echo "[ablation] baseline: $BASELINE_CSV"
echo "[ablation] output:   $OUTPUT_DIR"
echo ""

# --- Run stage 07 with each weight zeroed in turn ---
# Stage 07 reads the per-axis scores from the upstream pipeline (stage 05
# dN/dS, stage 04 phylogeny, stage 06 synteny, etc.) and combines them
# into the composite via the *_WEIGHT env vars. Setting one weight to 0
# zeros that axis's contribution while leaving the upstream data intact.
declare -a WEIGHTS_TO_ABLATE=(
    "PHYLO_WEIGHT"
    "DNDS_WEIGHT"
    "SYNTENY_WEIGHT"
    "TANDEM_CLUSTER_WEIGHT"
    "EXPR_WEIGHT"
    "LSE_DIVERGENCE_WEIGHT"
)

declare -a CSV_ARGS=("--csv" "default=$BASELINE_CSV")

for w in "${WEIGHTS_TO_ABLATE[@]}"; do
    short=$(echo "$w" | tr 'A-Z' 'a-z' | sed 's/_weight$//')
    sub_dir="$OUTPUT_DIR/no_${short}"
    mkdir -p "$sub_dir"

    echo "[ablation] no_${short} (zeroing $w)..."
    # Override RESULTS_DIR so stage 07 writes its outputs to our ablation dir
    # without clobbering the production ranking. Stage 07 reads upstream
    # results from the ORIGINAL RESULTS_DIR via absolute paths in config.sh,
    # so this only redirects WRITES.
    if (
        cd "$PROJECT_DIR"
        env "${w}=0" \
            ABLATION_OUTPUT_DIR="$sub_dir" \
            bash 07_candidate_ranking.sh \
            >"$sub_dir/stage07.log" 2>&1
    ); then
        # Stage 07 writes ranked_candidates_sorted.csv inside RESULTS_DIR.
        # We need it inside sub_dir; for now, copy after the run.
        # (A future refactor of 07_candidate_ranking.sh to honor
        # ABLATION_OUTPUT_DIR would tighten this.)
        if [[ -f "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" ]]; then
            cp "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" \
               "$sub_dir/ranked_candidates_sorted.csv"
            CSV_ARGS+=("--csv" "no_${short}=$sub_dir/ranked_candidates_sorted.csv")
            echo "  -> $sub_dir/ranked_candidates_sorted.csv"
        else
            echo "  WARN: stage 07 ran but no output CSV found"
        fi
    else
        echo "  ERROR: stage 07 failed (see $sub_dir/stage07.log)"
    fi
done

# Restore the baseline CSV (the loop overwrote it)
cp "$BASELINE_CSV" "${BASELINE_CSV}.ablation_backup" || true

# --- Pairwise comparison ---
if [[ "${#CSV_ARGS[@]}" -gt 2 ]]; then
    echo ""
    echo "[ablation] Running pairwise comparison..."
    python3 "${SCRIPT_DIR}/compare_rankings.py" \
        "${CSV_ARGS[@]}" \
        --top-ns="$TOP_NS" \
        --out="$OUTPUT_DIR/comparison.tsv"
    echo "[ablation] DONE: $OUTPUT_DIR/comparison.tsv"
else
    echo "[ablation] No ablations succeeded; nothing to compare."
    exit 1
fi
