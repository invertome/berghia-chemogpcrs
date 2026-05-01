#!/bin/bash
# run_tandem_detection.sh — Run intra-genome tandem-cluster detection on the
# Berghia GPCR candidates (bead -ar8). Idempotent.
#
# Reads paths from config.sh (sourced by the calling pipeline stage):
#   GENOME_GFF              — Berghia GFF3 (default at config.sh)
#   TANDEM_WINDOW_KB        — sliding window (default 100)
#   TANDEM_MIN_SIZE         — minimum cluster size (default 3)
#   TANDEM_CLUSTERS_FILE    — output CSV path (default: $RESULTS_DIR/synteny/tandem_clusters.csv)
#
# Usage:
#   bash scripts/run_tandem_detection.sh <candidates_file>
# where candidates_file is the FASTA or .txt produced by step 02 / 02a.

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <candidates_file>" >&2
    exit 1
fi

CANDIDATES="$1"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Source config if not already
if [ -z "${GENOME_GFF:-}" ]; then
    # shellcheck disable=SC1091
    source "$BASE_DIR/config.sh"
fi

if [ ! -e "$GENOME_GFF" ]; then
    echo "[run_tandem_detection] GFF not found at $GENOME_GFF" >&2
    echo "                      Run scripts/fetch_berghia_genome.sh first." >&2
    exit 2
fi

if [ ! -e "$CANDIDATES" ]; then
    echo "[run_tandem_detection] Candidates file not found: $CANDIDATES" >&2
    exit 3
fi

OUT="${TANDEM_CLUSTERS_FILE:-${RESULTS_DIR:-results}/synteny/tandem_clusters.csv}"
mkdir -p "$(dirname "$OUT")"

python3 "$SCRIPT_DIR/compute_tandem_clusters.py" \
    --gff "$GENOME_GFF" \
    --candidates "$CANDIDATES" \
    --window-kb "${TANDEM_WINDOW_KB:-100}" \
    --min-size "${TANDEM_MIN_SIZE:-3}" \
    --out "$OUT"

echo "[run_tandem_detection] OK -> $OUT" >&2
