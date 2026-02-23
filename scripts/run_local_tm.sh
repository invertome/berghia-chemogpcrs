#!/bin/bash
# run_local_tm.sh - Wrapper for predict_tm_helices.py matching DeepTMHMM interface
# Usage: run_local_tm.sh -f input.fasta -o output_dir
# Produces output_dir/prediction in the same format as DeepTMHMM

set -euo pipefail

FASTA=""
OUTDIR="tm_output"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--fasta) FASTA="$2"; shift 2 ;;
        -o|-od|--output-dir) OUTDIR="$2"; shift 2 ;;
        *) shift ;;
    esac
done

if [ -z "$FASTA" ] || [ ! -f "$FASTA" ]; then
    echo "Error: FASTA file not found: $FASTA" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$OUTDIR"

echo "Running local TM helix prediction on $(grep -c '^>' "$FASTA") sequences..." >&2

python3 "${SCRIPT_DIR}/predict_tm_helices.py" "$FASTA" "${OUTDIR}/prediction"

echo "TM prediction completed: ${OUTDIR}/prediction" >&2
