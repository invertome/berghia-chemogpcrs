#!/bin/bash
# run_deeptmhmm.sh - Unified TM prediction wrapper
# Tries: (1) local deeptmhmm binary, (2) biolib CLI, (3) local Kyte-Doolittle fallback
# Output: DeepTMHMM-compatible prediction file
#
# Usage: bash run_deeptmhmm.sh -f INPUT.fasta -o OUTPUT_DIR

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FASTA=""
OUTPUT_DIR=""

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        -f) FASTA="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        *) shift ;;
    esac
done

if [ -z "$FASTA" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 -f INPUT.fasta -o OUTPUT_DIR" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Strategy 1: Local DeepTMHMM binary
if command -v deeptmhmm &>/dev/null; then
    echo "Using local DeepTMHMM binary" >&2
    deeptmhmm -f "$FASTA" -o "$OUTPUT_DIR"
    exit $?
fi

# Strategy 2: biolib CLI
if command -v biolib &>/dev/null; then
    echo "Using biolib DeepTMHMM" >&2
    BIOLIB_OUT="${OUTPUT_DIR}/biolib_output"
    mkdir -p "$BIOLIB_OUT"
    # biolib outputs to CWD, so cd into output dir
    (cd "$BIOLIB_OUT" && biolib run DTU/DeepTMHMM --fasta "$FASTA") 2>/dev/null
    if [ -f "$BIOLIB_OUT/predicted_topologies.3line" ]; then
        # Parse biolib 3-line format into simple prediction TSV
        python3 -c "
import sys
out_dir = '${OUTPUT_DIR}'
with open('${BIOLIB_OUT}/predicted_topologies.3line') as f:
    lines = f.readlines()
with open(f'{out_dir}/prediction', 'w') as out:
    i = 0
    while i < len(lines):
        header = lines[i].strip().lstrip('>')
        topology = lines[i+1].strip() if i+1 < len(lines) else ''
        # Count TM regions (consecutive 'M' stretches in topology string)
        in_tm = False
        tm_count = 0
        for c in topology:
            if c in ('H', 'h', 'M') and not in_tm:
                tm_count += 1
                in_tm = True
            elif c not in ('H', 'h', 'M'):
                in_tm = False
        pred_type = 'TM' if tm_count > 0 else 'GLOB'
        out.write(f'{header}\t{pred_type}\t1.0\t\t{tm_count}\n')
        i += 3  # Skip sequence line
"
        if [ -f "$OUTPUT_DIR/prediction" ]; then
            echo "biolib DeepTMHMM completed successfully" >&2
            exit 0
        fi
    fi
    echo "Warning: biolib DeepTMHMM output not in expected format, falling back" >&2
fi

# Strategy 3: Local Kyte-Doolittle predictor (fast fallback)
echo "WARNING: Using local Kyte-Doolittle TM predictor (less accurate than DeepTMHMM)" >&2
echo "Install DeepTMHMM or biolib for production accuracy" >&2
python3 "${SCRIPT_DIR}/predict_tm_helices.py" "$FASTA" "$OUTPUT_DIR/prediction"
exit $?
