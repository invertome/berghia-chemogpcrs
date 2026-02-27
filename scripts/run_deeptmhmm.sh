#!/bin/bash
# run_deeptmhmm.sh - Unified TM prediction wrapper
# Tries (in order):
#   1. Apptainer SIF container (DEEPTMHMM_SIF or containers/deeptmhmm.sif)
#   2. Local venv installation (tools/deeptmhmm with .venv)
#   3. System deeptmhmm binary
#   4. biolib CLI
#   5. Kyte-Doolittle fallback (with topology output)
#
# Output: DeepTMHMM-compatible files:
#   prediction                      - TSV: seqid, pred_type, confidence, SP, n_tm_regions
#   predicted_topologies.3line      - Per-residue topology per sequence
#                                     Format: >header\nsequence\ntopology_labels\n
#
# Usage: bash run_deeptmhmm.sh -f INPUT.fasta -o OUTPUT_DIR

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
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

# Resolve to absolute path for container bind mounts
FASTA=$(realpath "$FASTA")
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

# --- Helper: Convert 3-line topology file to prediction TSV ---
# DeepTMHMM 3-line format:
#   >header [| TYPE]
#   amino_acid_sequence
#   topology_labels (I/O/M/S/B/P characters)
generate_prediction_tsv() {
    local threeline_file="$1"
    local output_file="$2"
    python3 -c "
import sys
with open('${threeline_file}') as f:
    lines = f.readlines()
with open('${output_file}', 'w') as out:
    i = 0
    while i < len(lines):
        if not lines[i].startswith('>'):
            i += 1
            continue
        # Line 1: >header [| TYPE]
        header = lines[i].strip().lstrip('>').split()[0]
        # Line 2: amino acid sequence (skip)
        # Line 3: topology labels — count consecutive M stretches for TM regions
        topology = lines[i+2].strip() if i+2 < len(lines) else ''
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
        i += 3
"
}

# --- Strategy 1: Apptainer SIF container ---
SIF="${DEEPTMHMM_SIF:-$BASE_DIR/containers/deeptmhmm.sif}"
if [ -f "$SIF" ]; then
    echo "Using Apptainer DeepTMHMM container: $SIF" >&2
    # DeepTMHMM's predict.py errors if output dir exists, use a temp subdir
    DTMHMM_OUT="${OUTPUT_DIR}/dtmhmm_run"
    rm -rf "$DTMHMM_OUT"

    apptainer run --bind "$(dirname "$FASTA")":"$(dirname "$FASTA"):ro" \
                  --bind "$OUTPUT_DIR":"$OUTPUT_DIR" \
                  "$SIF" --fasta "$FASTA" --output-dir "$DTMHMM_OUT"
    EXITCODE=$?

    if [ $EXITCODE -eq 0 ] && [ -f "$DTMHMM_OUT/predicted_topologies.3line" ]; then
        cp "$DTMHMM_OUT/predicted_topologies.3line" "$OUTPUT_DIR/predicted_topologies.3line"
        [ -f "$DTMHMM_OUT/TMRs.gff3" ] && cp "$DTMHMM_OUT/TMRs.gff3" "$OUTPUT_DIR/TMRs.gff3"
        generate_prediction_tsv "$OUTPUT_DIR/predicted_topologies.3line" "$OUTPUT_DIR/prediction"
        rm -rf "$DTMHMM_OUT"
        echo "Apptainer DeepTMHMM completed successfully" >&2
        exit 0
    else
        echo "Warning: Apptainer DeepTMHMM failed (exit=$EXITCODE), falling back" >&2
        rm -rf "$DTMHMM_OUT"
    fi
fi

# --- Strategy 2: Local DeepTMHMM installation (tools/deeptmhmm with venv) ---
DTMHMM_INSTALL="${DEEPTMHMM_INSTALL:-$BASE_DIR/tools/deeptmhmm}"
if [ -f "$DTMHMM_INSTALL/predict.py" ] && [ -f "$DTMHMM_INSTALL/.venv/bin/python3" ]; then
    echo "Using local DeepTMHMM installation: $DTMHMM_INSTALL" >&2
    DTMHMM_OUT="${OUTPUT_DIR}/dtmhmm_run"
    rm -rf "$DTMHMM_OUT"

    (cd "$DTMHMM_INSTALL" && "$DTMHMM_INSTALL/.venv/bin/python3" predict.py \
        --fasta "$FASTA" --output-dir "$DTMHMM_OUT")
    EXITCODE=$?

    if [ $EXITCODE -eq 0 ] && [ -f "$DTMHMM_OUT/predicted_topologies.3line" ]; then
        cp "$DTMHMM_OUT/predicted_topologies.3line" "$OUTPUT_DIR/predicted_topologies.3line"
        [ -f "$DTMHMM_OUT/TMRs.gff3" ] && cp "$DTMHMM_OUT/TMRs.gff3" "$OUTPUT_DIR/TMRs.gff3"
        generate_prediction_tsv "$OUTPUT_DIR/predicted_topologies.3line" "$OUTPUT_DIR/prediction"
        rm -rf "$DTMHMM_OUT"
        echo "Local DeepTMHMM installation completed successfully" >&2
        exit 0
    else
        echo "Warning: Local DeepTMHMM installation failed (exit=$EXITCODE), falling back" >&2
        rm -rf "$DTMHMM_OUT"
    fi
fi

# --- Strategy 3: System DeepTMHMM binary ---
if command -v deeptmhmm &>/dev/null; then
    echo "Using local DeepTMHMM binary" >&2
    DTMHMM_OUT="${OUTPUT_DIR}/dtmhmm_run"
    rm -rf "$DTMHMM_OUT"

    deeptmhmm -f "$FASTA" -o "$DTMHMM_OUT"
    EXITCODE=$?

    if [ $EXITCODE -eq 0 ] && [ -d "$DTMHMM_OUT" ]; then
        [ -f "$DTMHMM_OUT/predicted_topologies.3line" ] && \
            cp "$DTMHMM_OUT/predicted_topologies.3line" "$OUTPUT_DIR/predicted_topologies.3line"
        [ -f "$DTMHMM_OUT/TMRs.gff3" ] && \
            cp "$DTMHMM_OUT/TMRs.gff3" "$OUTPUT_DIR/TMRs.gff3"
        if [ -f "$OUTPUT_DIR/predicted_topologies.3line" ]; then
            generate_prediction_tsv "$OUTPUT_DIR/predicted_topologies.3line" "$OUTPUT_DIR/prediction"
        fi
        rm -rf "$DTMHMM_OUT"
        echo "Local DeepTMHMM completed successfully" >&2
        exit 0
    fi
    echo "Warning: Local DeepTMHMM failed, falling back" >&2
    rm -rf "$DTMHMM_OUT"
fi

# --- Strategy 4: biolib CLI ---
if command -v biolib &>/dev/null; then
    echo "Using biolib DeepTMHMM" >&2
    BIOLIB_OUT="${OUTPUT_DIR}/biolib_output"
    mkdir -p "$BIOLIB_OUT"
    (cd "$BIOLIB_OUT" && biolib run DTU/DeepTMHMM --fasta "$FASTA") 2>/dev/null
    if [ -f "$BIOLIB_OUT/predicted_topologies.3line" ]; then
        cp "$BIOLIB_OUT/predicted_topologies.3line" "$OUTPUT_DIR/predicted_topologies.3line"
        [ -f "$BIOLIB_OUT/TMRs.gff3" ] && \
            cp "$BIOLIB_OUT/TMRs.gff3" "$OUTPUT_DIR/TMRs.gff3"
        generate_prediction_tsv "$OUTPUT_DIR/predicted_topologies.3line" "$OUTPUT_DIR/prediction"
        rm -rf "$BIOLIB_OUT"
        if [ -f "$OUTPUT_DIR/prediction" ]; then
            echo "biolib DeepTMHMM completed successfully" >&2
            exit 0
        fi
    fi
    echo "Warning: biolib DeepTMHMM output not in expected format, falling back" >&2
    rm -rf "$BIOLIB_OUT"
fi

# --- Strategy 5: Local Kyte-Doolittle predictor (fast fallback) ---
echo "WARNING: Using local Kyte-Doolittle TM predictor (less accurate than DeepTMHMM)" >&2
echo "Build the Apptainer container or install biolib for production accuracy" >&2
python3 "${SCRIPT_DIR}/predict_tm_helices.py" "$FASTA" "$OUTPUT_DIR/prediction" \
    --topology-output "$OUTPUT_DIR/predicted_topologies.3line"
exit $?
