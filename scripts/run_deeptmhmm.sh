#!/bin/bash
# run_deeptmhmm.sh - Unified TM prediction wrapper
# Bead -6nh: prevent silent fallback to Kyte-Doolittle, log which predictor
#            actually ran, and prefer TMbed over DeepTMHMM where available.
#
# Strategies (in order):
#   0. TMbed (Bernhofer & Rost 2022; current SOTA, +9pp recall vs DeepTMHMM)
#      Requires: pip install tmbed; first run downloads ProtT5 model (~3GB)
#   1. Apptainer SIF container (DEEPTMHMM_SIF or containers/deeptmhmm.sif)
#   2. Local venv installation (tools/deeptmhmm with .venv)
#   3. System deeptmhmm binary
#   4. Phobius (HMM-based, local, bioconda)
#   5. Kyte-Doolittle fallback — DISABLED by default (set
#      ALLOW_KD_FALLBACK=1 to opt in; produces less accurate results,
#      should NEVER be the silent default in a production pipeline run).
#
# Output: DeepTMHMM-compatible files:
#   prediction                      - TSV: seqid, pred_type, confidence, SP, n_tm_regions
#   predicted_topologies.3line      - Per-residue topology per sequence
#   predictor_used.txt              - Which strategy succeeded (provenance)
#
# Usage: bash run_deeptmhmm.sh -f INPUT.fasta -o OUTPUT_DIR

set -uo pipefail

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

FASTA=$(realpath "$FASTA")
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

ALLOW_KD_FALLBACK="${ALLOW_KD_FALLBACK:-0}"
TM_PREDICTOR_PRIMARY="${TM_PREDICTOR_PRIMARY:-tmbed}"  # tmbed | deeptmhmm

# --- Helper: Convert 3-line topology file to prediction TSV ---
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
        header = lines[i].strip().lstrip('>').split()[0]
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

record_predictor() {
    local name="$1"
    echo "$name" > "$OUTPUT_DIR/predictor_used.txt"
    echo "TM predictor used: $name (logged to $OUTPUT_DIR/predictor_used.txt)" >&2
}

# --- Strategy 0: TMbed (preferred when available) ---
if [ "$TM_PREDICTOR_PRIMARY" = "tmbed" ] && command -v tmbed &>/dev/null; then
    echo "Using TMbed (Bernhofer & Rost 2022)" >&2
    TMBED_OUT="${OUTPUT_DIR}/tmbed_out"
    mkdir -p "$TMBED_OUT"
    if tmbed predict -f "$FASTA" -p "$TMBED_OUT/predictions.3line" -f-format=fasta --out-format=3 2>"$TMBED_OUT/tmbed.log"; then
        cp "$TMBED_OUT/predictions.3line" "$OUTPUT_DIR/predicted_topologies.3line"
        generate_prediction_tsv "$OUTPUT_DIR/predicted_topologies.3line" "$OUTPUT_DIR/prediction"
        record_predictor "tmbed"
        rm -rf "$TMBED_OUT"
        exit 0
    fi
    echo "Warning: TMbed failed, trying DeepTMHMM strategies" >&2
    rm -rf "$TMBED_OUT"
fi

# --- Strategy 1: Apptainer SIF container ---
SIF="${DEEPTMHMM_SIF:-$BASE_DIR/containers/deeptmhmm.sif}"
if [ -f "$SIF" ]; then
    echo "Using Apptainer DeepTMHMM container: $SIF" >&2
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
        record_predictor "deeptmhmm-apptainer"
        rm -rf "$DTMHMM_OUT"
        exit 0
    fi
    echo "Warning: Apptainer DeepTMHMM failed (exit=$EXITCODE)" >&2
    rm -rf "$DTMHMM_OUT"
fi

# --- Strategy 2: Local DeepTMHMM venv installation ---
DTMHMM_INSTALL="${DEEPTMHMM_INSTALL:-$BASE_DIR/tools/deeptmhmm}"
if [ -f "$DTMHMM_INSTALL/predict.py" ] && [ -f "$DTMHMM_INSTALL/.venv/bin/python3" ]; then
    echo "Using local DeepTMHMM venv: $DTMHMM_INSTALL" >&2
    DTMHMM_OUT="${OUTPUT_DIR}/dtmhmm_run"
    rm -rf "$DTMHMM_OUT"
    (cd "$DTMHMM_INSTALL" && "$DTMHMM_INSTALL/.venv/bin/python3" predict.py \
        --fasta "$FASTA" --output-dir "$DTMHMM_OUT")
    EXITCODE=$?
    if [ $EXITCODE -eq 0 ] && [ -f "$DTMHMM_OUT/predicted_topologies.3line" ]; then
        cp "$DTMHMM_OUT/predicted_topologies.3line" "$OUTPUT_DIR/predicted_topologies.3line"
        [ -f "$DTMHMM_OUT/TMRs.gff3" ] && cp "$DTMHMM_OUT/TMRs.gff3" "$OUTPUT_DIR/TMRs.gff3"
        generate_prediction_tsv "$OUTPUT_DIR/predicted_topologies.3line" "$OUTPUT_DIR/prediction"
        record_predictor "deeptmhmm-venv"
        rm -rf "$DTMHMM_OUT"
        exit 0
    fi
    echo "Warning: local DeepTMHMM venv failed (exit=$EXITCODE)" >&2
    rm -rf "$DTMHMM_OUT"
fi

# --- Strategy 3: System DeepTMHMM binary ---
if command -v deeptmhmm &>/dev/null; then
    echo "Using system DeepTMHMM binary" >&2
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
            record_predictor "deeptmhmm-system"
            rm -rf "$DTMHMM_OUT"
            exit 0
        fi
    fi
    echo "Warning: system DeepTMHMM failed" >&2
    rm -rf "$DTMHMM_OUT"
fi

# --- Strategy 4: Phobius (HMM-based; bioconda; fast local fallback) ---
# DROPPED: biolib CLI (it makes a remote API call to biolib.com — slow,
# rate-limited, network-dependent). For batch HPC use we want only LOCAL
# predictors. Phobius (Käll et al. 2007) is HMM-based — older but reliable
# and trivial to install via bioconda.
if command -v phobius &>/dev/null; then
    echo "Using Phobius (HMM-based, local)" >&2
    PHOBIUS_OUT="${OUTPUT_DIR}/phobius.txt"
    if phobius -short "$FASTA" > "$PHOBIUS_OUT" 2>"${OUTPUT_DIR}/phobius.log"; then
        # Phobius -short emits one line per sequence: NAME TM SP PREDICTION
        # Convert to DeepTMHMM-compatible 3-line + prediction TSV.
        python3 - "$PHOBIUS_OUT" "$OUTPUT_DIR/prediction" "$OUTPUT_DIR/predicted_topologies.3line" <<'PYEOF'
import sys
short, pred_tsv, three_line = sys.argv[1], sys.argv[2], sys.argv[3]
with open(short) as f, open(pred_tsv, 'w') as out_pred, open(three_line, 'w') as out_3l:
    for line in f:
        line = line.strip()
        if not line or line.startswith('SEQENCE'):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        name, tm_count, sp = parts[0], parts[1], parts[2]
        topology_str = ' '.join(parts[3:])
        # DeepTMHMM-compatible 5-col TSV
        try:
            n_tm = int(tm_count)
        except ValueError:
            n_tm = 0
        pred_type = 'TM' if n_tm > 0 else 'GLOB'
        out_pred.write(f"{name}\t{pred_type}\t1.0\t\t{n_tm}\n")
        # 3-line approximation: header / sequence (unknown) / topology-string
        out_3l.write(f">{name}\nN\n{topology_str}\n")
PYEOF
        if [ -s "$OUTPUT_DIR/prediction" ]; then
            record_predictor "phobius"
            exit 0
        fi
    fi
    echo "Warning: Phobius failed; trying next strategy" >&2
fi

# --- Strategy 5: Kyte-Doolittle (DISABLED unless explicitly opted in) ---
if [ "$ALLOW_KD_FALLBACK" = "1" ]; then
    echo "WARNING: ALLOW_KD_FALLBACK=1 — using legacy Kyte-Doolittle predictor." >&2
    echo "         KD is significantly less accurate than DeepTMHMM/TMbed." >&2
    echo "         This should NEVER be the production setting (bead -6nh)." >&2
    if [ -f "${SCRIPT_DIR}/legacy/predict_tm_helices.py" ]; then
        KD_SCRIPT="${SCRIPT_DIR}/legacy/predict_tm_helices.py"
    else
        KD_SCRIPT="${SCRIPT_DIR}/predict_tm_helices.py"
    fi
    python3 "$KD_SCRIPT" "$FASTA" "$OUTPUT_DIR/prediction" \
        --topology-output "$OUTPUT_DIR/predicted_topologies.3line"
    EXITCODE=$?
    if [ $EXITCODE -eq 0 ]; then
        record_predictor "kyte-doolittle (LEGACY, opt-in via ALLOW_KD_FALLBACK=1)"
    fi
    exit $EXITCODE
fi

# --- No predictor available ---
echo "" >&2
echo "ERROR: No TM predictor available. Tried (in order):" >&2
[ "$TM_PREDICTOR_PRIMARY" = "tmbed" ] && echo "  - TMbed (preferred, local)" >&2
echo "  - Apptainer DeepTMHMM container ($SIF)" >&2
echo "  - Local DeepTMHMM venv ($DTMHMM_INSTALL)" >&2
echo "  - System 'deeptmhmm' binary" >&2
echo "  - Phobius (HMM-based, local)" >&2
echo "" >&2
echo "Install one of the following (LOCAL — no remote APIs):" >&2
echo "  # TMbed (preferred, current SOTA, +9pp recall):" >&2
echo "  git clone https://github.com/BernhoferM/TMbed && pip install ./TMbed" >&2
echo "  # Phobius (fast HMM-based fallback, bioconda):" >&2
echo "  conda install -c bioconda phobius" >&2
echo "  # DeepTMHMM Apptainer (build SIF — see docs/installation_hpc.md)" >&2
echo "" >&2
echo "If you absolutely must run without a real predictor for testing," >&2
echo "set ALLOW_KD_FALLBACK=1 (Kyte-Doolittle is significantly less accurate" >&2
echo "and should NEVER be used in a production pipeline run; bead -6nh)." >&2
exit 2
