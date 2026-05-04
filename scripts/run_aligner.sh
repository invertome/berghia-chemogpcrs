#!/bin/bash
# run_aligner.sh — Regime-based MSA wrapper.
#
# Picks the scientifically-best aligner per input size:
#
#   N <  ALIGNER_LINSI_THRESHOLD  (default 200)   → MAFFT L-INS-i
#                                                   (--localpair --maxiterate 1000)
#                                                   Highest-accuracy small alignment;
#                                                   appropriate for per-orthogroup
#                                                   GPCR families where every column
#                                                   matters for downstream dN/dS.
#
#   N <  ALIGNER_FAMSA_THRESHOLD  (default 1000)  → MAFFT --auto
#                                                   MAFFT picks L-INS-i / G-INS-i /
#                                                   FFT-NS-i / FFT-NS-2 based on
#                                                   sequence properties.
#
#   N >= ALIGNER_FAMSA_THRESHOLD                  → FAMSA 2
#                                                   (Deorowicz 2024 NAR 52:e30)
#                                                   On HomFam ~3pp better SP-score
#                                                   than MAFFT --auto at this scale,
#                                                   ~5-10× faster than MAFFT --retree 2
#                                                   on large protein families.
#
# Override via ALIGNER_BACKEND={auto,mafft_linsi,mafft_auto,mafft_retree2,famsa}.
#
# Falls back to MAFFT --retree 2 if the preferred backend is unavailable
# (e.g. FAMSA not installed). Logs which backend ran, what N was, and
# elapsed wall time to <output_dir>/aligner_used.txt.
#
# Usage:
#   bash run_aligner.sh --input=<fasta> --output=<aligned.fa> [--threads=N]
#   bash run_aligner.sh --input=<fasta> --output=<aligned.fa> --force-backend=famsa
#
# Exit codes: 0 success, 2 missing tool with no fallback, 3 alignment failed.

set -uo pipefail

INPUT=""
OUTPUT=""
THREADS="${ALIGNER_THREADS:-${CPUS:-4}}"
FORCE_BACKEND="${ALIGNER_BACKEND:-auto}"

# Backend binaries (overridable for HPC)
MAFFT_BIN="${MAFFT:-mafft}"
FAMSA_BIN="${FAMSA:-famsa}"

# Thresholds
LINSI_THRESHOLD="${ALIGNER_LINSI_THRESHOLD:-200}"
FAMSA_THRESHOLD="${ALIGNER_FAMSA_THRESHOLD:-1000}"

usage() {
    cat >&2 <<EOF
run_aligner.sh — regime-based MSA wrapper (bead -align)

Required:
  --input=<fasta>           Input FASTA
  --output=<file>           Output aligned FASTA

Optional:
  --threads=<N>             Threads (default: \$ALIGNER_THREADS or \$CPUS or 4)
  --force-backend=<name>    Override auto-selection. One of:
                              auto (default), mafft_linsi, mafft_auto,
                              mafft_retree2, famsa
  --linsi-threshold=<N>     N below which MAFFT L-INS-i is used (default 200)
  --famsa-threshold=<N>     N at/above which FAMSA 2 is used (default 1000)

Env overrides:
  MAFFT, FAMSA              binary names/paths
  ALIGNER_BACKEND           same as --force-backend
  ALIGNER_LINSI_THRESHOLD   same as --linsi-threshold
  ALIGNER_FAMSA_THRESHOLD   same as --famsa-threshold

Selection (when --force-backend=auto):
  N < LINSI_THRESHOLD          MAFFT L-INS-i (highest accuracy small)
  LINSI <= N < FAMSA           MAFFT --auto (MAFFT picks mode)
  N >= FAMSA_THRESHOLD         FAMSA 2 (best for large)

Fallback chain when preferred backend missing:
  L-INS-i / --auto fail / MAFFT missing  ->  MAFFT --retree 2  ->  exit 2
  FAMSA missing                          ->  MAFFT --auto       ->  --retree 2 -> exit 2
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input=*)             INPUT="${1#*=}"; shift ;;
        --output=*)            OUTPUT="${1#*=}"; shift ;;
        --threads=*)           THREADS="${1#*=}"; shift ;;
        --force-backend=*)     FORCE_BACKEND="${1#*=}"; shift ;;
        --linsi-threshold=*)   LINSI_THRESHOLD="${1#*=}"; shift ;;
        --famsa-threshold=*)   FAMSA_THRESHOLD="${1#*=}"; shift ;;
        -h|--help)             usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "ERROR: --input and --output are required" >&2
    usage
    exit 1
fi
if [ ! -f "$INPUT" ]; then
    echo "ERROR: input not found: $INPUT" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"

# Count sequences
N=$(grep -c '^>' "$INPUT" 2>/dev/null || echo 0)

# Decide backend
choose_backend() {
    if [ "$FORCE_BACKEND" != "auto" ]; then
        echo "$FORCE_BACKEND"
        return
    fi
    if [ "$N" -lt "$LINSI_THRESHOLD" ]; then
        echo "mafft_linsi"
    elif [ "$N" -lt "$FAMSA_THRESHOLD" ]; then
        echo "mafft_auto"
    else
        echo "famsa"
    fi
}

BACKEND=$(choose_backend)

# Tool availability checks + fallback chain
have_mafft() { command -v "$MAFFT_BIN" &>/dev/null; }
have_famsa() { command -v "$FAMSA_BIN" &>/dev/null; }

if [ "$BACKEND" = "famsa" ] && ! have_famsa; then
    echo "WARN: FAMSA ($FAMSA_BIN) not found; falling back to MAFFT --auto" >&2
    BACKEND="mafft_auto"
fi
if [[ "$BACKEND" == mafft_* ]] && ! have_mafft; then
    echo "ERROR: MAFFT ($MAFFT_BIN) not found and chosen backend was $BACKEND" >&2
    exit 2
fi

# Run the chosen aligner
START=$(date +%s)
case "$BACKEND" in
    mafft_linsi)
        # L-INS-i: localpair + iterative refinement, gold-standard accuracy
        # for small alignments. Time per N^2 — practical up to ~500.
        "$MAFFT_BIN" --localpair --maxiterate 1000 --thread "$THREADS" \
            "$INPUT" > "$OUTPUT" 2>"${OUTPUT}.aligner.log"
        RC=$?
        ;;
    mafft_auto)
        # --auto: MAFFT picks L-INS-i / G-INS-i / FFT-NS-i / FFT-NS-2
        # based on N + sequence-similarity heuristics.
        "$MAFFT_BIN" --auto --thread "$THREADS" \
            "$INPUT" > "$OUTPUT" 2>"${OUTPUT}.aligner.log"
        RC=$?
        ;;
    mafft_retree2)
        "$MAFFT_BIN" --retree 2 --thread "$THREADS" \
            "$INPUT" > "$OUTPUT" 2>"${OUTPUT}.aligner.log"
        RC=$?
        ;;
    famsa)
        # FAMSA 2: -t threads. Reads FASTA, writes FASTA. Includes built-in
        # tree construction; we don't need separate guide tree.
        "$FAMSA_BIN" -t "$THREADS" "$INPUT" "$OUTPUT" 2>"${OUTPUT}.aligner.log"
        RC=$?
        ;;
    *)
        echo "ERROR: unknown backend '$BACKEND'" >&2
        exit 1
        ;;
esac
END=$(date +%s)
ELAPSED=$((END - START))

if [ $RC -ne 0 ] || [ ! -s "$OUTPUT" ]; then
    echo "ERROR: $BACKEND failed (rc=$RC) on $INPUT" >&2
    tail -10 "${OUTPUT}.aligner.log" >&2 || true
    # Fallback to MAFFT --retree 2 if we haven't already tried it
    if [ "$BACKEND" != "mafft_retree2" ] && have_mafft; then
        echo "WARN: falling back to MAFFT --retree 2" >&2
        "$MAFFT_BIN" --retree 2 --thread "$THREADS" \
            "$INPUT" > "$OUTPUT" 2>>"${OUTPUT}.aligner.log"
        if [ $? -ne 0 ] || [ ! -s "$OUTPUT" ]; then
            exit 3
        fi
        BACKEND="mafft_retree2 (fallback)"
    else
        exit 3
    fi
fi

# Provenance — bead -ryr style
PROV="${OUTPUT}.aligner_used.txt"
{
    echo "tool=run_aligner.sh"
    echo "backend=${BACKEND}"
    echo "n_sequences=${N}"
    echo "linsi_threshold=${LINSI_THRESHOLD}"
    echo "famsa_threshold=${FAMSA_THRESHOLD}"
    echo "threads=${THREADS}"
    echo "elapsed_seconds=${ELAPSED}"
    echo "input=${INPUT}"
    echo "output=${OUTPUT}"
    echo "run_at=$(date -Iseconds)"
} > "$PROV"

echo "Aligner: ${BACKEND} on N=${N} seqs in ${ELAPSED}s -> $OUTPUT" >&2
exit 0
