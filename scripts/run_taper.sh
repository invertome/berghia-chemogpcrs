#!/bin/bash
# run_taper.sh — TAPER per-sequence outlier-residue masking wrapper.
#
# TAPER (Zhang, Zhao, Braun & Mirarab 2021, Methods Ecol Evol 13:91)
# uses a 2D outlier-detection algorithm to find residues that are
# species-specific outliers given the overall divergence pattern. Masks
# them with `X`. Designed to catch small per-sequence error stretches
# (sequencing errors, late-stage assembly artifacts) that PREQUAL +
# CLOAK miss.
#
# Implementation: TAPER is a Julia script (`correction_multi.jl`)
# distributed via github chaoszhang/TAPER. We invoke it via the env's
# julia binary.
#
# Usage:
#   bash run_taper.sh --input=<aln.fa> --output=<masked_aln.fa>

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

INPUT=""
OUTPUT=""
TAPER_JL="${TAPER:-}"
JULIA_BIN="${JULIA:-julia}"

usage() {
    cat >&2 <<EOF
run_taper.sh — TAPER residue-outlier masking wrapper

Required:
  --input=<path>     Aligned protein FASTA
  --output=<path>    Output FASTA with outlier residues masked as X

Optional:
  --taper-jl=<path>  Path to correction_multi.jl (default: \$TAPER)
  --julia-bin=<path> Julia binary (default: \$JULIA or 'julia')
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input=*)     INPUT="${1#*=}"; shift ;;
        --output=*)    OUTPUT="${1#*=}"; shift ;;
        --taper-jl=*)  TAPER_JL="${1#*=}"; shift ;;
        --julia-bin=*) JULIA_BIN="${1#*=}"; shift ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -z "$INPUT"  ]] && { echo "ERROR: --input required"  >&2; exit 1; }
[[ -z "$OUTPUT" ]] && { echo "ERROR: --output required" >&2; exit 1; }
[[ ! -f "$INPUT" ]] && { echo "ERROR: input not found: $INPUT" >&2; exit 1; }

if [[ -z "$TAPER_JL" ]] || [[ ! -f "$TAPER_JL" ]]; then
    echo "ERROR: TAPER correction_multi.jl not found ($TAPER_JL). Set \$TAPER to its path." >&2
    exit 2
fi

if ! command -v "$JULIA_BIN" >/dev/null 2>&1; then
    echo "ERROR: julia not found ($JULIA_BIN). Install: mamba install -c conda-forge julia" >&2
    exit 2
fi

mkdir -p "$(dirname "$OUTPUT")"

# TAPER prints to stdout. Capture stderr separately.
TAPER_LOG="${OUTPUT}.taper.log"
"$JULIA_BIN" "$TAPER_JL" "$INPUT" > "$OUTPUT" 2>"$TAPER_LOG"
TAPER_RC=$?

if [[ $TAPER_RC -ne 0 ]] || [[ ! -s "$OUTPUT" ]]; then
    echo "ERROR: TAPER failed (exit=$TAPER_RC)" >&2
    tail -20 "$TAPER_LOG" >&2 || true
    exit 3
fi

# Provenance
{
    echo "tool=taper"
    echo "julia_version=$("$JULIA_BIN" --version 2>/dev/null | head -1)"
    echo "input=$INPUT"
    echo "output=$OUTPUT"
    n_in=$(grep -c '^>' "$INPUT" || true)
    n_in=${n_in:-0}
    n_out=$(grep -c '^>' "$OUTPUT" || true)
    n_out=${n_out:-0}
    echo "input_seqs=$n_in"
    echo "output_seqs=$n_out"
    echo "run_at=$(date -Iseconds)"
} > "${OUTPUT}.predictor_used.txt"

echo "TAPER: $OUTPUT (${n_in} -> ${n_out} sequences; provenance: ${OUTPUT}.predictor_used.txt)" >&2
exit 0
