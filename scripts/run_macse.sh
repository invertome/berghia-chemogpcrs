#!/bin/bash
# run_macse.sh — MACSE v2 codon-aware MSA.
#
# Bead -i61. Wrapper around MACSE v2 (Ranwez et al. 2018, MBE 35:2582).
# Codon-aware multiple sequence alignment — essential when input CDS may
# have frameshifts (e.g. miniprot-recovered CDS for non-Berghia species in
# stage 05). Naïve back-translation of frame-broken CDS via pal2nal alone
# inflates dN/dS estimates.
#
# Wraps the upstream MACSE_v2.X.X.jar — does NOT roll our own codon
# alignment. See https://www.agap-ge2pop.org/macse-v2/.
#
# Usage:
#   bash run_macse.sh --input=<cds.fna> --output=<aligned_codon.fna> \
#                     [--protein-output=<aligned_aa.fa>]
#
# Outputs:
#   <output>                       MACSE codon alignment (NT)
#   <protein-output>               aligned protein (optional)
#   <output>.predictor_used.txt    provenance

set -uo pipefail

INPUT=""
OUTPUT=""
PROTEIN_OUTPUT=""
MACSE_JAR="${MACSE_JAR:-${BASE_DIR:-.}/tools/macse_v2.07.jar}"
JAVA_BIN="${JAVA_BIN:-java}"
PROGRAM="${MACSE_PROGRAM:-alignSequences}"

usage() {
    cat >&2 <<EOF
run_macse.sh — MACSE v2 wrapper (bead -i61)

Required:
  --input=<file>            Input CDS FASTA (unaligned)
  --output=<file>            Output codon-alignment NT FASTA

Optional:
  --protein-output=<file>    Output aligned-protein FASTA (default: <output>.aa)
  --program=<name>           MACSE subprogram (default: alignSequences;
                             also: refineAlignment, exportAlignment, ...)

Env:
  MACSE_JAR  Path to macse_v2.XX.jar (default: \$BASE_DIR/tools/macse_v2.07.jar)
  JAVA_BIN   Java binary (default: java)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input=*)          INPUT="${1#*=}"; shift ;;
        --output=*)         OUTPUT="${1#*=}"; shift ;;
        --protein-output=*) PROTEIN_OUTPUT="${1#*=}"; shift ;;
        --program=*)        PROGRAM="${1#*=}"; shift ;;
        -h|--help)          usage; exit 0 ;;
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
if [ ! -f "$MACSE_JAR" ]; then
    echo "ERROR: MACSE jar not found: $MACSE_JAR" >&2
    echo "       See docs/installation_hpc.md for download URL." >&2
    exit 2
fi
if ! command -v "$JAVA_BIN" &>/dev/null; then
    echo "ERROR: Java binary '$JAVA_BIN' not found." >&2
    exit 2
fi

OUT_DIR=$(dirname "$OUTPUT")
mkdir -p "$OUT_DIR"
[ -z "$PROTEIN_OUTPUT" ] && PROTEIN_OUTPUT="${OUTPUT%.fna}.aa.fa"

# MACSE writes _NT.fasta + _AA.fasta with -out_NT / -out_AA flags.
# Use -prog alignSequences (default).
"$JAVA_BIN" -jar "$MACSE_JAR" \
    -prog "$PROGRAM" \
    -seq "$INPUT" \
    -out_NT "$OUTPUT" \
    -out_AA "$PROTEIN_OUTPUT" \
    > "${OUTPUT}.macse.log" 2>&1
RC=$?

if [ $RC -ne 0 ] || [ ! -s "$OUTPUT" ]; then
    echo "ERROR: MACSE failed (rc=$RC); see ${OUTPUT}.macse.log" >&2
    tail -10 "${OUTPUT}.macse.log" >&2 || true
    exit 3
fi

MACSE_VERSION=$(basename "$MACSE_JAR" | sed -n 's/macse_v\(.*\)\.jar/\1/p')
[ -z "$MACSE_VERSION" ] && MACSE_VERSION="unknown"
{
    echo "tool=macse_v2"
    echo "jar=${MACSE_JAR}"
    echo "version=${MACSE_VERSION}"
    echo "program=${PROGRAM}"
    echo "input=${INPUT}"
    echo "output_nt=${OUTPUT}"
    echo "output_aa=${PROTEIN_OUTPUT}"
    echo "run_at=$(date -Iseconds)"
} > "${OUTPUT}.predictor_used.txt"

echo "MACSE codon alignment: $OUTPUT (provenance: ${OUTPUT}.predictor_used.txt)" >&2
exit 0
