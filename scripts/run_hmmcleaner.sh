#!/bin/bash
# run_hmmcleaner.sh — HmmCleaner per-sequence segment cleaning.
#
# Bead -i61. Wrapper around HmmCleaner.pl (Di Franco et al. 2019,
# BMC Evol Biol 19:21). Detects misaligned/contaminated segments WITHIN
# individual sequences (column trimming alone, e.g. ClipKit, doesn't
# catch these). Two-stage filtering = HmmCleaner (segment) + ClipKit
# (column) is the modern best practice for paralog-rich gene families.
#
# Wraps the upstream HmmCleaner.pl Perl script — does NOT roll our own
# segment detection. See https://metacpan.org/dist/Bio-MUST-Apps-HmmCleaner.
#
# Usage:
#   bash run_hmmcleaner.sh --input=<aln.fa> --output=<cleaned_aln.fa> [--threshold=N]
#
# Outputs:
#   <output>                       cleaned alignment FASTA
#   <output>.score                 per-sequence cleaning report
#   <output>.predictor_used.txt    provenance (tool + version)

set -uo pipefail

INPUT=""
OUTPUT=""
THRESHOLD="${HMMCLEANER_THRESHOLD:-}"
HMMCLEANER_BIN="${HMMCLEANER:-HmmCleaner.pl}"

usage() {
    cat >&2 <<EOF
run_hmmcleaner.sh — HmmCleaner wrapper (bead -i61)

Required:
  --input=<file>     Input MAFFT/MSA alignment (FASTA)
  --output=<file>    Output cleaned alignment path

Optional:
  --threshold=<N>    HmmCleaner threshold (default: tool default)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input=*)     INPUT="${1#*=}"; shift ;;
        --output=*)    OUTPUT="${1#*=}"; shift ;;
        --threshold=*) THRESHOLD="${1#*=}"; shift ;;
        -h|--help)     usage; exit 0 ;;
        *)             echo "Unknown argument: $1" >&2; usage; exit 1 ;;
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
if ! command -v "$HMMCLEANER_BIN" &>/dev/null; then
    echo "ERROR: HmmCleaner.pl not found ($HMMCLEANER_BIN). " \
         "See docs/installation_hpc.md for install instructions." >&2
    exit 2
fi

OUT_DIR=$(dirname "$OUTPUT")
mkdir -p "$OUT_DIR"

# HmmCleaner writes <input_basename>_hmm.fasta (and .score) alongside the
# input by default. Run in a temp workdir so we don't pollute the source.
WORKDIR=$(mktemp -d "${TMPDIR:-/tmp}/hmmcleaner_XXXXXX")
trap "rm -rf '$WORKDIR'" EXIT

cp -L "$INPUT" "$WORKDIR/input.fasta"
pushd "$WORKDIR" >/dev/null

THRESHOLD_ARG=()
if [ -n "$THRESHOLD" ]; then
    THRESHOLD_ARG=(--threshold "$THRESHOLD")
fi

if ! "$HMMCLEANER_BIN" "${THRESHOLD_ARG[@]}" input.fasta > hmmcleaner.log 2>&1; then
    popd >/dev/null
    echo "ERROR: HmmCleaner.pl failed (see $WORKDIR/hmmcleaner.log)" >&2
    tail -10 "$WORKDIR/hmmcleaner.log" >&2 || true
    exit 3
fi

# Locate cleaned output
CLEANED=$(ls input_hmm.fasta 2>/dev/null | head -1)
SCORE=$(ls input_hmm.score 2>/dev/null | head -1)
[ -z "$CLEANED" ] && CLEANED=$(ls input_hmm.fa 2>/dev/null | head -1)

popd >/dev/null

if [ -z "$CLEANED" ] || [ ! -s "$WORKDIR/$CLEANED" ]; then
    echo "ERROR: HmmCleaner produced no cleaned alignment in $WORKDIR" >&2
    ls -la "$WORKDIR" >&2
    exit 4
fi

cp "$WORKDIR/$CLEANED" "$OUTPUT"
[ -n "$SCORE" ] && [ -f "$WORKDIR/$SCORE" ] && cp "$WORKDIR/$SCORE" "${OUTPUT}.score"

# Provenance
HMMCLEANER_VERSION=$("$HMMCLEANER_BIN" --version 2>/dev/null | head -1 || echo "unknown")
{
    echo "tool=hmmcleaner"
    echo "binary=${HMMCLEANER_BIN}"
    echo "version=${HMMCLEANER_VERSION}"
    echo "threshold=${THRESHOLD:-default}"
    echo "input=${INPUT}"
    echo "output=${OUTPUT}"
    echo "run_at=$(date -Iseconds)"
} > "${OUTPUT}.predictor_used.txt"

echo "HmmCleaner: $OUTPUT (provenance: ${OUTPUT}.predictor_used.txt)" >&2
exit 0
