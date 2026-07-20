#!/bin/bash
# run_cloak.sh — CLOAK consensus mask wrapper.
#
# CLOAK (Cleaning on Alignment Konsensus, Chatur — Wheeler lab — bioRxiv
# Dec 2025) takes an ensemble of variant alignments of the same set of
# sequences and outputs a consensus alignment with disagreed columns
# masked. Performs best of all benchmarked filters for single-gene tree
# inference (Wheeler 2025).
#
# Implementation: cloak.py is a single Python script (stdlib only, no
# pip deps). It reads a directory of FASTA alignments and writes
# `<directory_basename>_result.fasta` to the current working directory.
# The first FASTA in alphabetical order serves as the column-structure
# reference (so name the canonical alignment `00_canonical.fa` to ensure
# the output is in canonical column order).
#
# Usage:
#   bash run_cloak.sh --ensemble-dir=<dir> --output=<consensus.fa>
#
# The ensemble directory must contain at least 2 FASTA files. The
# canonical alignment must be named so it sorts FIRST (e.g. canonical.fa
# vs variant_*.fa, OR 00_canonical.fa).

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

ENSEMBLE_DIR=""
OUTPUT=""
CLOAK_PY="${CLOAK:-}"

usage() {
    cat >&2 <<EOF
run_cloak.sh — CLOAK consensus-mask wrapper

Required:
  --ensemble-dir=<path>  Directory of variant alignments (>= 2 FASTAs)
  --output=<path>        Output consensus FASTA

Optional:
  --cloak-py=<path>      Path to cloak.py (default: \$CLOAK)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ensemble-dir=*) ENSEMBLE_DIR="${1#*=}"; shift ;;
        --output=*)       OUTPUT="${1#*=}"; shift ;;
        --cloak-py=*)     CLOAK_PY="${1#*=}"; shift ;;
        -h|--help)        usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -z "$ENSEMBLE_DIR" ]] && { echo "ERROR: --ensemble-dir required" >&2; exit 1; }
[[ -z "$OUTPUT"       ]] && { echo "ERROR: --output required"       >&2; exit 1; }
[[ ! -d "$ENSEMBLE_DIR" ]] && { echo "ERROR: ensemble dir not found: $ENSEMBLE_DIR" >&2; exit 1; }

if [[ -z "$CLOAK_PY" ]] || [[ ! -f "$CLOAK_PY" ]]; then
    echo "ERROR: cloak.py not found ($CLOAK_PY). Set \$CLOAK to the path of cloak.py." >&2
    exit 2
fi

# Count ensemble members
N_VARIANTS=$(find "$ENSEMBLE_DIR" -maxdepth 1 -name '*.fa' -type f | wc -l)
if [[ "$N_VARIANTS" -lt 2 ]]; then
    echo "ERROR: CLOAK requires >= 2 alignments in ensemble dir; found $N_VARIANTS" >&2
    exit 3
fi

# CLOAK writes its output and an `ordered_fasta_files_*` directory to CWD.
# Run in a sandbox tempdir so concurrent OG jobs don't collide.
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/cloak_XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT

# Stage the ensemble into the workdir under a name we control. CLOAK uses
# the directory basename to name its output, so we use a stable name.
STAGED_DIR="$WORKDIR/cloak_input"
mkdir -p "$STAGED_DIR"

# Force the canonical alignment to sort first (so CLOAK's output column
# structure follows the canonical alignment we'll use downstream).
# Named 00_canonical.fa so it sorts before 01_..., 02_... variants
# (ASCII: '0' < '1', and the two-digit prefix keeps the underscore from
# pushing 0_canonical after 01_ variants — underscore is ASCII 95, digits
# 49-57, so "0_" > "01" in ASCII ordering).
if [[ -f "$ENSEMBLE_DIR/canonical.fa" ]]; then
    cp "$ENSEMBLE_DIR/canonical.fa" "$STAGED_DIR/00_canonical.fa"
fi
i=1
for v in "$ENSEMBLE_DIR"/variant_*.fa; do
    [[ -f "$v" ]] || continue
    cp "$v" "$STAGED_DIR/$(printf '%02d' $i)_$(basename "$v")"
    i=$((i+1))
done

# Run cloak.py. It writes <basename>_result.fasta to CWD.
( cd "$WORKDIR" && python3 "$CLOAK_PY" "$STAGED_DIR" > cloak.log 2>&1 )
CLOAK_RC=$?

CLOAK_RESULT="$WORKDIR/cloak_input_result.fasta"
if [[ $CLOAK_RC -ne 0 ]] || [[ ! -s "$CLOAK_RESULT" ]]; then
    echo "ERROR: CLOAK failed (exit=$CLOAK_RC)" >&2
    tail -20 "$WORKDIR/cloak.log" >&2 || true
    exit 4
fi

mkdir -p "$(dirname "$OUTPUT")"
cp "$CLOAK_RESULT" "$OUTPUT"

# Provenance
{
    echo "tool=cloak"
    echo "n_alignments=$N_VARIANTS"
    echo "ensemble_dir=$ENSEMBLE_DIR"
    echo "output=$OUTPUT"
    n_seqs=$(grep -c '^>' "$OUTPUT" || true)
    n_seqs=${n_seqs:-0}
    echo "n_seqs=$n_seqs"
    echo "run_at=$(date -Iseconds)"
} > "${OUTPUT}.predictor_used.txt"

echo "CLOAK: $OUTPUT ($n_seqs seqs from $N_VARIANTS-variant ensemble)" >&2
exit 0
