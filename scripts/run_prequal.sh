#!/bin/bash
# run_prequal.sh — PREQUAL pre-alignment quality filter wrapper.
#
# PREQUAL (Whelan, Irisarri & Burki 2018, Bioinformatics 34:3929) reads a
# FASTA of unaligned protein sequences and replaces residues whose
# pairwise homology is uncertain with `X`. The masked sequences then go
# into the aligner; downstream tools (ClipKit, MACSE/pal2nal) treat `X`
# correctly so masked positions become NNN codons in the codon
# alignment and are ignored by HyPhy as missing data.
#
# Why this is in the pipeline (replaces HmmCleaner, May 2026):
#   HmmCleaner's Perl XS dependency stack does not build under the conda
#   env's Perl ABI on Unity (no sudo, no system headers). PREQUAL is the
#   peer-reviewed, conservative, conda-installable alternative. It runs
#   before MAFFT/FAMSA so masking is at the residue level on the input.
#
# Usage:
#   bash run_prequal.sh --input=<seqs.fa> --output=<masked.fa>
#
# Output: masked FASTA (same headers as input, X-replaced residues).

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

INPUT=""
OUTPUT=""

usage() {
    cat >&2 <<EOF
run_prequal.sh — PREQUAL wrapper

Required:
  --input=<path>    Unaligned protein FASTA
  --output=<path>   Output FASTA with low-confidence residues masked as X

Optional:
  --prequal-bin=<path>  PREQUAL binary (default: \$PREQUAL or 'prequal')
EOF
}

PREQUAL_BIN="${PREQUAL:-prequal}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input=*)        INPUT="${1#*=}"; shift ;;
        --output=*)       OUTPUT="${1#*=}"; shift ;;
        --prequal-bin=*)  PREQUAL_BIN="${1#*=}"; shift ;;
        -h|--help)        usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -z "$INPUT"  ]] && { echo "ERROR: --input required"  >&2; exit 1; }
[[ -z "$OUTPUT" ]] && { echo "ERROR: --output required" >&2; exit 1; }
[[ ! -f "$INPUT" ]] && { echo "ERROR: input not found: $INPUT" >&2; exit 1; }

if ! command -v "$PREQUAL_BIN" >/dev/null 2>&1; then
    echo "ERROR: prequal binary not found: $PREQUAL_BIN" >&2
    echo "Install: mamba install -c bioconda prequal" >&2
    exit 2
fi

# PREQUAL writes <input>.filtered alongside the input. Run in a tempdir so
# we don't pollute the input directory and so concurrent jobs don't collide
# on the same auto-named output.
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/prequal_XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT

cp -L "$INPUT" "$WORKDIR/seqs.fa"

# PREQUAL default -filterthresh = 0.994 (very strict — only retains
# residues with >=99.4% homology posterior). On divergent chemoreceptor
# data (Berghia + ~2000 mollusc refs), this dropped 71.5% of sequences
# in the all_berghia_refs tree (2829 -> 805) and was the root cause of
# the apparent "PREQUAL failed" WARN cascade in stage 04 array tasks
# (job 57844031, 2026-05-18). The original paper (Whelan 2018) recommends
# 0.95 for typical sets and looser for divergent families. Use 0.90 by
# default; override via PREQUAL_FILTERTHRESH env var.
PREQUAL_FILTERTHRESH="${PREQUAL_FILTERTHRESH:-0.90}"

( cd "$WORKDIR" && "$PREQUAL_BIN" -filterthresh "$PREQUAL_FILTERTHRESH" seqs.fa > prequal.log 2>&1 )
PREQUAL_RC=$?

if [[ $PREQUAL_RC -ne 0 ]] || [[ ! -s "$WORKDIR/seqs.fa.filtered" ]]; then
    echo "ERROR: PREQUAL failed (exit=$PREQUAL_RC)" >&2
    tail -20 "$WORKDIR/prequal.log" >&2 || true
    exit 3
fi

mkdir -p "$(dirname "$OUTPUT")"
cp "$WORKDIR/seqs.fa.filtered" "$OUTPUT"

# Provenance
{
    echo "tool=prequal"
    echo "version=$("$PREQUAL_BIN" -h 2>&1 | grep -oE 'PREQUAL[^[:space:]]*' | head -1 || echo unknown)"
    echo "input=$INPUT"
    echo "output=$OUTPUT"
    n_in=$(grep -c '^>' "$INPUT" 2>/dev/null || echo 0)
    n_out=$(grep -c '^>' "$OUTPUT" 2>/dev/null || echo 0)
    echo "input_seqs=$n_in"
    echo "output_seqs=$n_out"
    echo "run_at=$(date -Iseconds)"
} > "${OUTPUT}.predictor_used.txt"

echo "PREQUAL: $OUTPUT (${n_in} -> ${n_out} sequences; provenance: ${OUTPUT}.predictor_used.txt)" >&2
exit 0
