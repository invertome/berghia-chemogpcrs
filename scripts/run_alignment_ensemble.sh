#!/bin/bash
# run_alignment_ensemble.sh — generate K alignments of the same input for CLOAK.
#
# CLOAK (Chatur, Wheeler lab, bioRxiv Dec 2025) detects alignment
# uncertainty by running multiple aligners/parameters on the same input
# and masking columns where the alignments disagree. This script
# produces the ensemble.
#
# Default ensemble (6 alignments — chosen for chemoreceptor GPCRs where
# parameter-sensitivity in ECL/ICL loops dominates uncertainty):
#   1. MAFFT canonical (regime-based, kept for downstream tree)
#   2. MAFFT --localpair --maxiterate 1000  (LINSI)
#   3. MAFFT --globalpair --maxiterate 1000 (GINSI)
#   4. MAFFT --genafpair --maxiterate 1000  (EINSI; best for multi-domain
#                                            sequences with long inter-
#                                            domain spacers — exactly
#                                            7TM proteins)
#   5. MAFFT --retree 1 --maxiterate 0      (single-tree progressive)
#   6. FAMSA                                (independent algorithm family)
#
# Usage:
#   bash run_alignment_ensemble.sh \
#     --input=<seqs.fa> \
#     --output-dir=<dir>      # gets <dir>/canonical.fa + <dir>/variant_{1..N}.fa
#     [--threads=<int>]
#     [--canonical-aligner=<auto|linsi|ginsi|einsi|famsa>]
#
# Output: a directory containing one FASTA per ensemble member, with
# `canonical.fa` named first so CLOAK's reference-alignment column
# structure matches the alignment we keep for downstream tree building.

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

INPUT=""
OUTPUT_DIR=""
THREADS="${CPUS:-4}"
CANONICAL_ALIGNER="auto"

usage() {
    cat >&2 <<EOF
run_alignment_ensemble.sh — generate ensemble of MAFFT variants + FAMSA

Required:
  --input=<path>          Unaligned protein FASTA
  --output-dir=<path>     Directory to write canonical.fa + variant_*.fa

Optional:
  --threads=<int>         Threads per aligner (default \$CPUS or 4)
  --canonical-aligner=<auto|linsi|ginsi|einsi|famsa>
                          Picks which member is the canonical (default auto:
                          regime-based — uses scripts/run_aligner.sh).
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input=*)              INPUT="${1#*=}"; shift ;;
        --output-dir=*)         OUTPUT_DIR="${1#*=}"; shift ;;
        --threads=*)            THREADS="${1#*=}"; shift ;;
        --canonical-aligner=*)  CANONICAL_ALIGNER="${1#*=}"; shift ;;
        -h|--help)              usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -z "$INPUT"      ]] && { echo "ERROR: --input required"      >&2; exit 1; }
[[ -z "$OUTPUT_DIR" ]] && { echo "ERROR: --output-dir required" >&2; exit 1; }
[[ ! -f "$INPUT" ]]    && { echo "ERROR: input not found: $INPUT" >&2; exit 1; }

MAFFT_BIN="${MAFFT:-mafft}"
FAMSA_BIN="${FAMSA:-famsa}"

if ! command -v "$MAFFT_BIN" >/dev/null 2>&1; then
    echo "ERROR: mafft not found: $MAFFT_BIN" >&2
    exit 2
fi

mkdir -p "$OUTPUT_DIR"

# --- Canonical alignment (regime-based via run_aligner.sh, OR a chosen mode) ---
case "$CANONICAL_ALIGNER" in
    auto)
        bash "${SCRIPT_DIR}/run_aligner.sh" \
            --input="$INPUT" \
            --output="${OUTPUT_DIR}/canonical.fa" \
            --threads="$THREADS"
        ;;
    linsi) "$MAFFT_BIN" --localpair  --maxiterate 1000 --thread "$THREADS" "$INPUT" > "${OUTPUT_DIR}/canonical.fa" 2>"${OUTPUT_DIR}/canonical.log" ;;
    ginsi) "$MAFFT_BIN" --globalpair --maxiterate 1000 --thread "$THREADS" "$INPUT" > "${OUTPUT_DIR}/canonical.fa" 2>"${OUTPUT_DIR}/canonical.log" ;;
    einsi) "$MAFFT_BIN" --genafpair  --maxiterate 1000 --thread "$THREADS" "$INPUT" > "${OUTPUT_DIR}/canonical.fa" 2>"${OUTPUT_DIR}/canonical.log" ;;
    famsa)
        if ! command -v "$FAMSA_BIN" >/dev/null 2>&1; then
            echo "ERROR: famsa not found for --canonical-aligner=famsa" >&2; exit 2
        fi
        "$FAMSA_BIN" -t "$THREADS" "$INPUT" "${OUTPUT_DIR}/canonical.fa" 2>"${OUTPUT_DIR}/canonical.log"
        ;;
    *) echo "Unknown canonical aligner: $CANONICAL_ALIGNER" >&2; exit 1 ;;
esac

if [[ ! -s "${OUTPUT_DIR}/canonical.fa" ]]; then
    echo "ERROR: canonical alignment empty: ${OUTPUT_DIR}/canonical.fa" >&2
    exit 3
fi

# --- Variant alignments ---
# Each variant runs in series; failures emit a warning and skip the variant
# rather than aborting (CLOAK degrades gracefully with fewer ensemble members,
# down to K=2).
declare -a VARIANTS=(
    "linsi:--localpair --maxiterate 1000"
    "ginsi:--globalpair --maxiterate 1000"
    "einsi:--genafpair --maxiterate 1000"
    "fftnsi:--retree 1 --maxiterate 0"
)

run_mafft_variant() {
    local name="$1"; shift
    local out="${OUTPUT_DIR}/variant_${name}.fa"
    local log="${OUTPUT_DIR}/variant_${name}.log"
    if "$MAFFT_BIN" "$@" --thread "$THREADS" "$INPUT" > "$out" 2>"$log" \
       && [[ -s "$out" ]]; then
        echo "  variant_${name}: OK ($(grep -c '^>' "$out") seqs)" >&2
    else
        echo "  variant_${name}: FAILED (continuing without it)" >&2
        rm -f "$out"
    fi
}

echo "[ensemble] Generating MAFFT variants..." >&2
for entry in "${VARIANTS[@]}"; do
    name="${entry%%:*}"
    args="${entry#*:}"
    # shellcheck disable=SC2086
    run_mafft_variant "$name" $args
done

# FAMSA variant (different algorithm family)
if command -v "$FAMSA_BIN" >/dev/null 2>&1; then
    if "$FAMSA_BIN" -t "$THREADS" "$INPUT" "${OUTPUT_DIR}/variant_famsa.fa" \
        2>"${OUTPUT_DIR}/variant_famsa.log" \
       && [[ -s "${OUTPUT_DIR}/variant_famsa.fa" ]]; then
        echo "  variant_famsa: OK ($(grep -c '^>' "${OUTPUT_DIR}/variant_famsa.fa") seqs)" >&2
    else
        echo "  variant_famsa: FAILED (continuing without it)" >&2
        rm -f "${OUTPUT_DIR}/variant_famsa.fa"
    fi
else
    echo "  variant_famsa: SKIPPED (famsa not installed)" >&2
fi

# --- Provenance ---
N_VARIANTS=$(find "$OUTPUT_DIR" -maxdepth 1 -name 'variant_*.fa' -type f | wc -l)
N_TOTAL=$((N_VARIANTS + 1))
{
    echo "tool=run_alignment_ensemble"
    echo "canonical_aligner=$CANONICAL_ALIGNER"
    echo "n_total_alignments=$N_TOTAL  (canonical + $N_VARIANTS variants)"
    echo "input=$INPUT"
    echo "output_dir=$OUTPUT_DIR"
    echo "threads=$THREADS"
    echo "run_at=$(date -Iseconds)"
} > "${OUTPUT_DIR}/predictor_used.txt"

echo "[ensemble] DONE: $N_TOTAL alignments in $OUTPUT_DIR" >&2
exit 0
