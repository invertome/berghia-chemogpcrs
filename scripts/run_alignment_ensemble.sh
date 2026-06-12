#!/bin/bash
# run_alignment_ensemble.sh — generate K alignments of the same input for CLOAK.
#
# CLOAK (Chatur, Wheeler lab, bioRxiv Dec 2025) detects alignment
# uncertainty by running multiple aligners/parameters on the same input
# and masking columns where the alignments disagree. This script
# produces the ensemble.
#
# Default ensemble (5 alignments — all MAFFT variants):
#   1. MAFFT canonical (regime-based, kept for downstream tree)
#   2. MAFFT --localpair --maxiterate 1000  (LINSI)
#   3. MAFFT --globalpair --maxiterate 1000 (GINSI)
#   4. MAFFT --genafpair --maxiterate 1000  (EINSI; best for multi-domain
#                                            sequences with long inter-
#                                            domain spacers — exactly
#                                            7TM proteins)
#   5. MAFFT --retree 1 --maxiterate 0      (single-tree progressive)
#
# WHY MAFFT-ONLY (lesson from large smoke test 56832359, 2026-05-06):
# An earlier version of this script included FAMSA as a 6th ensemble
# member ("for cross-algorithm diversity"). On 150 Conus consors seqs,
# FAMSA produced a 3772-column alignment vs MAFFT variants' 1627-1999
# columns (2.3x longer). CLOAK's column-pairing consensus across that
# heterogeneous ensemble inflated the output to 6350 columns instead
# of masking down toward ~1600 — the OPPOSITE of what we wanted.
#
# Wheeler 2025 (CLOAK paper) benchmarked the algorithm on Muscle5
# stratified ensembles — same-algorithm parameter perturbations with
# similar column structures. CLOAK assumes that. Cross-algorithm
# ensembles (MAFFT + FAMSA) violate the assumption because progressive-
# vs-iterative aligners produce fundamentally different column counts.
# We get our diversity via MAFFT parameter perturbation instead, which
# is exactly what CLOAK was designed for.
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
run_alignment_ensemble.sh — generate ensemble of MAFFT variants for CLOAK

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

# --- mafft --dash: structural prior flag (locked decision 2026-05-28) ---
# When MAFFT_DASH=1 (default), prepend --dash to all MAFFT invocations.
# --dash uses PDB structural priors to guide alignment of GPCR TM helices.
# --originalseqonly is REQUIRED alongside --dash: it lets the DASH structural
# homologs GUIDE the alignment but EXCLUDES them from the output. Without it,
# MAFFT injects 'DASH|<pdb>_<chain>||...' homolog rows that pollute the per-class
# pool and tree (they broke the class-A IQ-TREE run: ~239 leaked DASH rows ->
# seed-tree/alignment name mismatch + OOM).
# FAMSA calls are not affected (FAMSA does not support --dash).
MAFFT_DASH_ARGS=""
if [[ "${MAFFT_DASH:-1}" == "1" ]]; then
    MAFFT_DASH_ARGS="--dash --originalseqonly"
fi

# --- Canonical alignment (regime-based via run_aligner.sh, OR a chosen mode) ---
case "$CANONICAL_ALIGNER" in
    auto)
        # run_aligner.sh respects MAFFT_DASH env var internally via the same mechanism
        bash "${SCRIPT_DIR}/run_aligner.sh" \
            --input="$INPUT" \
            --output="${OUTPUT_DIR}/canonical.fa" \
            --threads="$THREADS"
        ;;
    linsi) "$MAFFT_BIN" ${MAFFT_DASH_ARGS} --localpair  --maxiterate 1000 --thread "$THREADS" "$INPUT" > "${OUTPUT_DIR}/canonical.fa" 2>"${OUTPUT_DIR}/canonical.log" ;;
    ginsi) "$MAFFT_BIN" ${MAFFT_DASH_ARGS} --globalpair --maxiterate 1000 --thread "$THREADS" "$INPUT" > "${OUTPUT_DIR}/canonical.fa" 2>"${OUTPUT_DIR}/canonical.log" ;;
    einsi) "$MAFFT_BIN" ${MAFFT_DASH_ARGS} --genafpair  --maxiterate 1000 --thread "$THREADS" "$INPUT" > "${OUTPUT_DIR}/canonical.fa" 2>"${OUTPUT_DIR}/canonical.log" ;;
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
    # shellcheck disable=SC2086
    if "$MAFFT_BIN" ${MAFFT_DASH_ARGS} "$@" --thread "$THREADS" "$INPUT" > "$out" 2>"$log" \
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

# NOTE: FAMSA was deliberately removed from the ensemble (see header comment).
# It produces alignments with ~2x the column count of MAFFT, which breaks
# CLOAK's same-column-structure assumption. FAMSA is still available as the
# canonical aligner via --canonical-aligner=famsa for use cases outside the
# CLOAK ensemble (e.g. very large datasets where MAFFT would be too slow).

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
