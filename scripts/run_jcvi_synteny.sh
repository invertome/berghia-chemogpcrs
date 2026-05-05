#!/bin/bash
# run_jcvi_synteny.sh - Bead -e59
# Wrapper around JCVI MCscan (Tang et al. 2024 iMeta 3:e211) for pairwise
# synteny between the Berghia genome and a close-relative mollusc genome.
#
# This is a thin orchestration layer:
#   - put Berghia and target proteomes/GFFs into JCVI's expected naming
#     convention (<prefix>.bed and <prefix>.pep co-located in workdir),
#   - call `python -m jcvi.compara.catalog ortholog`,
#   - report which target was used (provenance file).
#
# The actual anchor detection is done entirely by JCVI; we do NOT roll our
# own. See https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version).
#
# Usage:
#   bash run_jcvi_synteny.sh \
#     --berghia-prefix=<prefix> \
#     --berghia-gff=<file> --berghia-pep=<file> \
#     --target=<species_proteome.fa> --target-gff=<file> \
#     [--target-prefix=<prefix>] \
#     --output-dir=<dir>

set -uo pipefail

# Force matplotlib's headless backend so JCVI's dotplot generation
# (jcvi.graphics.dotplot, called during compara.catalog ortholog) does not
# attempt to import a display backend on Unity compute nodes. Without this
# the ortholog step crashes with "_tkinter.TclError: no display name and no
# $DISPLAY environment variable" — which would abort the whole stage even
# though the scientifically-relevant anchors file is produced before the
# dotplot step.
export MPLBACKEND=Agg

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# ---- parse args ----
BERGHIA_PREFIX=""
BERGHIA_GFF=""
BERGHIA_PEP=""
TARGET=""
TARGET_GFF=""
TARGET_PREFIX=""
OUTPUT_DIR=""
CSCORE="${JCVI_CSCORE:-0.7}"
N_HITS="${JCVI_N_HITS:-1}"

usage() {
    cat >&2 <<EOF
run_jcvi_synteny.sh — JCVI MCscan wrapper (bead -e59)

Required arguments:
  --berghia-prefix=<str>   Short name used for Berghia files (e.g. berghia)
  --berghia-gff=<path>     Berghia GFF3 (e.g. \$GENOME_GFF)
  --berghia-pep=<path>     Berghia protein FASTA (e.g. \$GENOME_PROTEIN)
  --target=<path>          Target species protein FASTA
  --target-gff=<path>      Target species GFF3
  --output-dir=<path>      Working / output directory

Optional:
  --target-prefix=<str>    Override target prefix (default: basename of --target)
  --cscore=<float>         JCVI C-score threshold (default 0.7; \$JCVI_CSCORE)
  --n-hits=<int>           JCVI top-N reciprocal best hits (default 1; \$JCVI_N_HITS)

Outputs (in --output-dir):
  <berghia>.<target>.anchors           collinear anchor blocks
  <berghia>.<target>.lifted.anchors    anchors lifted across blocks (denser)
  predictor_used.txt                   tool + version + target species
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --berghia-prefix=*) BERGHIA_PREFIX="${1#*=}"; shift ;;
        --berghia-gff=*)    BERGHIA_GFF="${1#*=}"; shift ;;
        --berghia-pep=*)    BERGHIA_PEP="${1#*=}"; shift ;;
        --target=*)         TARGET="${1#*=}"; shift ;;
        --target-gff=*)     TARGET_GFF="${1#*=}"; shift ;;
        --target-prefix=*)  TARGET_PREFIX="${1#*=}"; shift ;;
        --output-dir=*)     OUTPUT_DIR="${1#*=}"; shift ;;
        --cscore=*)         CSCORE="${1#*=}"; shift ;;
        --n-hits=*)         N_HITS="${1#*=}"; shift ;;
        -h|--help)          usage; exit 0 ;;
        *)                  echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

if [ -z "$BERGHIA_PREFIX" ] || [ -z "$BERGHIA_GFF" ] || [ -z "$BERGHIA_PEP" ] \
   || [ -z "$TARGET" ] || [ -z "$TARGET_GFF" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: missing required argument(s)" >&2
    usage
    exit 1
fi

for f in "$BERGHIA_GFF" "$BERGHIA_PEP" "$TARGET" "$TARGET_GFF"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: input file not found: $f" >&2
        exit 1
    fi
done

if [ -z "$TARGET_PREFIX" ]; then
    TARGET_PREFIX=$(basename "$TARGET")
    TARGET_PREFIX="${TARGET_PREFIX%.*}"
    # Strip a trailing .proteins, .pep, .faa if present
    TARGET_PREFIX="${TARGET_PREFIX%.proteins}"
    TARGET_PREFIX="${TARGET_PREFIX%.pep}"
    TARGET_PREFIX="${TARGET_PREFIX%.faa}"
fi

mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

PYTHON_BIN="${JCVI_PYTHON:-python3}"

# ---- check JCVI is importable ----
if ! "$PYTHON_BIN" -c 'import jcvi' 2>/dev/null; then
    echo "ERROR: jcvi is not installed for $PYTHON_BIN. Run: pip install jcvi" >&2
    exit 2
fi

JCVI_VERSION=$("$PYTHON_BIN" -c 'import jcvi, sys; sys.stdout.write(getattr(jcvi, "__version__", "unknown"))' 2>/dev/null || echo "unknown")

# ---- copy inputs into workdir under JCVI naming convention ----
# JCVI expects <prefix>.bed and <prefix>.pep (or .cds) co-located in workdir.
cd "$OUTPUT_DIR"

# Build Berghia BED from GFF.
# `jcvi.formats.gff bed` reads a GFF3 and writes a BED with one line per
# `--type` feature, using `--key` as the gene name. mRNA + Name is the
# RefSeq convention (GCF_034508935.2 has both `gene` and `mRNA` features
# and uses the Name= attribute for stable IDs).
"$PYTHON_BIN" -m jcvi.formats.gff bed \
    --type=mRNA --key=Name --primary_only \
    "$BERGHIA_GFF" -o "${BERGHIA_PREFIX}.bed" \
    2> "${OUTPUT_DIR}/jcvi_bed_${BERGHIA_PREFIX}.log"
JCVI_RC=$?
if [ $JCVI_RC -ne 0 ] || [ ! -s "${BERGHIA_PREFIX}.bed" ]; then
    echo "ERROR: failed to build BED from $BERGHIA_GFF (exit=$JCVI_RC)" >&2
    cat "${OUTPUT_DIR}/jcvi_bed_${BERGHIA_PREFIX}.log" >&2 || true
    exit 3
fi

"$PYTHON_BIN" -m jcvi.formats.gff bed \
    --type=mRNA --key=Name --primary_only \
    "$TARGET_GFF" -o "${TARGET_PREFIX}.bed" \
    2> "${OUTPUT_DIR}/jcvi_bed_${TARGET_PREFIX}.log"
JCVI_RC=$?
if [ $JCVI_RC -ne 0 ] || [ ! -s "${TARGET_PREFIX}.bed" ]; then
    echo "ERROR: failed to build BED from $TARGET_GFF (exit=$JCVI_RC)" >&2
    cat "${OUTPUT_DIR}/jcvi_bed_${TARGET_PREFIX}.log" >&2 || true
    exit 3
fi

# Stage protein FASTAs as <prefix>.pep (JCVI looks for .pep / .cds beside .bed).
# Using cp -L to dereference symlinks (genomes/ has symlink chains).
cp -L "$BERGHIA_PEP" "${BERGHIA_PREFIX}.pep"
cp -L "$TARGET"     "${TARGET_PREFIX}.pep"

# ---- run jcvi.compara.catalog ortholog ----
# This runs the LAST (or BLAST) all-vs-all, scores hits with the
# C-score, and writes <a>.<b>.anchors / .lifted.anchors / .last / .last.filtered.
ANCHORS_OUT="${BERGHIA_PREFIX}.${TARGET_PREFIX}.anchors"
LIFTED_OUT="${BERGHIA_PREFIX}.${TARGET_PREFIX}.lifted.anchors"

"$PYTHON_BIN" -m jcvi.compara.catalog ortholog \
    --no_strip_names \
    --cscore="$CSCORE" \
    "$BERGHIA_PREFIX" "$TARGET_PREFIX" \
    > "${OUTPUT_DIR}/jcvi_ortholog.log" 2>&1
ORTH_RC=$?

if [ $ORTH_RC -ne 0 ] && [ ! -s "$ANCHORS_OUT" ]; then
    echo "ERROR: jcvi.compara.catalog ortholog failed (exit=$ORTH_RC) with no anchor output." >&2
    tail -40 "${OUTPUT_DIR}/jcvi_ortholog.log" >&2 || true
    exit 4
fi

# Anchors exist but the ortholog step exited non-zero — this typically means
# the dotplot rendering failed (e.g. matplotlib backend issue) AFTER the
# anchors were produced. Log a warning and continue: the report-image hook in
# stage 09 will simply skip the missing PDF, but the scientific signal (the
# anchor file feeding the synteny axis of ranking) is intact.
if [ $ORTH_RC -ne 0 ] && [ -s "$ANCHORS_OUT" ]; then
    echo "WARNING: jcvi.compara.catalog ortholog returned non-zero (exit=$ORTH_RC) but anchors are present;" >&2
    echo "         likely a dotplot-rendering failure. Continuing." >&2
    tail -10 "${OUTPUT_DIR}/jcvi_ortholog.log" >&2 || true
fi

if [ ! -s "$ANCHORS_OUT" ]; then
    echo "WARNING: anchor file is empty (no synteny detected): $ANCHORS_OUT" >&2
fi

# ---- record provenance ----
PROV_FILE="${OUTPUT_DIR}/predictor_used.txt"
{
    echo "tool=jcvi.compara.catalog.ortholog"
    echo "jcvi_version=${JCVI_VERSION}"
    echo "berghia_prefix=${BERGHIA_PREFIX}"
    echo "target_prefix=${TARGET_PREFIX}"
    echo "target_proteome=${TARGET}"
    echo "target_gff=${TARGET_GFF}"
    echo "cscore=${CSCORE}"
    echo "n_hits=${N_HITS}"
    echo "anchors_file=${ANCHORS_OUT}"
    echo "lifted_anchors_file=${LIFTED_OUT}"
    echo "run_at=$(date -Iseconds)"
} > "$PROV_FILE"

echo "JCVI synteny complete: $ANCHORS_OUT (provenance: $PROV_FILE)" >&2
exit 0
