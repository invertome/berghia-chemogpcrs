#!/bin/bash
# run_treeshrink.sh — TreeShrink (Mai & Mirarab 2018, BMC Genomics 19:272)
#
# Bead -iof. Removes outlier-long branches from per-OG and global gene
# trees before downstream selection analysis. Important for paralog-rich
# Class A GPCR families where some leaves are rogue artefacts.
#
# Wraps the upstream `run_treeshrink.py` shim (installed via
# `pip install treeshrink`); does NOT roll our own outlier detection.
#
# Modes:
#   --single-tree=<file>    Apply to a single newick tree
#   --input-dir=<dir>       Apply to every .treefile in a directory
#
# Outputs (per input tree):
#   <tree>.original.treefile   the input copied verbatim (rollback)
#   <tree>                     replaced with the cleaned (shrunk) tree
#   <out_dir>/<base>/output.* TreeShrink working files (kept for debugging)

set -uo pipefail

SINGLE_TREE=""
INPUT_DIR=""
OUTPUT_DIR=""
QUANTILE="${TREESHRINK_QUANTILE:-0.05}"
TREESHRINK_BIN="${TREESHRINK:-run_treeshrink.py}"

usage() {
    cat >&2 <<EOF
run_treeshrink.sh — TreeShrink wrapper (bead -iof)

Mode (one of):
  --single-tree=<file>   Apply TreeShrink to a single tree
  --input-dir=<dir>      Apply to every *.treefile in a directory

Required:
  --output-dir=<path>    Working / output directory

Optional:
  --quantile=<float>     Per-clade quantile threshold (default 0.05; \$TREESHRINK_QUANTILE)

Tool:
  TREESHRINK env var picks the binary (default: run_treeshrink.py).
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --single-tree=*) SINGLE_TREE="${1#*=}"; shift ;;
        --input-dir=*)   INPUT_DIR="${1#*=}"; shift ;;
        --output-dir=*)  OUTPUT_DIR="${1#*=}"; shift ;;
        --quantile=*)    QUANTILE="${1#*=}"; shift ;;
        -h|--help)       usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

if [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: --output-dir is required" >&2
    usage
    exit 1
fi
if [ -z "$SINGLE_TREE" ] && [ -z "$INPUT_DIR" ]; then
    echo "ERROR: must provide --single-tree=<file> or --input-dir=<dir>" >&2
    usage
    exit 1
fi

if ! command -v "$TREESHRINK_BIN" &>/dev/null; then
    echo "ERROR: TreeShrink not found ($TREESHRINK_BIN). Install via 'pip install treeshrink'." >&2
    exit 2
fi

mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

shrink_one() {
    local tree="$1"
    if [ ! -f "$tree" ]; then
        echo "WARN: tree not found: $tree" >&2
        return 1
    fi
    local base
    base=$(basename "$tree" .treefile)
    base="${base%.tree}"   # also strip .tree if user named it that way
    local workdir="$OUTPUT_DIR/$base"
    mkdir -p "$workdir"

    # Copy input verbatim into workdir so rollback is always possible.
    cp -L "$tree" "$workdir/input.tree"

    "$TREESHRINK_BIN" \
        -t "$workdir/input.tree" \
        -q "$QUANTILE" \
        -o "$workdir" \
        > "$workdir/treeshrink.log" 2>&1
    local rc=$?

    local shrunk_tree
    # TreeShrink writes <input_basename>_treeshrink/output.tree by default
    shrunk_tree=$(find "$workdir" -name "output.tree" -print -quit 2>/dev/null)
    if [ -z "$shrunk_tree" ] || [ ! -s "$shrunk_tree" ]; then
        # Try alternate path: workdir/input_treeshrink/output.tree
        shrunk_tree="$workdir/input_treeshrink/output.tree"
        [ -f "$shrunk_tree" ] || shrunk_tree=""
    fi

    if [ $rc -ne 0 ] || [ -z "$shrunk_tree" ]; then
        echo "WARN: TreeShrink failed for $tree (rc=$rc); leaving original in place." >&2
        tail -5 "$workdir/treeshrink.log" >&2 || true
        return 1
    fi

    # Save the original alongside (for rollback) and replace the canonical
    # tree path with the shrunk tree.
    cp -f "$tree" "${tree}.original.treefile"
    cp -f "$shrunk_tree" "$tree"
    echo "TreeShrink: $tree (q=$QUANTILE; original at ${tree}.original.treefile)" >&2
    return 0
}

if [ -n "$SINGLE_TREE" ]; then
    shrink_one "$SINGLE_TREE" || exit $?
else
    n_total=0
    n_ok=0
    for tree in "$INPUT_DIR"/*.treefile; do
        [ -f "$tree" ] || continue
        n_total=$((n_total + 1))
        if shrink_one "$tree"; then
            n_ok=$((n_ok + 1))
        fi
    done
    echo "TreeShrink: $n_ok/$n_total trees cleaned" >&2
fi

# Provenance summary
PROV_FILE="$OUTPUT_DIR/predictor_used.txt"
{
    echo "tool=treeshrink"
    echo "binary=${TREESHRINK_BIN}"
    echo "quantile=${QUANTILE}"
    echo "run_at=$(date -Iseconds)"
} > "$PROV_FILE"
exit 0
