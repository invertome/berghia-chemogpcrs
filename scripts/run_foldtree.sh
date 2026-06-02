#!/bin/bash
# run_foldtree.sh — adapter: stage-08 CLI -> DessimozLab/fold_tree CLI (bead -68w).
#
# Stage 08 invokes:  ${FOLDTREE} --input_dir <dir> --output <tre> --method <m>
# but the real tool (https://github.com/DessimozLab/fold_tree) is a Snakemake
# pipeline exposing:  foldtree --folder <dir with structs/ + identifiers.txt> --custom-structs
#
# This adapter reconciles them:
#   * maps --method to a fold_tree metric (mattypes = foldtree|alntmscore|lddt);
#     stage 08's default "upgma" -> "foldtree" (the structural metric);
#   * stages structures into <run>/structs/ as .pdb — fold_tree's structs2fasta.py
#     globs structs/*.pdb ONLY, so AF3 mmCIF (*.cif) is converted via `gemmi convert`;
#   * writes an empty identifiers.txt (required by --custom-structs);
#   * runs foldtree in its dedicated env with --no-filter (the input mixes AF3
#     predictions with experimental GPCRdb structures whose B-factors are NOT
#     pLDDT, so pLDDT filtering would wrongly drop them);
#   * copies the selected metric tree to --output.
#
# Env vars: FOLDTREE_ENV (default foldtree), FOLDTREE_CONDA_PREFIX (default
# ~/.foldtree), FOLDTREE_RUN_DIR (default <output dir>/foldtree_run), CPUS.

set -eo pipefail

INPUT_DIR=""; OUTPUT=""; METHOD="foldtree"; DRY_RUN=0
while [ $# -gt 0 ]; do
    case "$1" in
        --input_dir) INPUT_DIR="$2"; shift 2;;
        --output)    OUTPUT="$2";    shift 2;;
        --method)    METHOD="$2";    shift 2;;
        --dry-run)   DRY_RUN=1;      shift;;
        *) echo "run_foldtree: unknown arg: $1" >&2; exit 2;;
    esac
done
[ -n "$INPUT_DIR" ] && [ -n "$OUTPUT" ] || {
    echo "run_foldtree: --input_dir and --output are required" >&2; exit 2; }
[ -d "$INPUT_DIR" ] || { echo "run_foldtree: input dir not found: $INPUT_DIR" >&2; exit 2; }

# --method -> fold_tree metric (mattypes from the workflow: foldtree|alntmscore|lddt)
case "$(printf '%s' "$METHOD" | tr '[:upper:]' '[:lower:]')" in
    foldtree|structural|upgma|nj|fident) METRIC="foldtree" ;;
    tm|tmscore|alntmscore)               METRIC="alntmscore" ;;
    lddt)                                METRIC="lddt" ;;
    *)                                   METRIC="foldtree" ;;
esac

CORES="${CPUS:-4}"
FOLDTREE_ENV="${FOLDTREE_ENV:-foldtree}"
RUN_DIR="${FOLDTREE_RUN_DIR:-$(dirname "$OUTPUT")/foldtree_run}"
CONDA_PREFIX_FT="${FOLDTREE_CONDA_PREFIX:-$HOME/.foldtree}"

# Fresh run dir: structs/ + empty identifiers.txt (required by --custom-structs).
rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR/structs"
: > "$RUN_DIR/identifiers.txt"

# Activate the foldtree env lazily — only when a CIF needs gemmi, or for the
# real run. Keeps the .pdb-only dry-run path conda-free (and unit-testable).
FT_ENV_ACTIVE=0
ensure_ft_env() {
    [ "$FT_ENV_ACTIVE" = 1 ] && return 0
    # shellcheck disable=SC1091
    source ~/.miniconda3/etc/profile.d/conda.sh
    conda activate "$FOLDTREE_ENV"
    FT_ENV_ACTIVE=1
}

n_pdb=0
shopt -s nullglob
for f in "$INPUT_DIR"/*; do
    [ -f "$f" ] || continue
    bn="$(basename "$f")"; stem="${bn%.*}"
    case "$bn" in
        *.pdb)
            cp "$f" "$RUN_DIR/structs/${stem}.pdb"; n_pdb=$((n_pdb+1)) ;;
        *.cif|*.mmcif)
            ensure_ft_env
            if gemmi convert "$f" "$RUN_DIR/structs/${stem}.pdb"; then
                n_pdb=$((n_pdb+1))
            else
                echo "run_foldtree: WARN cif->pdb failed for $bn (skipped)" >&2
            fi ;;
        *) : ;;
    esac
done
shopt -u nullglob
echo "run_foldtree: prepared $n_pdb pdb structures in $RUN_DIR/structs (metric=$METRIC)"

FT_CMD=(foldtree --folder "$RUN_DIR" --custom-structs --cores "$CORES"
        --no-filter --conda-prefix "$CONDA_PREFIX_FT" -p)

if [ "$DRY_RUN" = 1 ]; then
    echo "DRY: metric=$METRIC structs=$n_pdb"
    echo "DRY: ${FT_CMD[*]}"
    exit 0
fi

[ "$n_pdb" -ge 3 ] || echo "run_foldtree: WARN only $n_pdb structures (<3) — tree may be trivial" >&2

ensure_ft_env
"${FT_CMD[@]}"

# Select the metric's tree (fallback chain if a late post-processing step failed).
SEL=""
for cand in \
    "$RUN_DIR/${METRIC}_struct_tree.PP.nwk.rooted.final" \
    "$RUN_DIR/${METRIC}_struct_tree.PP.nwk.rooted" \
    "$RUN_DIR/${METRIC}_struct_tree.PP.nwk" \
    "$RUN_DIR/${METRIC}_struct_tree.nwk"; do
    [ -s "$cand" ] && { SEL="$cand"; break; }
done
[ -n "$SEL" ] || { echo "run_foldtree: no ${METRIC} tree produced in $RUN_DIR" >&2; exit 1; }

mkdir -p "$(dirname "$OUTPUT")"
cp "$SEL" "$OUTPUT"
echo "run_foldtree: $(basename "$SEL") -> $OUTPUT"
