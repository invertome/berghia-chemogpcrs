#!/bin/bash
# run_alphafold.sh — AlphaFold 3 wrapper for stage 08 (bead -05o).
#
# Accepts AF2-style args (--fasta_paths, --output_dir) for backwards
# compatibility with the stage 08 call site, then internally:
#   1. Loads the alphafold3 Lmod module on Unity (Apptainer-based)
#   2. Converts the input FASTA to AF3-required JSON via `af3 convert`
#   3. Runs `af3 run` (which wraps run_alphafold.py inside the SIF) with
#      the correct --json_path, --output_dir, --model_dir, --db_dir args
#
# Inputs (AF2 args, preserved from the legacy stage 08 invocation):
#   --fasta_paths=PATH        Path to the input FASTA (required)
#   --output_dir=PATH         Output directory (required)
#   --max_template_date=...   IGNORED (AF3 uses its own date cutoff)
#   --use_gpu=...             IGNORED (AF3 requires GPU; aborts if absent)
#   --model_preset=...        IGNORED (AF3 only has one preset)
#
# Required env vars (defaults from config.sh):
#   ALPHAFOLD3_DB_DIR    — sequence DBs (default /datasets/bio/alphafold3)
#   ALPHAFOLD3_MODEL_DIR — AF3 weights (NO DEFAULT; Google-gated)
#                          Request: https://forms.gle/svvpY4u2jsHEwWYS6
#
# Why Apptainer (not pip-install): the alphafold3 Lmod module on Unity
# defines `run_alphafold.py` and `af3` as shell functions that invoke
# `apptainer exec --nv unity-alphafold_latest.sif ...`. The container
# bundles JAX, CUDA libs, and run_alphafold.py — there is no host-side
# python install of alphafold3 to invoke directly.
#
# Bead -05o, 2026-05-20.

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

FASTA=""
OUTPUT_DIR=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta_paths=*)        FASTA="${1#*=}"; shift ;;
        --fasta_paths)          FASTA="$2"; shift 2 ;;
        --output_dir=*)         OUTPUT_DIR="${1#*=}"; shift ;;
        --output_dir)           OUTPUT_DIR="$2"; shift 2 ;;
        --max_template_date=*|--use_gpu=*|--model_preset=*)
            # AF2-era flags, no AF3 equivalent — silently ignored.
            shift ;;
        --max_template_date|--use_gpu|--model_preset)
            shift 2 ;;
        *)
            echo "[run_alphafold] WARN: unrecognized arg '$1' (ignored)" >&2
            shift ;;
    esac
done

if [ -z "$FASTA" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 --fasta_paths=INPUT.fasta --output_dir=OUT_DIR" >&2
    exit 2
fi
[ -s "$FASTA" ] || { echo "[run_alphafold] ERROR: FASTA missing or empty: $FASTA" >&2; exit 2; }
mkdir -p "$OUTPUT_DIR"

# --- Database + model paths (sourced from config.sh; allow direct env override) ---
DB_DIR="${ALPHAFOLD3_DB_DIR:-/datasets/bio/alphafold3}"
MODEL_DIR="${ALPHAFOLD3_MODEL_DIR:-}"

if [ -z "$MODEL_DIR" ] || [ ! -d "$MODEL_DIR" ]; then
    cat >&2 <<EOF
[run_alphafold] ERROR: AlphaFold 3 model weights not configured.

  ALPHAFOLD3_MODEL_DIR is unset or not a directory: '$MODEL_DIR'

  AF3 weights are Google-gated and not bundled with the Unity dataset
  ($DB_DIR contains only sequence/structure databases).

  To enable stage 08:
    1. Request weights via https://forms.gle/svvpY4u2jsHEwWYS6
       (review the terms of use:
        https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)
    2. Once approved, download to a Unity scratch path, e.g.:
        /scratch3/workspace/\$USER-jorge/alphafold3_models/
    3. Export ALPHAFOLD3_MODEL_DIR=<that path> in config.sh or your env.
EOF
    exit 3
fi
[ -d "$DB_DIR" ] || { echo "[run_alphafold] ERROR: DB_DIR missing: $DB_DIR" >&2; exit 3; }

# --- GPU check ---
# AF3 runs only on GPU; CPU inference is impractical (hours per protein).
# If we're in a SLURM job, verify a GPU was actually allocated.
GPU_ALLOCATED=0
[ -n "${SLURM_JOB_GPUS:-}${SLURM_GPUS:-}${SLURM_GPUS_ON_NODE:-}${CUDA_VISIBLE_DEVICES:-}" ] && GPU_ALLOCATED=1
if [ "$GPU_ALLOCATED" != "1" ] && [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "[run_alphafold] ERROR: SLURM job has no GPU allocation. AF3 requires GPU." >&2
    echo "                Submit stage 08 with --partition=gpu --gres=gpu:1" >&2
    exit 4
fi

# --- Load alphafold3 module ---
# The module is Lmod and only available in login-shell context; source the
# Unity profile so `module` is defined inside this sbatch step.
if ! command -v module &>/dev/null; then
    # shellcheck disable=SC1091
    [ -r /etc/profile.d/lmod.sh ] && source /etc/profile.d/lmod.sh
    [ -r /etc/profile ] && source /etc/profile
fi

if ! command -v module &>/dev/null; then
    echo "[run_alphafold] ERROR: cannot find Lmod 'module' command after sourcing /etc/profile*." >&2
    echo "                Are you running on Unity?" >&2
    exit 5
fi

# Apptainer cache dir — Unity admins explicitly warn against the default
# (~/.apptainer/cache) which fills home quotas.
mkdir -p "$APPTAINER_CACHEDIR" 2>/dev/null || true
export APPTAINER_CACHEDIR

# shellcheck disable=SC1091
module load alphafold3/latest 2>&1 | grep -v '^$' >&2 || true

# Sanity check: `af3` should now be defined as a function.
if ! type af3 >/dev/null 2>&1; then
    echo "[run_alphafold] ERROR: 'af3' command not found after 'module load alphafold3/latest'." >&2
    exit 5
fi

# --- 1) Convert FASTA → AF3 JSON ---
JSON_DIR="$OUTPUT_DIR/_af3_input_json"
mkdir -p "$JSON_DIR"
echo "[run_alphafold] Converting FASTA → JSON via 'af3 convert'..." >&2
af3 convert --fasta "$FASTA" --output_dir "$JSON_DIR"

# `af3 convert` writes one JSON per sequence. Detect them.
shopt -s nullglob
JSON_FILES=("$JSON_DIR"/*.json)
shopt -u nullglob
if [ "${#JSON_FILES[@]}" -eq 0 ]; then
    echo "[run_alphafold] ERROR: 'af3 convert' produced no JSON files in $JSON_DIR" >&2
    exit 6
fi
echo "[run_alphafold] Converted ${#JSON_FILES[@]} sequences → JSON" >&2

# --- 2) Run AF3 per JSON ---
for json in "${JSON_FILES[@]}"; do
    name="$(basename "$json" .json)"
    sub_out="$OUTPUT_DIR/$name"
    if [ -d "$sub_out" ] && [ -f "$sub_out/ranking_scores.csv" ]; then
        echo "[run_alphafold] $name: skipping (output exists at $sub_out)" >&2
        continue
    fi
    mkdir -p "$sub_out"
    echo "[run_alphafold] $name: running AF3..." >&2
    af3 run \
        --json_path "$json" \
        --output_dir "$sub_out" \
        --model_dir "$MODEL_DIR" \
        --db_dir "$DB_DIR"
done

echo "[run_alphafold] Done. Outputs under $OUTPUT_DIR" >&2
