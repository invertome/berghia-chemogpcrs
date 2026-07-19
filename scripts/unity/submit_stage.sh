#!/bin/bash
# submit_stage.sh — submit a numbered pipeline stage on Unity.
#
# Background: every numbered stage script (01_…sh, 02_…sh, …) carries
# `#SBATCH` directives with bash `${VAR}` substitutions and a
# `#SBATCH $(scale_resources)` command-substitution line. sbatch parses
# directives BEFORE running the script body, so it sees those tokens
# verbatim and rejects them ("Invalid --time specification", "Invalid
# directive found in batch script: $(scale_resources)").
#
# Workaround captured here once: source config.sh + functions.sh, then
# `envsubst` the script and strip the `$(…)` directive line. The result
# is a sbatch-clean copy under /tmp; we feed that to sbatch and forward
# any extra CLI flags (--partition, --gres, --dependency, …) through.
#
# Usage:
#   bash scripts/unity/submit_stage.sh STAGE_SCRIPT [extra sbatch flags...]
#
# Examples (from the workdir):
#   bash scripts/unity/submit_stage.sh 01_reference_processing.sh \
#        --time=24:00:00 --cpus-per-task=8 --mem=96G
#
#   bash scripts/unity/submit_stage.sh 02_chemogpcrs_identification.sh \
#        --partition=gpu --gres=gpu:1 --time=12:00:00 --mem=64G \
#        --cpus-per-task=8 --dependency=afterok:57426390
#
# The script itself is not a Slurm job — it just preprocesses + submits
# from the login node and prints the new job id on stdout.

set -eo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 STAGE_SCRIPT [extra sbatch flags...]" >&2
    exit 2
fi

STAGE_SCRIPT="$1"; shift
STAGE_BASENAME=$(basename "$STAGE_SCRIPT" .sh)

if [[ ! -f "$STAGE_SCRIPT" ]]; then
    echo "ERROR: stage script not found: $STAGE_SCRIPT" >&2
    exit 2
fi

# Source config + functions so that envsubst can substitute `${VAR}`
# tokens with their actual values.
# shellcheck disable=SC1091
source config.sh
# shellcheck disable=SC1091
source functions.sh

# Hand envsubst a complete environment — the SBATCH directives reference
# multiple project variables.
export DEFAULT_TIME LOGS_DIR CPUS DEFAULT_MEM RESULTS_DIR
export SLURM_EMAIL="${SLURM_EMAIL:-${USER}@umass.edu}"

OUT_DIR="${TMPDIR:-/tmp}/berghia_sbatch"
mkdir -p "$OUT_DIR"
PREPROCESSED="${OUT_DIR}/${STAGE_BASENAME}.sbatch.sh"

# envsubst handles ${VAR}; sed strips any directive that's a $(...)
# command substitution since sbatch can't run those at parse time.
# Resource flags should be passed via CLI overrides instead.
envsubst < "$STAGE_SCRIPT" | sed '/^#SBATCH \\?\$(/d' > "$PREPROCESSED"

echo "[submit_stage] preprocessed -> $PREPROCESSED" >&2
echo "[submit_stage] sbatch flags: $*" >&2

# --- Bead fxx: size --array from the manifest, not the hardcoded range --------
# Stages 04/05 carry `#SBATCH --array=0-999%50`, which cannot cover a full-557
# manifest (orthogroups at index >=1000 fall outside the array). sbatch CLI
# flags take precedence over #SBATCH directives, so compute the real range from
# OG_COUNT and inject it. Skipped when the caller passes an explicit --array.
ARRAY_FLAG=()
if grep -qE '^#SBATCH[[:space:]]+--array=' "$STAGE_SCRIPT" \
   && ! printf '%s\n' "$@" | grep -q -- '--array='; then
    _manifest="${RESULTS_DIR}/orthogroup_manifest.tsv"
    _og_count=$(get_orthogroup_count "$_manifest")
    # Preserve the stage's intended concurrency limit (%K) unless overridden.
    _throttle="${ARRAY_THROTTLE:-$(grep -oP '^#SBATCH\s+--array=\S*%\K\d+' "$STAGE_SCRIPT" | head -1)}"
    _throttle="${_throttle:-50}"
    _detected_max=$(scontrol show config 2>/dev/null | grep -oP 'MaxArraySize\s*=\s*\K\d+' | head -1)
    _max_array="${MAX_ARRAY_SIZE:-${_detected_max:-10001}}"

    if [ "${_og_count:-0}" -gt 0 ] 2>/dev/null; then
        if _spec=$(python3 "${SCRIPTS_DIR}/array_plan.py" --n-tasks "$_og_count" \
                       --throttle "$_throttle" --max-array-size "$_max_array"); then
            ARRAY_FLAG=(--array="$_spec")
            echo "[submit_stage] sized array from manifest: --array=${_spec} (OG_COUNT=${_og_count}, MaxArraySize=${_max_array})" >&2
        else
            echo "[submit_stage] ERROR: manifest (${_og_count} orthogroups) exceeds MaxArraySize=${_max_array}; refusing to submit a truncated array." >&2
            exit 3
        fi
    else
        echo "[submit_stage] WARN: orthogroup count unknown/zero (${_manifest}); leaving the script's --array directive in place" >&2
    fi
fi

# Submit and surface the job id on stdout for capture by callers.
sbatch_out=$(sbatch "${ARRAY_FLAG[@]}" "$@" "$PREPROCESSED")
echo "$sbatch_out"
echo "$sbatch_out" | grep -oP 'Submitted batch job \K\d+'
