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

# Project variables the submitted job should inherit in its environment
# (sbatch propagates the submitting environment by default). This is NOT the
# envsubst substitution set. See the SHELL-FORMAT derivation below, which is
# computed from the stage's own #SBATCH directives.
export DEFAULT_TIME LOGS_DIR CPUS DEFAULT_MEM RESULTS_DIR
export SLURM_EMAIL="${SLURM_EMAIL:-${USER}@umass.edu}"

OUT_DIR="${TMPDIR:-/tmp}/berghia_sbatch"
mkdir -p "$OUT_DIR"
PREPROCESSED="${OUT_DIR}/${STAGE_BASENAME}.sbatch.sh"

# --- Bead mqme: substitute ONLY the variables the #SBATCH directives use ------
# `envsubst` with NO SHELL-FORMAT argument substitutes EVERY $VAR/${VAR} in the
# whole file and replaces unset ones with the empty string. That guts the stage
# BODY: loop variables, locals, and (decisively) $SLURM_ARRAY_TASK_ID, which
# Slurm sets at RUNTIME and is therefore necessarily unset here. Stage 04's
# `if [ -n "$SLURM_ARRAY_TASK_ID" ]` became `if [ -n "" ]`, false for every array
# task, so every task skipped all per-orthogroup work and exited 0. Measured on
# stage 04: 224 corrupted lines.
#
# So derive the substitution set from the stage's OWN #SBATCH directives and
# hand it to envsubst explicitly; every other token in the file passes through
# untouched. `$(scale_resources)` is not a variable reference and never matches.
mapfile -t SBATCH_VARS < <(
    grep '^#SBATCH' "$STAGE_SCRIPT" \
        | grep -oE '\$\{[A-Za-z_][A-Za-z0-9_]*\}|\$[A-Za-z_][A-Za-z0-9_]*' \
        | tr -d '${}' \
        | sort -u
)

# Each name must resolve to a NON-EMPTY value. An unset one would reintroduce
# exactly the silent-empty-substitution failure this fix exists to remove, only
# now inside a resource directive (`--time=` with no value), so fail loudly.
SHELL_FORMAT=""
for _var in "${SBATCH_VARS[@]:-}"; do
    [ -n "$_var" ] || continue
    if [ -z "${!_var:-}" ]; then
        echo "ERROR: #SBATCH directive in ${STAGE_SCRIPT} references \${${_var}}," \
             "which is unset or empty after sourcing config.sh + functions.sh." >&2
        echo "       Refusing to substitute it with an empty string." >&2
        exit 4
    fi
    export "${_var?}"
    SHELL_FORMAT+="\$${_var} "
done
echo "[submit_stage] envsubst shell-format: ${SHELL_FORMAT:-<none>}" >&2

# Substitution is applied ONLY to `#SBATCH` directive lines. The script body is
# copied through byte-for-byte, so no body construct can be rewritten even if it
# happens to mention one of the directive variables. Directives that are a
# $(...) command substitution are dropped, since sbatch parses directives before
# running the script and cannot evaluate those; pass resource flags via CLI
# overrides instead.
{
    while IFS= read -r _line || [ -n "$_line" ]; do
        if [[ "$_line" =~ ^#SBATCH[[:space:]] ]]; then
            [[ "$_line" =~ ^#SBATCH[[:space:]]+\$\( ]] && continue
            printf '%s\n' "$_line" | envsubst "$SHELL_FORMAT"
        else
            printf '%s\n' "$_line"
        fi
    done < "$STAGE_SCRIPT"
} > "$PREPROCESSED"

# Nothing may reach sbatch with an unresolved reference in a directive: sbatch
# would either reject the file or silently take a malformed value.
if grep -nE '^#SBATCH.*\$' "$PREPROCESSED" >&2; then
    echo "ERROR: unresolved substitution left in a #SBATCH directive (above)." >&2
    exit 4
fi

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
