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

# Submit and surface the job id on stdout for capture by callers.
sbatch_out=$(sbatch "$@" "$PREPROCESSED")
echo "$sbatch_out"
echo "$sbatch_out" | grep -oP 'Submitted batch job \K\d+'
