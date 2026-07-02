#!/bin/bash
# Unit test for berghia_candidate_fasta() — the genome-track-aware selection of
# the Berghia candidate FASTA that stage 03 (orthology clustering) consumes.
#
# Contract (Task 7, genome-track reconciliation):
#   (i)   RUN_GENOME_TRACK=1 + reconciled file present -> reconciled path
#   (ii)  RUN_GENOME_TRACK=1 + reconciled file absent  -> legacy path
#   (iii) RUN_GENOME_TRACK=0                            -> legacy path
#         (regression guard: byte-identical to pre-genome-track behavior)
#
# Hermetic: RESULTS_DIR is a throwaway temp dir; the reconciled file's presence
# is toggled with touch/rm. berghia_candidate_fasta only reads RESULTS_DIR and
# RUN_GENOME_TRACK and only echoes a path, so no heavy binaries are invoked.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

# 1. Syntax of the helper's home + the consumer.
bash -n "$PROJECT_DIR/functions.sh"
bash -n "$PROJECT_DIR/03_orthology_clustering.sh"
echo "OK: functions.sh + 03_orthology_clustering.sh parse"

TMP="$(mktemp -d "${TMPDIR:-/tmp}/test_candsrc_XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

export RESULTS_DIR="$TMP/results"
mkdir -p "$RESULTS_DIR/reconciliation" "$RESULTS_DIR/chemogpcrs"

RECONCILED="$RESULTS_DIR/reconciliation/reconciled_candidates.faa"
LEGACY="$RESULTS_DIR/chemogpcrs/chemogpcrs_berghia.fa"

# Evaluate berghia_candidate_fasta from functions.sh with a given
# RUN_GENOME_TRACK, in an isolated subshell. Sourcing functions.sh installs an
# EXIT trap (finalize_pipeline); we disarm it (trap - EXIT) so only the helper's
# echoed path reaches stdout and this test's own cleanup trap is never clobbered.
# Pass a toggle value ("0"/"1"), or the sentinel "unset" to genuinely leave
# RUN_GENOME_TRACK undefined (exercises the ${RUN_GENOME_TRACK:-1} default).
eval_helper() {
    (
        if [ "$1" = "unset" ]; then
            unset RUN_GENOME_TRACK
        else
            export RUN_GENOME_TRACK="$1"
        fi
        # shellcheck disable=SC1091
        source "$PROJECT_DIR/functions.sh"
        trap - EXIT
        berghia_candidate_fasta
    )
}

assert_eq() {
    local label="$1" expected="$2" got="$3"
    if [ "$got" != "$expected" ]; then
        echo "FAIL: $label" >&2
        echo "  expected: $expected" >&2
        echo "  got:      $got" >&2
        exit 1
    fi
    echo "OK: $label"
}

# (i) genome track ON + reconciled file present -> reconciled path
: > "$RECONCILED"
assert_eq "RUN_GENOME_TRACK=1 + reconciled present -> reconciled" \
    "$RECONCILED" "$(eval_helper 1)"

# (ii) genome track ON + reconciled file absent -> legacy path
rm -f "$RECONCILED"
assert_eq "RUN_GENOME_TRACK=1 + reconciled absent -> legacy" \
    "$LEGACY" "$(eval_helper 1)"

# (iii) genome track OFF -> legacy path. A reconciled file is present here (02c
# still writes a byte-identical pass-through copy at toggle=0), so this also
# guards that the toggle — not mere file presence — drives the selection.
: > "$RECONCILED"
assert_eq "RUN_GENOME_TRACK=0 -> legacy (reconciled present ignored)" \
    "$LEGACY" "$(eval_helper 0)"

# (iv) production default: RUN_GENOME_TRACK UNSET + reconciled present ->
# reconciled. Proves the helper's ${RUN_GENOME_TRACK:-1} default enables the
# track (config.sh ships RUN_GENOME_TRACK=1; this pins the same track-on behavior
# for the unset case, e.g. functions.sh sourced without config.sh).
: > "$RECONCILED"
assert_eq "RUN_GENOME_TRACK unset (default) + reconciled present -> reconciled" \
    "$RECONCILED" "$(eval_helper unset)"

echo "OK"
