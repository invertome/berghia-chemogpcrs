#!/bin/bash
# Regression: create_checkpoint must work when results/checkpoints/ doesn't
# exist yet. Before the fix, atomic_write's `echo > $tempfile` failed silently
# because the parent dir was missing; the JSON .checkpoint file was never
# written. Resumability fell back on the legacy step_completed_*.txt marker
# (functions.sh:178 OR's the two), so the bug went unnoticed in production
# until stage 01 job 57594440 emitted stderr noise on 2026-05-12.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

TMPDIR=$(mktemp -d /tmp/test_create_checkpoint_XXXX)
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

export RESULTS_DIR="$TMPDIR"
export LOGS_DIR="$TMPDIR/logs"
mkdir -p "$LOGS_DIR"
export PIPELINE_RUN_ID="test_run"

# DELIBERATELY do NOT mkdir "$RESULTS_DIR/checkpoints" — the regression
# is that create_checkpoint must handle the missing-dir case itself.

# shellcheck disable=SC1091
source "$PROJECT_DIR/functions.sh"

# Capture stderr so we can fail if atomic_write spews "No such file or directory"
ERR_LOG="$TMPDIR/create_checkpoint.err"
create_checkpoint "my_test_step" 2>"$ERR_LOG"

if grep -q 'No such file or directory' "$ERR_LOG"; then
    echo "FAIL: create_checkpoint emitted 'No such file or directory' on stderr:" >&2
    cat "$ERR_LOG" >&2
    exit 1
fi

CKPT="$RESULTS_DIR/checkpoints/my_test_step.checkpoint"
if [[ ! -f "$CKPT" ]]; then
    echo "FAIL: checkpoint file not created at $CKPT" >&2
    exit 1
fi

python3 - "$CKPT" <<'PYEOF'
import json, sys
with open(sys.argv[1]) as f:
    data = json.load(f)
assert data["step"] == "my_test_step", f"step name mismatch: {data}"
assert data["exit_code"] == 0, f"exit_code mismatch: {data}"
print("OK: create_checkpoint writes JSON checkpoint even when dir missing")
PYEOF

# Second call must also succeed (mkdir -p is idempotent)
create_checkpoint "my_test_step_2" 2>>"$ERR_LOG"
if [[ ! -f "$RESULTS_DIR/checkpoints/my_test_step_2.checkpoint" ]]; then
    echo "FAIL: second checkpoint file missing" >&2
    exit 1
fi

echo "OK: create_checkpoint regression test passed"
