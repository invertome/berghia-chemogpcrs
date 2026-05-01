#!/bin/bash
# Test that record_provenance actually writes step records to the JSON file.
# Bead -ryr.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

TMPDIR=$(mktemp -d /tmp/test_provenance_XXXX)
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

# Set up minimum environment for functions.sh to source cleanly
export RESULTS_DIR="$TMPDIR"
export LOGS_DIR="$TMPDIR/logs"
mkdir -p "$LOGS_DIR" "$RESULTS_DIR/provenance"
export PIPELINE_RUN_ID="test_run"
export PROVENANCE_FILE="$RESULTS_DIR/provenance/run_${PIPELINE_RUN_ID}.json"

# Initialize a minimal provenance file
cat > "$PROVENANCE_FILE" <<EOF
{
    "run_id": "${PIPELINE_RUN_ID}",
    "started_at": "2026-01-01T00:00:00",
    "steps": []
}
EOF

# Source functions.sh (it should be safe with these env vars set; only the
# init_pipeline function performs file I/O, and we won't call it).
# shellcheck disable=SC1091
source "$PROJECT_DIR/functions.sh"

echo "input data" > "$TMPDIR/input.txt"
echo "output data" > "$TMPDIR/output.txt"

# First step
record_provenance "test_step_1" "echo hello world" "$TMPDIR/input.txt" "$TMPDIR/output.txt"

python3 - "$PROVENANCE_FILE" <<'PYEOF'
import json, sys
with open(sys.argv[1]) as f:
    prov = json.load(f)
assert "steps" in prov, f"missing 'steps' in {prov}"
assert len(prov["steps"]) == 1, f"expected 1 step, got {len(prov['steps'])}"
step = prov["steps"][0]
assert step["name"] == "test_step_1", f"name mismatch: {step['name']}"
assert "input_checksums" in step
assert "output_checksums" in step
print("OK: provenance step 1 recorded")
PYEOF

# Second step — must append, not overwrite
record_provenance "test_step_2" "ls -la" "$TMPDIR/input.txt" "$TMPDIR/output.txt"

python3 - "$PROVENANCE_FILE" <<'PYEOF'
import json, sys
with open(sys.argv[1]) as f:
    prov = json.load(f)
assert len(prov["steps"]) == 2, f"expected 2 steps, got {len(prov['steps'])}"
assert prov["steps"][0]["name"] == "test_step_1"
assert prov["steps"][1]["name"] == "test_step_2"
print("OK: provenance correctly appended a second step")
PYEOF

echo "PASS test_provenance.sh"
