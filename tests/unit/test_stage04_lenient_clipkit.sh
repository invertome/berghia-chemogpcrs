#!/bin/bash
# Smoke test for bead -7cy stage 04 change: per-OG path must call ClipKit
# twice — once with -m kpic-smart-gap (strict, → _trimmed.fa) and once with
# -m smart-gap (lenient, → _trimmed_lenient.fa). Both outputs are run
# through drop_near_all_gap_rows.py. The global + per-LSE paths must NOT
# be touched (MEME only runs per-OG).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
STAGE_SCRIPT="$PROJECT_DIR/04_phylogenetic_analysis.sh"

# 1. Syntax check — catches typos / unclosed quotes from the edit
bash -n "$STAGE_SCRIPT"

# 2. The strict kpic-smart-gap call exists (existing behavior)
grep -qE '\$\{base\}_clipkit\b.*kpic-smart-gap.*\$\{base\}_trimmed\.fa' "$STAGE_SCRIPT" \
    || { echo "FAIL: strict kpic-smart-gap call missing from per-OG path" >&2; exit 1; }

# 3. The new lenient smart-gap call exists with the right flag + output path
grep -qE '\$\{base\}_clipkit_lenient\b.*-m smart-gap.*\$\{base\}_trimmed_lenient\.fa' "$STAGE_SCRIPT" \
    || { echo "FAIL: lenient smart-gap call missing or mis-flagged" >&2; exit 1; }

# 4. drop_near_all_gap_rows runs on the lenient output too
grep -qE 'drop_near_all_gap_rows\.py' "$STAGE_SCRIPT" | head -1 >/dev/null
grep -qE '\$\{base\}_trimmed_lenient\.fa' "$STAGE_SCRIPT" \
    || { echo "FAIL: lenient trimmed alignment path not present" >&2; exit 1; }
grep -qE 'drop_gappy_lenient' "$STAGE_SCRIPT" \
    || { echo "FAIL: drop_near_all_gap_rows.py log path for lenient missing" >&2; exit 1; }

# 5. Global + per-LSE paths must NOT have lenient ClipKit calls (MEME is
#    per-OG only; touching the others would be over-trimming work).
if grep -nE "(all_berghia_refs_clipkit_lenient|lse_.*_clipkit_lenient)" "$STAGE_SCRIPT"; then
    echo "FAIL: lenient ClipKit accidentally leaked into global or per-LSE paths" >&2
    exit 1
fi

echo "OK: stage 04 lenient ClipKit smoke test passed"
