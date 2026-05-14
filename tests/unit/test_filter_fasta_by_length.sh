#!/bin/bash
# filter_fasta_by_length: drop protein sequences longer than MAX_AA_LENGTH
# from a FASTA before feeding it to TMbed. Bead -m1f follow-up: prevents the
# quadratic-memory slow-tail that timed out stage 02 job 57653730 at 99.88%
# (transcript-assembly outliers thousands of aa long take 130+ sec/seq on a
# 2080 Ti; real chemoreceptor GPCRs are 300-500 aa).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

TMPDIR=$(mktemp -d /tmp/test_filter_fasta_XXXX)
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

export RESULTS_DIR="$TMPDIR"
export LOGS_DIR="$TMPDIR/logs"
mkdir -p "$LOGS_DIR"

# shellcheck disable=SC1091
source "$PROJECT_DIR/functions.sh"

# Test 1: a mixed FASTA — short (50 aa), medium (1500 aa), oversize (3000 aa).
# With max_len=1500, expect short + medium kept, oversize dropped.
python3 - <<'PYEOF' > "$TMPDIR/mixed.aa"
print(">short_50aa")
print("M" * 50)
print(">medium_1500aa")
print("M" * 1500)
print(">oversize_3000aa")
print("M" * 3000)
PYEOF

filter_fasta_by_length "$TMPDIR/mixed.aa" "$TMPDIR/filtered.aa" 1500 \
    2>"$TMPDIR/filter.log"

# Verify exactly two sequences kept (short + medium), oversize dropped.
KEPT_HEADERS=$(grep -c '^>' "$TMPDIR/filtered.aa" || true)
if [[ "$KEPT_HEADERS" != "2" ]]; then
    echo "FAIL: expected 2 sequences after filter, got $KEPT_HEADERS" >&2
    echo "filtered.aa contents:" >&2
    cat "$TMPDIR/filtered.aa" >&2
    exit 1
fi

if ! grep -q '^>short_50aa$' "$TMPDIR/filtered.aa"; then
    echo "FAIL: short_50aa missing from filtered output" >&2
    exit 1
fi
if ! grep -q '^>medium_1500aa$' "$TMPDIR/filtered.aa"; then
    echo "FAIL: medium_1500aa missing from filtered output (boundary case: <=, not <)" >&2
    exit 1
fi
if grep -q '^>oversize_3000aa$' "$TMPDIR/filtered.aa"; then
    echo "FAIL: oversize_3000aa was kept; should have been dropped" >&2
    exit 1
fi

# stderr should report kept/total count for provenance
if ! grep -q 'kept 2.*3' "$TMPDIR/filter.log"; then
    echo "FAIL: expected 'kept 2/3' provenance line on stderr, got:" >&2
    cat "$TMPDIR/filter.log" >&2
    exit 1
fi

echo "OK: mixed-length filter (kept 2/3, boundary inclusive)"

# Test 2: empty input -> empty output, no crash.
: > "$TMPDIR/empty.aa"
filter_fasta_by_length "$TMPDIR/empty.aa" "$TMPDIR/empty_out.aa" 1500 \
    2>"$TMPDIR/empty.log"
if [[ -s "$TMPDIR/empty_out.aa" ]]; then
    echo "FAIL: empty input produced non-empty output" >&2
    exit 1
fi
echo "OK: empty input produces empty output"

# Test 3: defaults from MAX_AA_LENGTH env var when third arg omitted.
export MAX_AA_LENGTH=100
python3 - <<'PYEOF' > "$TMPDIR/env_test.aa"
print(">a_80")
print("M" * 80)
print(">b_150")
print("M" * 150)
PYEOF
filter_fasta_by_length "$TMPDIR/env_test.aa" "$TMPDIR/env_out.aa" \
    2>"$TMPDIR/env.log"
unset MAX_AA_LENGTH
if [[ "$(grep -c '^>' "$TMPDIR/env_out.aa")" != "1" ]]; then
    echo "FAIL: MAX_AA_LENGTH=100 should keep only a_80; got:" >&2
    cat "$TMPDIR/env_out.aa" >&2
    exit 1
fi
if ! grep -q '^>a_80$' "$TMPDIR/env_out.aa"; then
    echo "FAIL: a_80 missing under MAX_AA_LENGTH=100" >&2
    exit 1
fi
echo "OK: MAX_AA_LENGTH env var honored"

# Test 4: missing input file -> nonzero exit, no output file written.
if filter_fasta_by_length "$TMPDIR/does_not_exist.aa" "$TMPDIR/no_input_out.aa" 1500 \
       2>"$TMPDIR/missing.log"; then
    echo "FAIL: filter_fasta_by_length should fail when input doesn't exist" >&2
    exit 1
fi
echo "OK: missing input raises error"

echo "PASS test_filter_fasta_by_length.sh"
