#!/bin/bash
# Integration smoke test for identify_gpcr_candidates() in functions.sh.
#
# Bead -cpq: this test exists because we hit caller-callee flag-name
# mismatches between functions.sh and scripts/identify_gpcrs.py three
# times in May 2026 — each time it took a stage 02 SLURM run to surface
# the failure. The class of bugs is "functions.sh passes --some-flag
# but identify_gpcrs.py's argparse doesn't define that flag (or just
# renamed it)". Pure pytest of merge_gpcr_evidence() missed it because
# the tests call the Python kwarg directly, never via CLI.
#
# What this test does (no real compute, runs anywhere):
#   1. bash -n syntax check on functions.sh
#   2. Static grep: identify_gpcrs.py invocation in functions.sh passes
#      ALL the expected --flags for the current architecture
#   3. Argparse-vs-caller cross-check: every --flag the bash caller
#      passes is ALSO defined in identify_gpcrs.py's argparse
#   4. classify_via_hmm.py and hmmsearch invocations use the right
#      argument forms (catches Unity-only env-activation breakage)
#   5. Stub-and-execute: mock hmmsearch + classify_via_hmm.py via PATH,
#      provide canned tblout/TSV outputs, call identify_gpcr_candidates
#      end-to-end, verify expected output FASTA + ID file get written.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
FUNCTIONS="$PROJECT_DIR/functions.sh"
IDENTIFY_PY="$PROJECT_DIR/scripts/identify_gpcrs.py"

# 1. Syntax
bash -n "$FUNCTIONS"
echo "OK: functions.sh parses"

# 2. Expected flag set in the identify_gpcrs.py invocation in functions.sh
for flag in '--classification-tsv' '--tiammat-tblout' '--curated-chemo-tblout' \
            '--ids-out'; do
    grep -qE -- "$flag " "$FUNCTIONS" \
        || { echo "FAIL: functions.sh does not pass $flag to identify_gpcrs.py" >&2; exit 1; }
done
echo "OK: functions.sh passes expected --flags"

# 3. Argparse-vs-caller cross-check: extract flags from identify_gpcrs.py's
#    argparse and the flags the bash caller passes; verify the caller's
#    flag set is a SUBSET of the argparse-defined flags.
python3 - "$FUNCTIONS" "$IDENTIFY_PY" <<'PYEOF'
import re, sys
functions_sh, identify_py = sys.argv[1], sys.argv[2]

# argparse-defined flags
with open(identify_py) as f:
    py_src = f.read()
defined = set(re.findall(r'ap\.add_argument\(\s*"(--[a-z][a-z0-9-]*)"', py_src))

# Caller's invocation block: find the `python3 ... identify_gpcrs.py ...`
# call in functions.sh and pull every `--xxx-yyy` token within it.
with open(functions_sh) as f:
    sh_src = f.read()
# Match from 'identify_gpcrs.py' to the closing `||` or `}` of the if-block
m = re.search(r'identify_gpcrs\.py(.*?)(?=\}\s*\n|\|\|\s*\{)', sh_src,
              re.DOTALL)
if not m:
    print("could not locate identify_gpcrs.py invocation block in functions.sh",
          file=sys.stderr)
    sys.exit(1)
passed = set(re.findall(r'(--[a-z][a-z0-9-]*)', m.group(1)))

# Whitelist non-identify_gpcrs flags that may appear in the block (none today)
WHITELIST = set()
extras = passed - defined - WHITELIST
if extras:
    print(f"FAIL: functions.sh passes flags NOT defined in identify_gpcrs.py "
          f"argparse: {sorted(extras)}", file=sys.stderr)
    print(f"  caller passes: {sorted(passed)}", file=sys.stderr)
    print(f"  argparse defines: {sorted(defined)}", file=sys.stderr)
    sys.exit(1)
print(f"OK: argparse cross-check ({len(passed)} caller flags, "
      f"all subset of {len(defined)} argparse-defined flags)")
PYEOF

# 4. hmmsearch + classify_via_hmm.py invocations: verify the bash function
#    uses the expected forms (catches accidental flag drift).
grep -qE 'classify_via_hmm\.py' "$FUNCTIONS" \
    || { echo "FAIL: classify_via_hmm.py not invoked in functions.sh" >&2; exit 1; }
grep -qE '\bhmmsearch\b' "$FUNCTIONS" \
    || { echo "FAIL: hmmsearch not invoked in functions.sh" >&2; exit 1; }
grep -qE -- '--tblout' "$FUNCTIONS" \
    || { echo "FAIL: --tblout not used in functions.sh" >&2; exit 1; }
echo "OK: classify_via_hmm.py + hmmsearch invocations present"

# 5. Stub-and-execute end-to-end smoke test.
TMPDIR=$(mktemp -d /tmp/test_identify_gpcr_candidates_XXXX)
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

# Minimal env setup
export RESULTS_DIR="$TMPDIR/results"
export LOGS_DIR="$TMPDIR/logs"
export SCRIPTS_DIR="$TMPDIR/stub_scripts"
export BASE_DIR="$TMPDIR"
export CPUS=2
export SEQTK=seqtk

mkdir -p "$RESULTS_DIR" "$LOGS_DIR" "$SCRIPTS_DIR" "$TMPDIR/bin"
mkdir -p "$RESULTS_DIR/classification/hmms/pfam_fallback"
mkdir -p "$RESULTS_DIR/hmms"
mkdir -p "$RESULTS_DIR/classification/hmms_curated_chemo"

# Tiny input FASTA: 3 sequences with predictable IDs
cat > "$TMPDIR/input.fa" <<'INPUT'
>seq_A
MASTRSEQUENCEALPHAONESOMERESIDUESHEREMOREMOREMOREMOREMOREMOREMOREMOREMORE
>seq_B
MASTRSEQUENCEBRAVOTWOONLYSOMEDIFFERENTRESIDUESHEREMOREMOREMOREMOREMOREMORE
>seq_C
MASTRSEQUENCECHARLIETHREEYETMOREDIFFERENTRESIDUESHEREMOREMOREMOREMOREMOREM
INPUT

mkdir -p "$BASE_DIR/references"
# Stub HMM files (non-empty so the existence-check passes)
echo "STUB HMM" > "$RESULTS_DIR/classification/hmms_curated_chemo/curated_chemo.hmm"
echo "STUB HMM" > "$BASE_DIR/references/tiammat_mollusca_gpcr.hmm"
export TIAMMAT_GPCR_HMM="$BASE_DIR/references/tiammat_mollusca_gpcr.hmm"
export CURATED_CHEMO_HMM="$RESULTS_DIR/classification/hmms_curated_chemo/curated_chemo.hmm"

# Stub hmmsearch: produce a tblout flagging seq_A as a hit
cat > "$TMPDIR/bin/hmmsearch" <<'HMMSEARCH'
#!/bin/bash
# Stub: parse out --tblout target and write a hmmsearch-format tblout
# claiming seq_A as a hit. Args order: --cpu N -E E --tblout PATH HMM_DB SEQ_DB
tblout=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tblout) tblout="$2"; shift 2 ;;
        *) shift ;;
    esac
done
[[ -n "$tblout" ]] || exit 0
mkdir -p "$(dirname "$tblout")"
cat > "$tblout" <<TBL
#                                                               --- full sequence ----
# target name        accession  query name           accession    E-value  score  bias
# ------------------- ---------- -------------------- ---------- --------- ------ -----
seq_A                -          stub_hmm             -           1.0e-50   180.0   0.0
TBL
HMMSEARCH

# Stub classify_via_hmm.py: produce a classification TSV with seq_B classified
cat > "$SCRIPTS_DIR/classify_via_hmm.py" <<'CLASSIFY'
#!/usr/bin/env python3
import sys, os
# Find --output-tsv arg
out = None
for i, a in enumerate(sys.argv):
    if a == "--output-tsv" and i + 1 < len(sys.argv):
        out = sys.argv[i + 1]
if out:
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        f.write("candidate_id\thmm_family\thmm_subfamily\tevalue\tscore\tevidence\n")
        f.write("seq_B\taminergic\t5HT\t1.0e-100\t340.5\thmm\n")
CLASSIFY
chmod +x "$SCRIPTS_DIR/classify_via_hmm.py"

# identify_gpcrs.py: real, not stubbed — it's pure-Python merge logic
# we want to test in vivo. Copy the real script into the stub dir so
# the SCRIPTS_DIR override picks it up.
cp "$PROJECT_DIR/scripts/identify_gpcrs.py" "$SCRIPTS_DIR/identify_gpcrs.py"

# Stub seqtk subseq: just extract by ID using awk
cat > "$TMPDIR/bin/seqtk" <<'SEQTK'
#!/bin/bash
# Stub seqtk: only `subseq INFILE IDFILE` is implemented
if [[ "$1" != "subseq" ]]; then
    echo "stub seqtk: only subseq supported" >&2
    exit 1
fi
awk -v idfile="$3" '
BEGIN { while ((getline line < idfile) > 0) ids[line] = 1 }
/^>/ {
    name = substr($1, 2)
    keep = (name in ids)
    if (keep) print
    next
}
{ if (keep) print }
' "$2"
SEQTK

chmod +x "$TMPDIR/bin/hmmsearch" "$TMPDIR/bin/seqtk"
export PATH="$TMPDIR/bin:$PATH"
export SEQTK="$TMPDIR/bin/seqtk"

# Source functions.sh and call identify_gpcr_candidates
# shellcheck disable=SC1091
source "$FUNCTIONS"

OUT_FA="$TMPDIR/results/output.fa"
CENSUS="$TMPDIR/results/census.tsv"
identify_gpcr_candidates "$TMPDIR/input.fa" "$OUT_FA" "$CENSUS" \
    > "$TMPDIR/identify.log" 2>&1 \
    || { echo "FAIL: identify_gpcr_candidates returned non-zero" >&2;
         tail -20 "$TMPDIR/identify.log" >&2; exit 1; }

# Verify outputs
[[ -s "$OUT_FA" ]] || { echo "FAIL: output FASTA not written" >&2; exit 1; }
n_out=$(grep -c '^>' "$OUT_FA")
# Expected: seq_A (from hmmsearch tiammat/curated stubs) + seq_B (from classify_via_hmm)
# Both should be in the output.
[[ "$n_out" -ge 2 ]] || { echo "FAIL: expected >=2 sequences in output, got $n_out" >&2; cat "$OUT_FA" >&2; exit 1; }
grep -q '^>seq_A$' "$OUT_FA" \
    || { echo "FAIL: seq_A missing from output (hmmsearch tiammat/curated should have caught it)" >&2; exit 1; }
grep -q '^>seq_B$' "$OUT_FA" \
    || { echo "FAIL: seq_B missing from output (classify_via_hmm stub should have classified it)" >&2; exit 1; }
echo "OK: end-to-end stub-execution produced $n_out sequences in output"

# Census TSV exists + has expected schema
[[ -s "$CENSUS" ]] || { echo "FAIL: census TSV not written" >&2; exit 1; }
head -1 "$CENSUS" | grep -qE '^seq_id\sfamily\ssubfamily\sevalue\ssource' \
    || { echo "FAIL: census header malformed: $(head -1 "$CENSUS")" >&2; exit 1; }
echo "OK: census TSV has expected schema"

echo "PASS test_identify_gpcr_candidates.sh"
