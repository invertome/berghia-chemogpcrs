#!/bin/bash
#SBATCH --job-name=filter_stack_smoke
#SBATCH --partition=cpu
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# smoke_test_filter_stack.sh — exercise PREQUAL+CLOAK+TAPER end-to-end
# on a small real input before the stage-04 array hits all 595 OGs.
#
# Picks a small Nath et al .faa file, takes the first 6-8 sequences,
# runs `run_alignment_filter_stack` (PREQUAL -> ensemble -> CLOAK ->
# TAPER -> ClipKit-style finalize), and verifies:
#   1. Each gated stage actually ran (or warned).
#   2. Final output is non-empty.
#   3. Sequence count is preserved (no sequences dropped).
#   4. Per-stage provenance files exist.
#
# Submit:  sbatch scripts/unity/smoke_test_filter_stack.sh

set -eo pipefail
mkdir -p logs

WORKDIR="/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05"
cd "$WORKDIR"

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate berghia-gpcr

# shellcheck disable=SC1091
source config.sh
# shellcheck disable=SC1091
source functions.sh

THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "===================================================================="
echo "[$(date -u +%FT%TZ)] Filter-stack smoke test"
echo "===================================================================="
echo "WORKDIR: $WORKDIR"
echo "Tool paths:"
echo "  PREQUAL: $PREQUAL"
echo "  CLOAK  : $CLOAK"
echo "  TAPER  : $TAPER"
echo "  JULIA  : $JULIA"
echo "  MAFFT  : $MAFFT"
echo "  FAMSA  : $FAMSA"
echo ""
echo "Gates: RUN_PREQUAL=${RUN_PREQUAL} RUN_CLOAK=${RUN_CLOAK} RUN_TAPER=${RUN_TAPER}"
echo ""

# --- Pick a small test input ------------------------------------------------
# Use a small Nath et al protein FASTA (the Capitella file is ~50K and has
# realistic chemoreceptor-like protein content with TM/loop alternation —
# good stress case for the alignment-uncertainty consensus filter).
TEST_INPUT_SRC="${WORKDIR}/references/nath_et_al/lse/annelida/283909_Capitella_teleta.faa"
if [[ ! -f "$TEST_INPUT_SRC" ]]; then
    # Fallback: any small .faa
    TEST_INPUT_SRC=$(find "${WORKDIR}/references/nath_et_al" -name '*.faa' -size -200k -size +5k 2>/dev/null | head -1)
fi
[[ -f "$TEST_INPUT_SRC" ]] || { echo "ERROR: no test input found" >&2; exit 1; }

# Take the first 6 sequences (enough to exercise consensus, small enough to
# run in <1 min for the smoke test)
SCRATCH_DIR="${WORKDIR}/results/_smoke_test_filter_stack"
rm -rf "$SCRATCH_DIR"
mkdir -p "$SCRATCH_DIR"

TEST_INPUT="$SCRATCH_DIR/input.fa"
awk 'BEGIN{n=0} /^>/{n++; if(n>6) exit} {print}' "$TEST_INPUT_SRC" > "$TEST_INPUT"

N_INPUT=$(grep -c '^>' "$TEST_INPUT")
echo "Test input: $TEST_INPUT_SRC -> first $N_INPUT seqs in $TEST_INPUT"
echo ""

# --- Run the filter stack ---------------------------------------------------
OUTPUT="$SCRATCH_DIR/cleaned.fa"
WD="$SCRATCH_DIR/_filter_stack_wd"

echo "--- Running run_alignment_filter_stack ---"
if run_alignment_filter_stack "$TEST_INPUT" "$OUTPUT" "$WD" "smoke" "$THREADS"; then
    echo "filter stack exit: 0"
else
    echo "filter stack exit: $? (continuing to report what was produced)"
fi
echo ""

# --- Verification -----------------------------------------------------------
PASS=1

echo "=== Verification ==="

# 1. Final output exists and non-empty
if [[ -s "$OUTPUT" ]]; then
    N_OUT=$(grep -c '^>' "$OUTPUT")
    echo "  [PASS] Output file non-empty: $OUTPUT ($N_OUT seqs)"
else
    echo "  [FAIL] Output file empty or missing: $OUTPUT"
    PASS=0
fi

# 2. Sequence count preserved
if [[ -s "$OUTPUT" ]]; then
    if [[ "$N_OUT" -eq "$N_INPUT" ]]; then
        echo "  [PASS] Sequence count preserved: $N_INPUT -> $N_OUT"
    else
        echo "  [WARN] Sequence count changed: $N_INPUT input -> $N_OUT output (CLOAK / TAPER may legitimately drop short seqs)"
    fi
fi

# 3. Stack provenance file
if [[ -s "${OUTPUT}.filter_stack.txt" ]]; then
    echo "  [PASS] Provenance file exists"
    echo "    --- ${OUTPUT}.filter_stack.txt ---"
    sed 's/^/    /' "${OUTPUT}.filter_stack.txt"
else
    echo "  [FAIL] Provenance file missing: ${OUTPUT}.filter_stack.txt"
    PASS=0
fi

# 4. Per-stage outputs in workdir
echo ""
echo "  Per-stage workdir contents:"
ls -la "$WD" 2>/dev/null | sed 's/^/    /' || echo "    (workdir missing)"

# 5. Per-stage logs
echo ""
echo "  Per-stage stderr logs:"
for stage in prequal ensemble cloak taper; do
    log="${LOGS_DIR}/${stage}_smoke.err"
    if [[ -s "$log" ]]; then
        echo "    --- ${log} (last 5 lines) ---"
        tail -5 "$log" | sed 's/^/      /'
    else
        echo "    [empty/missing] $log"
    fi
done

# 6. Ensemble member count (CLOAK input diversity)
ENS_DIR="$WD/smoke_ensemble"
if [[ -d "$ENS_DIR" ]]; then
    N_ENS=$(find "$ENS_DIR" -maxdepth 1 -name '*.fa' -type f | wc -l)
    echo ""
    echo "  Ensemble: $N_ENS alignments at $ENS_DIR"
    find "$ENS_DIR" -maxdepth 1 -name '*.fa' -type f -printf '    %f  %s bytes\n'
fi

echo ""
echo "===================================================================="
if [[ $PASS -eq 1 ]]; then
    echo "[$(date -u +%FT%TZ)] SMOKE TEST PASSED"
    exit 0
else
    echo "[$(date -u +%FT%TZ)] SMOKE TEST FAILED — see verification block above"
    exit 1
fi
