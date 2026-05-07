#!/bin/bash
#SBATCH --job-name=filter_stack_large
#SBATCH --partition=cpu
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# smoke_test_filter_stack_large.sh — exercise the filter stack at
# realistic scale before stage 04 fires the 595-OG array.
#
# The small smoke (smoke_test_filter_stack.sh, 6 seqs) confirmed wiring;
# this one tests scale: 150 sequences from a real Nath et al .faa
# (Conus consors, gastropoda — same broad family as our chemoreceptor
# orthogroups). Catches issues that only show up at scale:
#   - MAFFT variant runtime (LINSI/GINSI/EINSI on hundreds of seqs)
#   - CLOAK consensus column-count handling
#   - TAPER on a many-seq alignment
#   - Per-stage memory profiles
#
# Submit:  sbatch scripts/unity/smoke_test_filter_stack_large.sh

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

THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "===================================================================="
echo "[$(date -u +%FT%TZ)] Filter-stack LARGE smoke test (~150 seqs)"
echo "===================================================================="
echo "Threads: $THREADS"
echo "Gates: RUN_PREQUAL=${RUN_PREQUAL} RUN_CLOAK=${RUN_CLOAK} RUN_TAPER=${RUN_TAPER}"
echo ""

# Use a real Nath protein FASTA with hundreds of sequences. Conus consors
# gastropod proteome is realistic for our use case (chemoreceptor-like
# divergence within a single genus expansion).
TEST_INPUT_SRC="${WORKDIR}/references/nath_et_al/lse/gastropoda/101297_Conus_consors.faa"
[[ -f "$TEST_INPUT_SRC" ]] || { echo "ERROR: test input missing: $TEST_INPUT_SRC" >&2; exit 1; }

SCRATCH_DIR="${WORKDIR}/results/_smoke_test_filter_stack_large"
rm -rf "$SCRATCH_DIR"
mkdir -p "$SCRATCH_DIR"

TEST_INPUT="$SCRATCH_DIR/input.fa"
awk 'BEGIN{n=0} /^>/{n++; if(n>150) exit} {print}' "$TEST_INPUT_SRC" > "$TEST_INPUT"
N_INPUT=$(grep -c '^>' "$TEST_INPUT")
INPUT_BYTES=$(stat -c '%s' "$TEST_INPUT")
echo "Test input: $TEST_INPUT_SRC -> first $N_INPUT seqs in $TEST_INPUT (${INPUT_BYTES} bytes)"
echo ""

# --- Run the filter stack with timing ---
OUTPUT="$SCRATCH_DIR/cleaned.fa"
WD="$SCRATCH_DIR/_filter_stack_wd"

echo "--- Running run_alignment_filter_stack (timed) ---"
START_TS=$(date +%s)
if run_alignment_filter_stack "$TEST_INPUT" "$OUTPUT" "$WD" "smoke_large" "$THREADS"; then
    RC=0
else
    RC=$?
fi
END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo ""
echo "filter stack exit: $RC  elapsed: ${ELAPSED}s"
echo ""

# --- Verification at scale ---
PASS=1
echo "=== Verification ==="

# 1. Output non-empty
if [[ -s "$OUTPUT" ]]; then
    N_OUT=$(grep -c '^>' "$OUTPUT")
    OUT_BYTES=$(stat -c '%s' "$OUTPUT")
    echo "  [PASS] Output non-empty: $N_OUT seqs, ${OUT_BYTES} bytes"
else
    echo "  [FAIL] Output empty/missing"
    PASS=0
fi

# 2. Sequence count preserved (the filter stack must not drop sequences)
if [[ -s "$OUTPUT" ]] && [[ "$N_OUT" -eq "$N_INPUT" ]]; then
    echo "  [PASS] Sequence count preserved: $N_INPUT -> $N_OUT"
elif [[ -s "$OUTPUT" ]]; then
    echo "  [WARN] Sequence count changed: $N_INPUT -> $N_OUT"
fi

# 3. Canonical sibling preserved (for 04b ECL analysis)
CANONICAL_SIBLING="${OUTPUT%.fa}_canonical.fa"
if [[ -s "$CANONICAL_SIBLING" ]]; then
    N_CANON=$(grep -c '^>' "$CANONICAL_SIBLING")
    echo "  [PASS] Canonical sibling: $N_CANON seqs ($(stat -c '%s' "$CANONICAL_SIBLING") bytes)"
else
    echo "  [FAIL] Canonical sibling missing"
    PASS=0
fi

# 4. Per-stage timing + sizes (for performance budgeting)
echo ""
echo "  Per-stage outputs (size, n_seqs):"
for stage_file in \
    "${WD}/smoke_large_prequal.fa:PREQUAL" \
    "${WD}/smoke_large_ensemble/canonical.fa:CANONICAL_ALN" \
    "${WD}/smoke_large_cloaked.fa:CLOAK" \
    "${WD}/smoke_large_tapered.fa:TAPER"; do
    path="${stage_file%%:*}"
    label="${stage_file##*:}"
    if [[ -s "$path" ]]; then
        size=$(stat -c '%s' "$path")
        n=$(grep -c '^>' "$path" 2>/dev/null || echo 0)
        printf "    %-15s %10d bytes %5d seqs\n" "$label" "$size" "$n"
    else
        printf "    %-15s [missing]\n" "$label"
    fi
done

# 5. Ensemble member count + sizes
ENS_DIR="$WD/smoke_large_ensemble"
echo ""
echo "  Ensemble members:"
if [[ -d "$ENS_DIR" ]]; then
    find "$ENS_DIR" -maxdepth 1 -name '*.fa' -type f -printf '    %-30f %10s bytes\n' | sort
fi

# 6. CLOAK column-mask quality: does the cloaked alignment have FEWER columns
#    than the canonical (proving CLOAK actually masked something)?
if [[ -s "${WD}/smoke_large_ensemble/canonical.fa" ]] && [[ -s "${WD}/smoke_large_cloaked.fa" ]]; then
    canonical_cols=$(awk '/^>/{if(s){print length(s); exit} s=""} !/^>/{s=s $0}' "${WD}/smoke_large_ensemble/canonical.fa")
    cloaked_cols=$(awk '/^>/{if(s){print length(s); exit} s=""} !/^>/{s=s $0}' "${WD}/smoke_large_cloaked.fa")
    echo ""
    echo "  Column counts: canonical=${canonical_cols}, cloaked=${cloaked_cols}"
    if [[ "$cloaked_cols" -lt "$canonical_cols" ]]; then
        masked=$((canonical_cols - cloaked_cols))
        echo "  [INFO] CLOAK masked $masked uncertain columns (~$(awk -v m=$masked -v c=$canonical_cols 'BEGIN{printf "%.1f%%", m*100/c}') of canonical alignment)"
    elif [[ "$cloaked_cols" -eq "$canonical_cols" ]]; then
        echo "  [WARN] CLOAK output has same column count as canonical — alignments fully agreed (unusual at this scale)"
    elif [[ "$cloaked_cols" -gt "$canonical_cols" ]]; then
        excess=$((cloaked_cols - canonical_cols))
        echo "  [INFO] CLOAK super-alignment +${excess} cols ($(awk -v e=$excess -v c=$canonical_cols 'BEGIN{printf "%.1f%%", e*100/c}') over canonical) — disputed pairings spread into singleton columns"
    fi
fi

# 7. End-to-end with downstream ClipKit (kpic-smart-gap, the stage-04 default).
#    Validates that ClipKit cleanly trims CLOAK's super-alignment back down
#    to a tractable column count for IQ-TREE.
echo ""
echo "  ClipKit kpic-smart-gap (stage-04 trim step):"
CLIPKIT_OUT="${SCRATCH_DIR}/clipkit_kpic_smart_gap.fa"
CLIPKIT_LOG="${SCRATCH_DIR}/clipkit.log"
if clipkit "$OUTPUT" -m kpic-smart-gap -o "$CLIPKIT_OUT" >"$CLIPKIT_LOG" 2>&1 && [[ -s "$CLIPKIT_OUT" ]]; then
    trimmed_cols=$(awk '/^>/{if(s){print length(s); exit} s=""} !/^>/{s=s $0}' "$CLIPKIT_OUT")
    n_trimmed=$(grep -c '^>' "$CLIPKIT_OUT")
    if [[ -n "$cloaked_cols" ]] && [[ "$cloaked_cols" -gt 0 ]]; then
        kept_pct=$(awk -v t=$trimmed_cols -v c=$cloaked_cols 'BEGIN{printf "%.1f%%", t*100/c}')
    else
        kept_pct="?"
    fi
    if [[ -n "$canonical_cols" ]] && [[ "$canonical_cols" -gt 0 ]]; then
        vs_canonical=$(awk -v t=$trimmed_cols -v c=$canonical_cols 'BEGIN{printf "%.0f%%", t*100/c}')
    else
        vs_canonical="?"
    fi
    echo "  [PASS] ClipKit trimmed: $trimmed_cols cols ($n_trimmed seqs)"
    echo "         (kept $kept_pct of CLOAK super-alignment; $vs_canonical of canonical $canonical_cols cols)"
    if [[ "$trimmed_cols" -ge 200 ]] && [[ "$trimmed_cols" -le $((canonical_cols * 2)) ]]; then
        echo "  [PASS] Trimmed column count is sensible (>=200, <=2x canonical) — IQ-TREE-ready"
    elif [[ "$trimmed_cols" -lt 200 ]]; then
        echo "  [WARN] Trimmed alignment has only $trimmed_cols cols — may be too short for ModelFinder"
    elif [[ "$trimmed_cols" -gt $((canonical_cols * 2)) ]]; then
        echo "  [WARN] Trimmed alignment is >2x canonical ($trimmed_cols vs $canonical_cols) — under-trimmed?"
    fi
else
    echo "  [FAIL] ClipKit kpic-smart-gap failed (see $CLIPKIT_LOG)"
    PASS=0
fi

echo ""
echo "===================================================================="
if [[ $PASS -eq 1 ]]; then
    echo "[$(date -u +%FT%TZ)] LARGE SMOKE TEST PASSED (${ELAPSED}s wall)"
    exit 0
else
    echo "[$(date -u +%FT%TZ)] LARGE SMOKE TEST FAILED — see above"
    exit 1
fi
