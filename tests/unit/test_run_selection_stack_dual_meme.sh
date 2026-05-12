#!/bin/bash
# Smoke test for bead -7cy step 4: run_selection_stack.sh must apply
# codon-level ClipKit (strict kpic-smart-gap + lenient smart-gap) on
# MACSE's codon output, run MEME on each, then invoke parse_meme.py
# in dual-mode with --codon. GARD / BUSTED / aBSREL stay on the raw
# MACSE codon (untouched).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
STACK_SCRIPT="$PROJECT_DIR/scripts/hpc/run_selection_stack.sh"

# 1. Syntax
bash -n "$STACK_SCRIPT"

# 2. Strict-trim ClipKit invocation (codon-aware -co, kpic-smart-gap)
grep -qE 'clipkit .*\$CODON_ALIGN.*-co.*-m kpic-smart-gap.*\$CODON_STRICT' "$STACK_SCRIPT" \
    || { echo "FAIL: strict codon-trim ClipKit invocation missing" >&2; exit 1; }

# 3. Lenient-trim ClipKit invocation (codon-aware -co, smart-gap)
grep -qE 'clipkit .*\$CODON_ALIGN.*-co.*-m smart-gap.*\$CODON_LENIENT' "$STACK_SCRIPT" \
    || { echo "FAIL: lenient codon-trim ClipKit invocation missing" >&2; exit 1; }

# 4. Both MEME calls present (multi-line --alignment args; check via count + alignment refs)
N_MEME=$(grep -c '^hyphy meme' "$STACK_SCRIPT" || true)
[ "$N_MEME" -ge 2 ] \
    || { echo "FAIL: expected ≥2 'hyphy meme' invocations, got $N_MEME" >&2; exit 1; }
grep -qE -- '--alignment "\$CODON_STRICT"' "$STACK_SCRIPT" \
    || { echo "FAIL: strict MEME --alignment \$CODON_STRICT not found" >&2; exit 1; }
grep -qE -- '--alignment "\$CODON_LENIENT"' "$STACK_SCRIPT" \
    || { echo "FAIL: lenient MEME --alignment \$CODON_LENIENT not found" >&2; exit 1; }

# 5. parse_meme dual-mode call with --codon and concordance outputs
grep -qE 'parse_meme\.py' "$STACK_SCRIPT" \
    || { echo "FAIL: parse_meme.py not invoked" >&2; exit 1; }
grep -qE -- '--lenient-json' "$STACK_SCRIPT" \
    || { echo "FAIL: parse_meme.py not called with --lenient-json" >&2; exit 1; }
grep -qE -- '--codon\b' "$STACK_SCRIPT" \
    || { echo "FAIL: parse_meme.py not called with --codon flag" >&2; exit 1; }
grep -qE -- '--out-og-concordance' "$STACK_SCRIPT" \
    || { echo "FAIL: parse_meme.py not called with --out-og-concordance" >&2; exit 1; }

# 6. GARD / BUSTED / aBSREL still consume raw $CODON_ALIGN (not trimmed).
#    Verify by checking each tool appears at the top-of-line + the only
#    --alignment "$CODON_ALIGN" references come BEFORE the new MEME block.
grep -q '^hyphy gard' "$STACK_SCRIPT" \
    || { echo "FAIL: GARD invocation missing" >&2; exit 1; }
grep -q '^hyphy absrel' "$STACK_SCRIPT" \
    || { echo "FAIL: aBSREL invocation missing" >&2; exit 1; }
N_BUSTED=$(grep -c '^hyphy busted' "$STACK_SCRIPT" || true)
[ "$N_BUSTED" -eq 2 ] \
    || { echo "FAIL: expected exactly 2 'hyphy busted' invocations, got $N_BUSTED" >&2; exit 1; }
# Confirm $CODON_ALIGN is still referenced as an alignment arg (untrimmed
# inputs preserved for the column-noise-robust tests).
grep -qE -- '--alignment "\$CODON_ALIGN"' "$STACK_SCRIPT" \
    || { echo "FAIL: \$CODON_ALIGN no longer used for GARD/BUSTED/aBSREL" >&2; exit 1; }

# 7. Cumulative concat loop includes meme_concordance
grep -qE 'for variant in busted_s busted_mh meme meme_concordance' "$STACK_SCRIPT" \
    || { echo "FAIL: cumulative concat loop missing meme_concordance variant" >&2; exit 1; }

echo "OK: run_selection_stack.sh dual-MEME smoke test passed"
