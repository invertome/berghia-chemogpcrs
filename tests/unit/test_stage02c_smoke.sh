#!/bin/bash
# Smoke test for 02c_genome_reconcile.sh — the RUN_GENOME_TRACK=0 pass-through.
#
# The genome track (RUN_GENOME_TRACK=1) can only be validated on Unity (it needs
# the RefSeq proteins/GFF, the EvidentialGene .mrna, and gmap/minimap2/miniprot/
# blastp). This test pins the OTHER contract: with the toggle OFF, stage 02c must
# reproduce legacy behavior exactly — copy the stage-02 transcriptome candidate
# set through to reconciled_candidates.faa UNCHANGED, invoking NO heavy binaries
# (the pass-through branch exits before any detector/aligner runs).
#
# Hermetic strategy: config.sh derives BASE_DIR from its OWN location
# (realpath(dirname(BASH_SOURCE))), so RESULTS_DIR can't be redirected by an env
# var. We therefore run the stage inside a temp "repo root" — copy config.sh +
# functions.sh + the stage into it — so BASE_DIR (hence RESULTS_DIR/LOGS_DIR)
# points at the throwaway dir and nothing touches the real results/.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
STAGE="$PROJECT_DIR/02c_genome_reconcile.sh"

# 1. Syntax
bash -n "$STAGE"
echo "OK: 02c_genome_reconcile.sh parses"

# 2. Hermetic pass-through run
TMP="$(mktemp -d "${TMPDIR:-/tmp}/test_stage02c_XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
cp "$PROJECT_DIR/config.sh" "$PROJECT_DIR/functions.sh" "$STAGE" "$TMP/"

# Stage-02 transcriptome candidate FASTA fixture (original EvidentialGene ids +
# aalen headers — arbitrary here, the pass-through copies bytes verbatim).
mkdir -p "$TMP/results/chemogpcrs"
CAND="$TMP/results/chemogpcrs/chemogpcrs_berghia.fa"
cat > "$CAND" <<'FA'
>BersteEVm000001t1 type=protein; aalen=350,95%,complete; clen=1053;
MASTALPHAONECHEMORECEPTORSEQUENCERESIDUESHEREMOREMOREMOREMOREMORE
>BersteEVm000002t3 type=protein; aalen=104,97%,partial3; clen=321;
MBETATWOCHEMORECEPTORSEQUENCERESIDUESHEREMOREMOREMOREMOREMOREMORE
FA

# Run only the toggle-OFF branch. No heavy binaries are reachable on this path.
if ! ( cd "$TMP" && RUN_GENOME_TRACK=0 bash 02c_genome_reconcile.sh ) \
        > "$TMP/run.log" 2>&1; then
    echo "FAIL: stage exited non-zero on RUN_GENOME_TRACK=0" >&2
    cat "$TMP/run.log" >&2
    exit 1
fi

OUT="$TMP/results/reconciliation/reconciled_candidates.faa"
if [ ! -s "$OUT" ]; then
    echo "FAIL: pass-through file not written: $OUT" >&2
    cat "$TMP/run.log" >&2
    exit 1
fi
echo "OK: pass-through wrote reconciled_candidates.faa"

# The pass-through must be byte-identical to the transcriptome candidate set.
if ! diff -q "$CAND" "$OUT" >/dev/null; then
    echo "FAIL: pass-through differs from the transcriptome candidate set" >&2
    diff "$CAND" "$OUT" >&2 || true
    exit 1
fi
echo "OK: pass-through is byte-identical to the transcriptome candidate set"

# Completion marker written (downstream/idempotency contract).
if [ ! -f "$TMP/results/step_completed_02c.txt" ]; then
    echo "FAIL: step_completed_02c.txt not written" >&2
    exit 1
fi
echo "OK: completion marker written"

echo "OK"
