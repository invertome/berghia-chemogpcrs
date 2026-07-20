#!/bin/bash
# run_selection_stack.sh — Modern HyPhy selection stack per orthogroup.
#
# Bead -urk. Replaces the aBSREL-only path in 05_selective_pressure_and_asr.sh
# with the recommended multi-method pipeline:
#
#   GARD (recombination screen) → BUSTED-S (gene-wide w/ synonymous-rate
#   variation) → BUSTED-MH (multi-nucleotide-substitution correction;
#   Lucaci 2023) → aBSREL (branch-level on confirmed signal) → MEME
#   (site-level for confirmed branches).
#
# Branch testing: aBSREL is run EXHAUSTIVELY over all branches. A targeted
# foreground pre-filter (Berghia + LSE-internal branches only) would have more
# power, but it requires the input tree to carry {FOREGROUND} labels and no
# component of this pipeline produces such a tree. Passing the optional 4th
# positional argument (BRANCH_LIST_FILE) is therefore rejected outright rather
# than silently degrading to an exhaustive run — see bead -0mhv(b). Multiple-
# testing correction over aBSREL branch p-values is applied downstream by
# rank_candidates.py (BH-FDR via statsmodels), not here.
#
# This is HPC-only — runtime is ~hours per OG depending on size.

#SBATCH --job-name=selection_stack
#SBATCH --output=logs/selection_stack_%A_%a.out
#SBATCH --error=logs/selection_stack_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$PROJECT_DIR/config.sh"
source "$PROJECT_DIR/functions.sh"

# Inputs (positional)
CODON_ALIGN="${1:?codon alignment FASTA}"
TREE="${2:?tree file}"
OG_BASE="${3:?orthogroup base name}"
BRANCH_LIST_FILE="${4:-}"   # Reserved; see the branch-testing note in the header.

OUT_DIR="${RESULTS_DIR}/selective_pressure"
mkdir -p "$OUT_DIR"

# Bead -ih5u: this script runs under `#SBATCH --array=0-999%50` in stage 05, so
# up to 50 copies write into the same shared OUT_DIR concurrently. Every path
# that is NOT namespaced by ${OG_BASE} needs a per-task-unique staging name.
# $$ alone is insufficient: OUT_DIR is on shared storage and PIDs are per-node,
# so two array tasks on different nodes can hold the same PID. Compose the
# SLURM array identity with the PID so the tag is unique cluster-wide.
TASK_TAG="${SLURM_JOB_ID:-nojob}.${SLURM_ARRAY_TASK_ID:-0}.$$"

# --- Per-method status tracking (bead -0mhv part a) --------------------------
# Previously every HyPhy call was suffixed `|| log --level=WARN`, so all five
# methods could fail and the script still exited 0 and logged "complete". The
# caller's own `|| log WARN` guard at 05:396 could therefore never fire and the
# dN/dS axis silently contributed nothing. We now record an explicit per-method
# status, publish it as a machine-readable per-OG TSV, and exit non-zero when
# no inference method produced usable output.
#
# "Success" requires BOTH exit 0 AND a non-empty output artifact: HyPhy can
# exit 0 having written nothing, and a zero-length JSON would sail through the
# `[ -f ... ]` parse guards below and yield an empty parse.
METHOD_NAMES=()
METHOD_STATUS=()

record_status() {
    METHOD_NAMES+=("$1")
    METHOD_STATUS+=("$2")
}

# check_method <name> <expected-output> <exit-code>
# Returns 0 on success so callers can branch on it; never aborts under `set -e`.
check_method() {
    local name="$1" outfile="$2" rc="$3"
    if [ "$rc" -ne 0 ]; then
        log --level=WARN "${name} failed for ${OG_BASE} (exit ${rc})"
        record_status "$name" "failed"
        return 1
    fi
    if [ ! -s "$outfile" ]; then
        log --level=WARN "${name} exited 0 for ${OG_BASE} but wrote no output to ${outfile}"
        record_status "$name" "no_output"
        return 1
    fi
    record_status "$name" "ok"
    return 0
}

# Reject the unimplemented branch pre-filter loudly rather than accepting an
# argument that cannot be honoured (bead -0mhv part b).
if [ -n "$BRANCH_LIST_FILE" ]; then
    log --level=ERROR "BRANCH_LIST_FILE was supplied ('${BRANCH_LIST_FILE}') but foreground branch pre-filtering is NOT implemented: aBSREL's --branches FOREGROUND requires a tree pre-labelled with {FOREGROUND} tags and nothing in this pipeline produces one. Refusing to run rather than silently testing all branches under a targeted-test label."
    exit 2
fi

# 1. GARD — recombination breakpoint screen
log "GARD: ${OG_BASE}"
rc=0
hyphy gard \
    --alignment "$CODON_ALIGN" \
    --output "$OUT_DIR/${OG_BASE}_gard.json" \
    > "$OUT_DIR/${OG_BASE}_gard.log" 2>&1 || rc=$?
check_method "gard" "$OUT_DIR/${OG_BASE}_gard.json" "$rc" \
    || log --level=WARN "continuing without recombination correction for ${OG_BASE}"

# 2. BUSTED-S — synonymous-rate variation
log "BUSTED-S: ${OG_BASE}"
rc=0
hyphy busted \
    --alignment "$CODON_ALIGN" \
    --tree "$TREE" \
    --srv Yes \
    --output "$OUT_DIR/${OG_BASE}_busted_s.json" \
    > "$OUT_DIR/${OG_BASE}_busted_s.log" 2>&1 || rc=$?
check_method "busted_s" "$OUT_DIR/${OG_BASE}_busted_s.json" "$rc" || true

# 3. BUSTED-MH — multi-nucleotide-substitution correction
log "BUSTED-MH: ${OG_BASE}"
rc=0
hyphy busted \
    --alignment "$CODON_ALIGN" \
    --tree "$TREE" \
    --srv Yes \
    --multiple-hits Double+Triple \
    --output "$OUT_DIR/${OG_BASE}_busted_mh.json" \
    > "$OUT_DIR/${OG_BASE}_busted_mh.log" 2>&1 || rc=$?
check_method "busted_mh" "$OUT_DIR/${OG_BASE}_busted_mh.json" "$rc" || true

# 4. aBSREL — branch-level, exhaustive over all branches (see header).
log "aBSREL: ${OG_BASE}"
rc=0
hyphy absrel \
    --alignment "$CODON_ALIGN" \
    --tree "$TREE" \
    --output "$OUT_DIR/${OG_BASE}_absrel.json" \
    > "$OUT_DIR/${OG_BASE}_absrel.log" 2>&1 || rc=$?
check_method "absrel" "$OUT_DIR/${OG_BASE}_absrel.json" "$rc" || true

# 5. MEME — site-level on the OG (HyPhy emits per-site posteriors).
#    Bead -7cy: MEME runs on a ClipKit-trimmed codon alignment, not on the
#    raw MACSE output. Two passes:
#      - strict  (kpic-smart-gap): drops singleton-variant + gappy columns
#      - lenient (smart-gap only): drops only gappy columns
#    parse_meme.py concordance-stratifies the per-site results; high-confidence
#    = MEME+ under BOTH passes feeds rank_candidates' positive_score boost.
#    GARD, BUSTED-S, BUSTED-MH, aBSREL stay on the raw MACSE output above —
#    they're column-noise robust and benefit from the wider alignment.
#
#    Bead -0mhv part c: a ClipKit failure MUST NOT fall back to copying the raw
#    alignment into both slots. That made the two "passes" the same file, so the
#    concordance check compared an alignment against itself, scored trivially
#    100% concordant, and fed a fabricated high-confidence signal into
#    rank_candidates' positive_score boost. It is also structurally invalid:
#    parse_meme.map_strict_to_lenient documents "strict is a column-subset of
#    lenient", and the raw untrimmed alignment is a SUPERSET of both. A failed
#    trim now leaves its output absent, which is reported and skips that pass.
CODON_STRICT="$OUT_DIR/${OG_BASE}_codon_strict.fa"
CODON_LENIENT="$OUT_DIR/${OG_BASE}_codon_lenient.fa"
log "ClipKit codon-aware trims: ${OG_BASE} (strict + lenient)"
rc=0
clipkit "$CODON_ALIGN" -co -m kpic-smart-gap -o "$CODON_STRICT" \
    > "$OUT_DIR/${OG_BASE}_clipkit_codon_strict.log" 2>&1 || rc=$?
if [ "$rc" -ne 0 ] || [ ! -s "$CODON_STRICT" ]; then
    log --level=ERROR "ClipKit strict-trim failed for ${OG_BASE} (exit ${rc}); MEME strict pass and the strict/lenient concordance check will be SKIPPED (no raw-alignment fallback — see bead -0mhv)"
    rm -f "$CODON_STRICT"          # drop any truncated partial output
    record_status "clipkit_strict" "failed"
else
    record_status "clipkit_strict" "ok"
fi
rc=0
clipkit "$CODON_ALIGN" -co -m smart-gap -o "$CODON_LENIENT" \
    > "$OUT_DIR/${OG_BASE}_clipkit_codon_lenient.log" 2>&1 || rc=$?
if [ "$rc" -ne 0 ] || [ ! -s "$CODON_LENIENT" ]; then
    log --level=ERROR "ClipKit lenient-trim failed for ${OG_BASE} (exit ${rc}); MEME lenient pass and the strict/lenient concordance check will be SKIPPED (no raw-alignment fallback — see bead -0mhv)"
    rm -f "$CODON_LENIENT"
    record_status "clipkit_lenient" "failed"
else
    record_status "clipkit_lenient" "ok"
fi

# NOTE ON INDENTATION: the `hyphy meme` invocations below sit at column 0 inside
# their guards on purpose. tests/unit/test_run_selection_stack_dual_meme.sh
# asserts the invocations are line-anchored (`grep -c '^hyphy meme'`), so they
# must not be indented into their enclosing `if`. Same for the GARD/BUSTED/
# aBSREL calls above, which are already top-level.

log "MEME (strict): ${OG_BASE}"
rc=0
if [ -f "$CODON_STRICT" ]; then
hyphy meme \
    --alignment "$CODON_STRICT" \
    --tree "$TREE" \
    --output "$OUT_DIR/${OG_BASE}_meme.json" \
    > "$OUT_DIR/${OG_BASE}_meme.log" 2>&1 || rc=$?
else
    log --level=WARN "MEME (strict) skipped for ${OG_BASE}: no strict-trimmed alignment"
    rc=1
fi
check_method "meme_strict" "$OUT_DIR/${OG_BASE}_meme.json" "$rc" || true

log "MEME (lenient): ${OG_BASE}"
rc=0
if [ -f "$CODON_LENIENT" ]; then
hyphy meme \
    --alignment "$CODON_LENIENT" \
    --tree "$TREE" \
    --output "$OUT_DIR/${OG_BASE}_meme_lenient.json" \
    > "$OUT_DIR/${OG_BASE}_meme_lenient.log" 2>&1 || rc=$?
else
    log --level=WARN "MEME (lenient) skipped for ${OG_BASE}: no lenient-trimmed alignment"
    rc=1
fi
check_method "meme_lenient" "$OUT_DIR/${OG_BASE}_meme_lenient.json" "$rc" || true

# Auto-parse BUSTED-S, BUSTED-MH, MEME JSONs to per-OG CSVs, appended to
# the cumulative results files that rank_candidates.py reads. Each parser
# is idempotent (atomic per-OG; concat happens here).
# Bead -ih5u: the parser stderr logs are per-OG, not shared. The previous
# shared `parse_meme.log` / `parse_busted_*.log` were appended to by all 50
# concurrent array tasks, interleaving (and, above PIPE_BUF, tearing) the
# diagnostics from different orthogroups into an unattributable stream.
SCRIPTS_DIR="${SCRIPTS_DIR:-${PROJECT_DIR}/scripts}"

if [ -s "$OUT_DIR/${OG_BASE}_busted_s.json" ]; then
    python3 "$SCRIPTS_DIR/parse_busted.py" \
        --json "$OUT_DIR/${OG_BASE}_busted_s.json" \
        --og-name "$OG_BASE" --variant S \
        --out "$OUT_DIR/${OG_BASE}_busted_s.csv" \
        2>>"$OUT_DIR/${OG_BASE}_parse_busted_s.log" || true
fi
if [ -s "$OUT_DIR/${OG_BASE}_busted_mh.json" ]; then
    python3 "$SCRIPTS_DIR/parse_busted.py" \
        --json "$OUT_DIR/${OG_BASE}_busted_mh.json" \
        --og-name "$OG_BASE" --variant MH \
        --out "$OUT_DIR/${OG_BASE}_busted_mh.csv" \
        2>>"$OUT_DIR/${OG_BASE}_parse_busted_mh.log" || true
fi
if [ -s "$OUT_DIR/${OG_BASE}_meme.json" ]; then
    # Bead -7cy: dual-mode when the lenient pass succeeded. Falls back to
    # strict-only parse when lenient is absent (preserves existing behavior).
    #
    # Bead -0mhv part c: dual-mode additionally requires the two trimmed
    # alignments to differ. Identical inputs make the concordance metric
    # degenerate — every strict-positive site maps onto itself, so
    # high_confidence == all positives and alignment_robustness_index == 1.0
    # by construction, which is not evidence of alignment robustness. That can
    # arise legitimately (both trim regimes removing the same columns), so it
    # is a skip-with-warning rather than an error, but it must never be scored.
    MEME_DUAL_OK=0
    if [ -s "$OUT_DIR/${OG_BASE}_meme_lenient.json" ] \
       && [ -s "$CODON_STRICT" ] && [ -s "$CODON_LENIENT" ]; then
        if cmp -s "$CODON_STRICT" "$CODON_LENIENT"; then
            log --level=WARN "strict and lenient trims are byte-identical for ${OG_BASE}; skipping the concordance pass (a self-comparison is trivially 100% concordant and would fabricate a high-confidence signal)"
            record_status "meme_concordance" "degenerate_identical_trims"
        else
            MEME_DUAL_OK=1
        fi
    else
        record_status "meme_concordance" "skipped_missing_inputs"
    fi

    if [ "$MEME_DUAL_OK" = "1" ]; then
        python3 "$SCRIPTS_DIR/parse_meme.py" \
            --json "$OUT_DIR/${OG_BASE}_meme.json" \
            --lenient-json "$OUT_DIR/${OG_BASE}_meme_lenient.json" \
            --strict-fa "$CODON_STRICT" \
            --lenient-fa "$CODON_LENIENT" \
            --codon \
            --og-name "$OG_BASE" \
            --out-og "$OUT_DIR/${OG_BASE}_meme.csv" \
            --out-sites "$OUT_DIR/${OG_BASE}_meme_sites.csv" \
            --out-og-concordance "$OUT_DIR/${OG_BASE}_meme_concordance.csv" \
            --out-sites-concordance "$OUT_DIR/${OG_BASE}_meme_concordance_sites.csv" \
            2>>"$OUT_DIR/${OG_BASE}_parse_meme.log" || true
        [ -s "$OUT_DIR/${OG_BASE}_meme_concordance.csv" ] \
            && record_status "meme_concordance" "ok" \
            || record_status "meme_concordance" "no_output"
    else
        # Stale concordance output from an earlier run would otherwise be
        # re-concatenated into the cumulative CSV by the loop below.
        rm -f "$OUT_DIR/${OG_BASE}_meme_concordance.csv" \
              "$OUT_DIR/${OG_BASE}_meme_concordance_sites.csv"
        python3 "$SCRIPTS_DIR/parse_meme.py" \
            --json "$OUT_DIR/${OG_BASE}_meme.json" \
            --og-name "$OG_BASE" \
            --out-og "$OUT_DIR/${OG_BASE}_meme.csv" \
            --out-sites "$OUT_DIR/${OG_BASE}_meme_sites.csv" \
            2>>"$OUT_DIR/${OG_BASE}_parse_meme.log" || true
    fi
fi

# Concatenate per-OG CSVs into the cumulative results consumed by
# rank_candidates.py (fast, atomic; safe to redo each invocation).
# Bead -7cy: added meme_concordance to the variants — same concat pattern.
#
# Bead -ih5u: the staging file MUST be per-task-unique. It was a fixed
# "${cum}.tmp" shared by all 50 concurrent array tasks: task A's `: >` truncated
# the file while task B was mid-append, and B then published the truncated
# result over the cumulative CSV. This mirrors the fix already applied to the
# aBSREL path in 05_selective_pressure_and_asr.sh:257 (`${cumulative}.tmp.$$`,
# bead 1z7) which was never propagated to the four variants below. The final
# `mv` is atomic within the directory, so readers never observe a partial file
# and the last task to finish publishes the complete cumulative.
for variant in busted_s busted_mh meme meme_concordance; do
    cum="$OUT_DIR/${variant}_results.csv"
    # meme_concordance has a different output filename suffix
    if [ "$variant" = "meme_concordance" ]; then
        cum="$OUT_DIR/meme_concordance.csv"
    fi
    tmp="${cum}.tmp.${TASK_TAG}"
    : > "$tmp"
    first=1
    for csv in "$OUT_DIR"/*_${variant}.csv; do
        [ -f "$csv" ] || continue
        # Skip per-site variants (different schema)
        case "$csv" in
            *_meme_sites.csv) continue ;;
            *_meme_concordance_sites.csv) continue ;;
        esac
        if [ "$first" = "1" ]; then
            cat "$csv" >> "$tmp"
            first=0
        else
            tail -n +2 "$csv" >> "$tmp"
        fi
    done
    if [ -s "$tmp" ]; then
        mv "$tmp" "$cum"
    else
        rm -f "$tmp"
    fi
done

# --- Publish the per-method status + decide the exit code -------------------
# The status TSV is per-OG (no cross-task contention) and is what a caller
# should read to learn which methods contributed for this orthogroup.
STATUS_TSV="$OUT_DIR/${OG_BASE}_selection_status.tsv"
status_tmp="${STATUS_TSV}.tmp.${TASK_TAG}"
{
    printf 'og_name\tmethod\tstatus\n'
    for i in "${!METHOD_NAMES[@]}"; do
        printf '%s\t%s\t%s\n' "$OG_BASE" "${METHOD_NAMES[$i]}" "${METHOD_STATUS[$i]}"
    done
} > "$status_tmp"
mv "$status_tmp" "$STATUS_TSV"

# Exit policy: GARD is advisory (a recombination screen; nothing downstream
# consumes gard.json) and the ClipKit/concordance rows are diagnostics, so
# neither gates the exit code. The stack has produced usable evidence iff at
# least one INFERENCE method succeeded. Zero successes means this orthogroup
# contributes nothing to the dN/dS axis, which the caller must be able to see —
# 05_selective_pressure_and_asr.sh:396 already has the `|| log WARN` guard that
# this exit code finally makes reachable.
n_ok=0
for i in "${!METHOD_NAMES[@]}"; do
    case "${METHOD_NAMES[$i]}" in
        busted_s|busted_mh|absrel|meme_strict)
            [ "${METHOD_STATUS[$i]}" = "ok" ] && n_ok=$((n_ok + 1))
            ;;
    esac
done

summary=""
for i in "${!METHOD_NAMES[@]}"; do
    summary+="${METHOD_NAMES[$i]}=${METHOD_STATUS[$i]} "
done

if [ "$n_ok" -eq 0 ]; then
    log --level=ERROR "Selection stack produced NO usable inference output for ${OG_BASE} [${summary%% }] — this orthogroup contributes nothing to the dN/dS axis"
    exit 1
fi

log "Selection stack complete for ${OG_BASE}: ${n_ok}/4 inference methods succeeded [${summary%% }]"
