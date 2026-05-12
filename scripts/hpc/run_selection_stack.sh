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
# Branch pre-filter: only test Berghia branches + LSE-internal branches
# (the targeted approach has dramatically more power than exhaustive
# correction; aBSREL paper notes this explicitly).
# BH-FDR (q=0.10) applied across the pre-filtered branch set.
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
BRANCH_LIST_FILE="${4:-}"   # Optional: file with one branch name per line for pre-filter

OUT_DIR="${RESULTS_DIR}/selective_pressure"
mkdir -p "$OUT_DIR"

# 1. GARD — recombination breakpoint screen
log "GARD: ${OG_BASE}"
hyphy gard \
    --alignment "$CODON_ALIGN" \
    --output "$OUT_DIR/${OG_BASE}_gard.json" \
    > "$OUT_DIR/${OG_BASE}_gard.log" 2>&1 \
    || log --level=WARN "GARD failed for ${OG_BASE} (continuing without recombination correction)"

# 2. BUSTED-S — synonymous-rate variation
log "BUSTED-S: ${OG_BASE}"
hyphy busted \
    --alignment "$CODON_ALIGN" \
    --tree "$TREE" \
    --srv Yes \
    --output "$OUT_DIR/${OG_BASE}_busted_s.json" \
    > "$OUT_DIR/${OG_BASE}_busted_s.log" 2>&1 \
    || log --level=WARN "BUSTED-S failed for ${OG_BASE}"

# 3. BUSTED-MH — multi-nucleotide-substitution correction
log "BUSTED-MH: ${OG_BASE}"
hyphy busted \
    --alignment "$CODON_ALIGN" \
    --tree "$TREE" \
    --srv Yes \
    --multiple-hits Double+Triple \
    --output "$OUT_DIR/${OG_BASE}_busted_mh.json" \
    > "$OUT_DIR/${OG_BASE}_busted_mh.log" 2>&1 \
    || log --level=WARN "BUSTED-MH failed for ${OG_BASE}"

# 4. aBSREL — branch-level (optionally restricted to a pre-filtered branch set)
log "aBSREL: ${OG_BASE}"
ABSREL_BRANCH_ARG=""
if [ -n "$BRANCH_LIST_FILE" ] && [ -f "$BRANCH_LIST_FILE" ]; then
    ABSREL_BRANCH_ARG="--branches FOREGROUND"
    # Note: HyPhy expects branch labels in the tree; pre-labeling the tree
    # with FOREGROUND tags is the user's responsibility (see scripts/
    # label_tree_branches.py — TODO).
fi
hyphy absrel \
    --alignment "$CODON_ALIGN" \
    --tree "$TREE" \
    --output "$OUT_DIR/${OG_BASE}_absrel.json" \
    $ABSREL_BRANCH_ARG \
    > "$OUT_DIR/${OG_BASE}_absrel.log" 2>&1 \
    || log --level=WARN "aBSREL failed for ${OG_BASE}"

# 5. MEME — site-level on the OG (HyPhy emits per-site posteriors).
#    Bead -7cy: MEME runs on a ClipKit-trimmed codon alignment, not on the
#    raw MACSE output. Two passes:
#      - strict  (kpic-smart-gap): drops singleton-variant + gappy columns
#      - lenient (smart-gap only): drops only gappy columns
#    parse_meme.py concordance-stratifies the per-site results; high-confidence
#    = MEME+ under BOTH passes feeds rank_candidates' positive_score boost.
#    GARD, BUSTED-S, BUSTED-MH, aBSREL stay on the raw MACSE output above —
#    they're column-noise robust and benefit from the wider alignment.
CODON_STRICT="$OUT_DIR/${OG_BASE}_codon_strict.fa"
CODON_LENIENT="$OUT_DIR/${OG_BASE}_codon_lenient.fa"
log "ClipKit codon-aware trims: ${OG_BASE} (strict + lenient)"
clipkit "$CODON_ALIGN" -co -m kpic-smart-gap -o "$CODON_STRICT" \
    > "$OUT_DIR/${OG_BASE}_clipkit_codon_strict.log" 2>&1 \
    || { log --level=WARN "ClipKit strict-trim failed for ${OG_BASE}; falling back to raw MACSE output"; cp "$CODON_ALIGN" "$CODON_STRICT"; }
clipkit "$CODON_ALIGN" -co -m smart-gap -o "$CODON_LENIENT" \
    > "$OUT_DIR/${OG_BASE}_clipkit_codon_lenient.log" 2>&1 \
    || { log --level=WARN "ClipKit lenient-trim failed for ${OG_BASE}; falling back to raw MACSE output"; cp "$CODON_ALIGN" "$CODON_LENIENT"; }

log "MEME (strict): ${OG_BASE}"
hyphy meme \
    --alignment "$CODON_STRICT" \
    --tree "$TREE" \
    --output "$OUT_DIR/${OG_BASE}_meme.json" \
    > "$OUT_DIR/${OG_BASE}_meme.log" 2>&1 \
    || log --level=WARN "MEME (strict) failed for ${OG_BASE}"

log "MEME (lenient): ${OG_BASE}"
hyphy meme \
    --alignment "$CODON_LENIENT" \
    --tree "$TREE" \
    --output "$OUT_DIR/${OG_BASE}_meme_lenient.json" \
    > "$OUT_DIR/${OG_BASE}_meme_lenient.log" 2>&1 \
    || log --level=WARN "MEME (lenient) failed for ${OG_BASE}"

log "Selection stack complete for ${OG_BASE}"

# Auto-parse BUSTED-S, BUSTED-MH, MEME JSONs to per-OG CSVs, appended to
# the cumulative results files that rank_candidates.py reads. Each parser
# is idempotent (atomic per-OG; concat happens here).
SCRIPTS_DIR="${SCRIPTS_DIR:-${PROJECT_DIR}/scripts}"

if [ -f "$OUT_DIR/${OG_BASE}_busted_s.json" ]; then
    python3 "$SCRIPTS_DIR/parse_busted.py" \
        --json "$OUT_DIR/${OG_BASE}_busted_s.json" \
        --og-name "$OG_BASE" --variant S \
        --out "$OUT_DIR/${OG_BASE}_busted_s.csv" \
        2>>"$OUT_DIR/parse_busted_s.log" || true
fi
if [ -f "$OUT_DIR/${OG_BASE}_busted_mh.json" ]; then
    python3 "$SCRIPTS_DIR/parse_busted.py" \
        --json "$OUT_DIR/${OG_BASE}_busted_mh.json" \
        --og-name "$OG_BASE" --variant MH \
        --out "$OUT_DIR/${OG_BASE}_busted_mh.csv" \
        2>>"$OUT_DIR/parse_busted_mh.log" || true
fi
if [ -f "$OUT_DIR/${OG_BASE}_meme.json" ]; then
    # Bead -7cy: dual-mode when the lenient pass succeeded. Falls back to
    # strict-only parse when lenient is absent (preserves existing behavior).
    if [ -f "$OUT_DIR/${OG_BASE}_meme_lenient.json" ] \
       && [ -f "$CODON_STRICT" ] && [ -f "$CODON_LENIENT" ]; then
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
            2>>"$OUT_DIR/parse_meme.log" || true
    else
        python3 "$SCRIPTS_DIR/parse_meme.py" \
            --json "$OUT_DIR/${OG_BASE}_meme.json" \
            --og-name "$OG_BASE" \
            --out-og "$OUT_DIR/${OG_BASE}_meme.csv" \
            --out-sites "$OUT_DIR/${OG_BASE}_meme_sites.csv" \
            2>>"$OUT_DIR/parse_meme.log" || true
    fi
fi

# Concatenate per-OG CSVs into the cumulative results consumed by
# rank_candidates.py (fast, atomic; safe to redo each invocation).
# Bead -7cy: added meme_concordance to the variants — same concat pattern.
for variant in busted_s busted_mh meme meme_concordance; do
    cum="$OUT_DIR/${variant}_results.csv"
    # meme_concordance has a different output filename suffix
    if [ "$variant" = "meme_concordance" ]; then
        cum="$OUT_DIR/meme_concordance.csv"
    fi
    : > "${cum}.tmp"
    first=1
    for csv in "$OUT_DIR"/*_${variant}.csv; do
        [ -f "$csv" ] || continue
        # Skip per-site variants (different schema)
        case "$csv" in
            *_meme_sites.csv) continue ;;
            *_meme_concordance_sites.csv) continue ;;
        esac
        if [ "$first" = "1" ]; then
            cat "$csv" >> "${cum}.tmp"
            first=0
        else
            tail -n +2 "$csv" >> "${cum}.tmp"
        fi
    done
    if [ -s "${cum}.tmp" ]; then
        mv "${cum}.tmp" "$cum"
    else
        rm -f "${cum}.tmp"
    fi
done
