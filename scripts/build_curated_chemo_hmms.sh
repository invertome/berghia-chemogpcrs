#!/bin/bash
# build_curated_chemo_hmms.sh — Align + hmmbuild per curated invertebrate
# chemoreceptor family.
#
# Reads per-family FASTAs at:
#   references/curated_chemoreceptors_lophotrochozoa/<family>/<family>_seqs.fa
#
# Writes per-family HMMs at:
#   ${RESULTS_DIR}/classification/hmms_curated_chemo/<family>/<family>.hmm
#
# Alignment: MAFFT-DASH E-INS-i (--genafpair --dash). E-INS-i is the
# right MAFFT strategy for GPCRs (sequences with shared conserved core
# domains separated by variable loops — 7 TM helices + ECL/ICL).
# --dash adds anchors from PDB structural homologs.
# Falls back to plain E-INS-i (no DASH) on DASH-server failure.
#
# Skips any family with <3 sequences (hmmbuild needs >=3 for a useful
# profile; <3 produces a degenerate HMM that flags everything).
#
# Bead -k0g, 2026-05-19.

set -eo pipefail

REFS_DIR="${REFS_DIR:-references/curated_chemoreceptors_lophotrochozoa}"
OUT_DIR="${OUT_DIR:-${RESULTS_DIR:-results}/classification/hmms_curated_chemo}"
CPUS="${CPUS:-${SLURM_CPUS_PER_TASK:-4}}"

mkdir -p "$OUT_DIR"

shopt -s nullglob
family_dirs=("$REFS_DIR"/*/)
shopt -u nullglob

if [ "${#family_dirs[@]}" -eq 0 ]; then
    echo "ERROR: no family directories in $REFS_DIR" >&2
    exit 1
fi

built=0
skipped_insufficient=0
total=0

for family_dir in "${family_dirs[@]}"; do
    family=$(basename "${family_dir%/}")

    # Multiple FASTA sources per family are supported via *_seqs*.fa glob.
    # Pattern intent: ${family}_seqs.fa (NCBI fetcher output) +
    # ${family}_seqs_<source>.fa (per-source contributions, e.g.,
    # audino2024). Combined into a single input for MAFFT.
    shopt -s nullglob
    seq_files=("$family_dir"/*_seqs*.fa)
    shopt -u nullglob
    if [ "${#seq_files[@]}" -eq 0 ]; then
        echo "[$family] no *_seqs*.fa file in $family_dir; skipping" >&2
        continue
    fi

    out_family_dir="$OUT_DIR/$family"
    mkdir -p "$out_family_dir"
    seqs="$out_family_dir/${family}_combined_seqs.fa"
    cat "${seq_files[@]}" > "$seqs"

    n_seqs=$(grep -c '^>' "$seqs" 2>/dev/null || true)
    n_seqs=${n_seqs:-0}
    total=$((total + 1))
    echo "[$family] combined $n_seqs sequences from ${#seq_files[@]} source(s): $(printf '%s ' "${seq_files[@]##*/}")" >&2

    if [ "$n_seqs" -lt 3 ]; then
        echo "[$family] only $n_seqs sequences (need >=3 for hmmbuild); skipping" >&2
        skipped_insufficient=$((skipped_insufficient + 1))
        continue
    fi

    aligned="$out_family_dir/${family}_aligned.fa"
    hmm="$out_family_dir/${family}.hmm"

    echo "[$family] aligning $n_seqs sequences via MAFFT-DASH E-INS-i" >&2

    # MAFFT-DASH first; fall back to plain E-INS-i on DASH failure.
    # --genafpair = E-INS-i strategy (designed for sequences with shared
    # core domains separated by variable regions; matches GPCR
    # architecture).
    # --originalseqonly: DASH structural homologs guide the alignment but are
    # excluded from output, so the HMM models only the curated references.
    if ! mafft --dash --originalseqonly --maxiterate 1000 --genafpair --thread "$CPUS" \
            "$seqs" > "$aligned" 2>"$out_family_dir/${family}_mafft.log"; then
        echo "[$family] DASH failed; falling back to plain E-INS-i" >&2
        mafft --maxiterate 1000 --genafpair --thread "$CPUS" \
            "$seqs" > "$aligned" 2>"$out_family_dir/${family}_mafft.log"
    fi

    if [ ! -s "$aligned" ]; then
        echo "[$family] alignment empty; check ${family}_mafft.log" >&2
        continue
    fi

    n_aligned=$(grep -c '^>' "$aligned")
    echo "[$family] aligned: $n_aligned seqs; running hmmbuild" >&2

    hmmbuild --cpu "$CPUS" --amino \
        -n "$family" \
        "$hmm" "$aligned" \
        > "$out_family_dir/${family}_hmmbuild.log" 2>&1

    if [ ! -s "$hmm" ]; then
        echo "[$family] hmmbuild produced empty HMM" >&2
        continue
    fi

    echo "[$family] HMM built: $hmm" >&2
    built=$((built + 1))
done

echo "" >&2
echo "Summary: built=$built  total_families=$total  skipped_insufficient=$skipped_insufficient" >&2

# Consolidate per-family HMMs into a single library file for easier
# hmmsearch invocation by stage 02. find+xargs is ARG_MAX-safe (we
# learned this lesson in bead -m1f, stage 01 conserved.hmm bug).
combined="$OUT_DIR/curated_chemo.hmm"
find "$OUT_DIR" -maxdepth 2 -name '*.hmm' -not -name 'curated_chemo.hmm' \
    -print0 | sort -z | xargs -0 cat > "$combined"
if [ -s "$combined" ]; then
    n_hmms=$(grep -c '^NAME' "$combined")
    echo "Consolidated $n_hmms HMMs at $combined" >&2
fi

[ "$built" -gt 0 ] || {
    echo "ERROR: no HMMs built. Check that manifest has >=3 entries per family." >&2
    exit 2
}
