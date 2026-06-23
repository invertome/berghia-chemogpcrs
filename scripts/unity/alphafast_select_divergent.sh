#!/bin/bash
# alphafast_select_divergent.sh — pick the prototype's 'divergent' and 'typical'
# Berghia correctness candidates by best-hit %identity to curated chemoreceptor
# references, among COMPLETE + Code candidates (real receptors). Replaces the
# bad Noncode-codepot proxy (bead q0o.1).
#
#   divergent = lowest best-hit %id to any known chemoreceptor (furthest from known)
#   typical   = median best-hit %id
#
#SBATCH --job-name=af_divsel
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=20:00
#SBATCH --output=logs/af_divsel-%j.out
#SBATCH --error=logs/af_divsel-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate base                       # mmseqs
set -u

WS="/scratch3/workspace/$USER-jorge"
R="$WS/chemogpcrs_2026-05"
CANDS="$R/results/chemogpcrs/chemogpcrs_berghia.fa"
OUT="$WS/alphafast_proto/divsel"; mkdir -p "$OUT" "$OUT/tmp"

# 1. high-quality candidate subset: fully complete ORF + Code coding potential
awk '/^>/{keep=($0 ~ /complete;/ && $0 ~ /codepot=Code/)} keep' "$CANDS" > "$OUT/cand_hq.fasta"
echo "complete+Code candidates: $(grep -c '^>' "$OUT/cand_hq.fasta")"

# 2. curated chemoreceptor reference set (recursive) + anchor set
: > "$OUT/refs.fasta"
find "$R/references/curated_chemoreceptors_lophotrochozoa" \( -name '*.fa' -o -name '*.fasta' \) -exec cat {} + >> "$OUT/refs.fasta" 2>/dev/null || true
[ -f "$R/references/anchors/anchor_set.fasta" ] && cat "$R/references/anchors/anchor_set.fasta" >> "$OUT/refs.fasta"
echo "reference seqs: $(grep -c '^>' "$OUT/refs.fasta")"

# 3. sensitive search candidates -> references
mmseqs easy-search "$OUT/cand_hq.fasta" "$OUT/refs.fasta" "$OUT/hits.m8" "$OUT/tmp" \
    --format-output "query,target,pident,alnlen,evalue,bits" -s 7.5 --max-seqs 25 -e 1e-3 >/dev/null 2>&1

# 4. best-hit %id per candidate; ascending -> lowest = most divergent
sort -k1,1 -k3,3nr "$OUT/hits.m8" | awk -F'\t' '!s[$1]++{print $1"\t"$3}' | sort -k2,2n > "$OUT/best_pident_sorted.tsv"
N=$(wc -l < "$OUT/best_pident_sorted.tsv")
echo "candidates with a reference hit: $N"
echo "=== 5 most divergent (lowest best-hit %id) ==="; head -5 "$OUT/best_pident_sorted.tsv"
echo "=== median (typical) ==="; sed -n "$(( (N+1)/2 ))p" "$OUT/best_pident_sorted.tsv"
DIV=$(head -1 "$OUT/best_pident_sorted.tsv" | cut -f1)
TYP=$(sed -n "$(( (N+1)/2 ))p" "$OUT/best_pident_sorted.tsv" | cut -f1)
{ echo "divergent_id=$DIV"; echo "typical_id=$TYP"; } | tee "$OUT/PICKS.txt"
