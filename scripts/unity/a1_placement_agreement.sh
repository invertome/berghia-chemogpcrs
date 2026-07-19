#!/bin/bash
#SBATCH --job-name=a1_placement_agreement
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# FULL-SCALE test of the amiu premise: does placing candidates on a fixed anchor
# backbone reproduce the distance-to-nearest-reference you get from a DE-NOVO
# joint tree? That confound feeds A1, which residualizes on rank(tree_distance),
# so RANK agreement (Spearman) is the operative statistic, not Pearson.
#
# FULL SCALE AND ML ONLY, DELIBERATELY. An earlier subset run (150 anchors /
# 60 candidates, FastTree) returned Spearman 0.735 — ambiguous against the
# pre-stated >=0.8 pass / <0.6 fail thresholds — but that design confounded the
# question with an inadequate method: FastTree is approximate, and a sparse
# backbone under-represents the reference structure a placement needs. This
# script uses the production sets and IQ-TREE-inferred trees on BOTH sides so the
# comparison measures placement, not the shortcut.
#
# Inputs (both produced by earlier jobs; self-gated):
#   DE-NOVO side : a1_tree_final.sh  -> class_A_a1.treefile   (790 cand + 953 anchors)
#   BACKBONE side: build_a1_backbone.sh -> a1_backbone.treefile
#                                       + backbone_reference_alignment.fa
#                                       + backbone_anchor_ids.txt
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "$REPO"
source config.sh
source functions.sh

A1_DIR="${A1_DIR:-${RESULTS_DIR}/phylogenies/protein/class_A_a1}"
BACKBONE_DIR="${BACKBONE_DIR:-${RESULTS_DIR}/phylogenies/protein/a1_backbone}"
OUT="${OUT:-${RESULTS_DIR}/phylogenies/protein/a1_agreement}"
CAND_FASTA="${CAND_FASTA:-${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia_classA.fa}"
DENOVO_TREE="${DENOVO_TREE:-${A1_DIR}/class_A_a1.treefile}"
REF_ALN="${BACKBONE_DIR}/backbone_reference_alignment.fa"
BACKBONE_TREE="${BACKBONE_DIR}/a1_backbone.treefile"
ANCHOR_IDS="${BACKBONE_DIR}/backbone_anchor_ids.txt"
THREADS="${THREADS:-32}"
mkdir -p "$OUT"

for f in "$DENOVO_TREE" "$REF_ALN" "$BACKBONE_TREE" "$ANCHOR_IDS" "$CAND_FASTA"; do
    [ -s "$f" ] || { echo "[agreement] ERROR: required input missing/empty: $f" >&2; exit 1; }
done
echo "[agreement] candidates=$(grep -c '^>' "$CAND_FASTA")  anchors=$(wc -l < "$ANCHOR_IDS")"

# --- DE-NOVO side: distances straight off the full joint ML tree ---
python3 "${SCRIPTS_DIR}/tree_distance_to_refs.py" \
    --tree "$DENOVO_TREE" --ref-ids "$ANCHOR_IDS" --out "${OUT}/denovo.tsv"

# --- PLACEMENT side: add candidates onto the FIXED backbone columns, then EPA-ng ---
# --keeplength preserves the reference columns so the alignment still matches the
# backbone tree (EPA-ng requires that correspondence).
mafft --add "$CAND_FASTA" --keeplength --thread "$THREADS" "$REF_ALN" \
    > "${OUT}/combined_aln.fa" 2> "${OUT}/mafft_add.log"

python3 - "$CAND_FASTA" "${OUT}/combined_aln.fa" "${OUT}/query_aln.fa" <<'PY'
import sys
cand_fa, combined, out = sys.argv[1:4]
ids = {l[1:].split()[0] for l in open(cand_fa) if l.startswith('>')}
keep = False
with open(out, 'w') as o:
    for line in open(combined):
        if line.startswith('>'):
            keep = line[1:].split()[0] in ids
        if keep:
            o.write(line)
PY

epa-ng --ref-msa "$REF_ALN" --tree "$BACKBONE_TREE" --query "${OUT}/query_aln.fa" \
    --model LG --outdir "$OUT" --redo > "${OUT}/epa.log" 2>&1

python3 "${SCRIPTS_DIR}/placement_distance.py" \
    --jplace "${OUT}/epa_result.jplace" --ref-ids "$ANCHOR_IDS" \
    --out "${OUT}/placement.tsv"

# --- Agreement (Spearman is the operative statistic: A1 uses rank) ---
python3 - "${OUT}/denovo.tsv" "${OUT}/placement.tsv" "${OUT}/agreement.json" <<'PY'
import csv, json, sys
from scipy.stats import spearmanr, pearsonr
def load(p): return {r['id']: float(r['tree_distance']) for r in csv.DictReader(open(p), delimiter='\t')}
d, pl, outp = load(sys.argv[1]), load(sys.argv[2]), sys.argv[3]
ids = sorted(set(d) & set(pl))
res = {"n_denovo": len(d), "n_placement": len(pl), "n_compared": len(ids)}
print(f"compared {len(ids)} candidates (denovo={len(d)}, placement={len(pl)})")
if len(ids) >= 5:
    x = [d[i] for i in ids]; y = [pl[i] for i in ids]
    rs, _ = spearmanr(x, y); rp, _ = pearsonr(x, y)
    res.update(spearman=round(float(rs), 4), pearson=round(float(rp), 4),
               denovo_median=round(sorted(x)[len(x)//2], 4),
               placement_median=round(sorted(y)[len(y)//2], 4))
    print(f"Spearman(de-novo, placement) = {rs:+.3f}   <-- operative (A1 uses rank)")
    print(f"Pearson (de-novo, placement) = {rp:+.3f}")
    print("VERDICT vs pre-stated thresholds: >=0.8 => placement is a sound")
    print("generalization; <0.6 => it does not reproduce the confound, revisit amiu.")
else:
    print("ERROR: too few shared ids — check id normalization on both sides", file=sys.stderr)
json.dump(res, open(outp, 'w'), indent=2)
PY
echo "[agreement] verdict -> ${OUT}/agreement.json"
