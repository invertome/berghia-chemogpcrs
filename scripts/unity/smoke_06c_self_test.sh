#!/bin/bash
#SBATCH --job-name=06c_smoke_self
#SBATCH --partition=cpu
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# 06c smoke test: run the non-chemoreceptor classifier against the
# curated reference set itself. Each reference must classify back to
# its known family. Validates HMM scan + placement + consensus.
# OG-vote is skipped (no Orthogroups.tsv on a self-test).

set -eo pipefail
cd /scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
mkdir -p logs

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate berghia-gpcr

THREADS="${SLURM_CPUS_PER_TASK:-4}"

REF_FASTA="references/non_chemo_gpcr/all_references.fasta"
REF_TSV="references/non_chemo_gpcr/all_references.tsv"
HMM_DIR="results/classification/hmms"
PFAM_DIR="results/classification/hmms/pfam_fallback"
LOO_METRICS="results/classification/loo/loo_metrics.tsv"
TREE_DIR="results/classification/trees"

OUT_DIR="results/classification/smoke_self_test"
mkdir -p "$OUT_DIR"

echo "[$(date +%T)] === 06c smoke test (self-classify references) ==="
N_REFS=$(grep -c '^>' "$REF_FASTA")
echo "  refs in: $N_REFS"

# 1. HMM scan
echo "[$(date +%T)] HMM scan…"
PFAM_ARGS=""
[ -d "$PFAM_DIR" ] && PFAM_ARGS="--pfam-fallback-dir $PFAM_DIR"
python3 scripts/classify_via_hmm.py \
    --candidate-fasta "$REF_FASTA" \
    --hmm-dir "$HMM_DIR" \
    $PFAM_ARGS \
    --loo-metrics "$LOO_METRICS" \
    --output-tsv "$OUT_DIR/hmm.tsv" \
    --threads "$THREADS"

# 2. Placement
echo "[$(date +%T)] Placement…"
python3 scripts/classify_via_placement.py \
    --candidate-fasta "$REF_FASTA" \
    --tree-dir "$TREE_DIR" \
    --output-tsv "$OUT_DIR/placement.tsv" \
    --threads "$THREADS" \
    || echo "[$(date +%T)] WARN placement failed; consensus will use HMM only"

# 3. Skip OG-vote (no orthogroups in a self-test): synth empty TSV
echo -e "candidate_id\tog_id\tog_vote_family\tog_vote_subfamily\tn_members\tn_annotated\tconsensus_fraction" > "$OUT_DIR/og.tsv"

# 4. Consensus
echo "[$(date +%T)] Consensus…"
PLACE_ARG=""
[ -f "$OUT_DIR/placement.tsv" ] && PLACE_ARG="--placement-tsv $OUT_DIR/placement.tsv"
# shellcheck disable=SC2086
python3 scripts/classify_consensus.py \
    --hmm-tsv "$OUT_DIR/hmm.tsv" \
    --og-tsv "$OUT_DIR/og.tsv" \
    $PLACE_ARG \
    --out "$OUT_DIR/consensus.tsv"

# 5. Compare to truth
echo "[$(date +%T)] Compare consensus to truth…"
python3 - <<'PY'
import csv
from collections import Counter

truth = {}
with open("references/non_chemo_gpcr/all_references.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        truth[r["accession"]] = (r.get("family", ""), r.get("subfamily", ""))

obs = {}
with open("results/classification/smoke_self_test/consensus.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        obs[r["candidate_id"]] = (
            r.get("classification", ""),
            r.get("classification_family", ""),
            r.get("classification_subfamily", ""),
            r.get("classification_confidence", ""),
        )

n_total = 0
n_fam_ok = 0
n_sub_ok = 0
n_subapplicable = 0
err_fam = Counter()
err_class = Counter()
non_chemo_n = 0
for acc, (exp_fam, exp_sub) in truth.items():
    n_total += 1
    if acc not in obs:
        err_class[("MISSING", exp_fam)] += 1
        continue
    cls, fam, sub, conf = obs[acc]
    if cls in ("non-chemoreceptor", "likely-non-chemoreceptor"):
        non_chemo_n += 1
    if fam == exp_fam:
        n_fam_ok += 1
    else:
        err_fam[(exp_fam, fam)] += 1
    if exp_sub:
        n_subapplicable += 1
        if sub == exp_sub:
            n_sub_ok += 1
print(f"total refs: {n_total}")
print(f"called as non-chemo (TP): {non_chemo_n}/{n_total} "
      f"= {non_chemo_n/n_total:.1%}")
print(f"family-correct: {n_fam_ok}/{n_total} = {n_fam_ok/n_total:.1%}")
if n_subapplicable:
    print(f"subfamily-correct (where applicable): "
          f"{n_sub_ok}/{n_subapplicable} = {n_sub_ok/n_subapplicable:.1%}")
if err_fam:
    print("\nTop family confusions (expected -> observed):")
    for (e, o), n in err_fam.most_common(15):
        print(f"  {e:<25} -> {o:<25} ({n})")
PY
echo "[$(date +%T)] === smoke test done ==="
