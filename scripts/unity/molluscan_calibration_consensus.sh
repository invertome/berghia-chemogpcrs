#!/bin/bash
#SBATCH --job-name=moll_cal_consensus
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Cross-model consensus of the molluscan-calibrated exclusion set.
#
# The production novelty channel is a CONSENSUS of proteinclip3b + protrek, so
# the decision-relevant output is not either model's exclusion list but their
# agreement. Reports the vertebrate-envelope and molluscan-calibrated sets per
# model, their intersections, and whether the two models agree on WHICH family
# each consensus exclusion belongs to -- agreeing that a candidate is "inside
# something" while disagreeing on what is a much weaker claim than it looks.
#
# CALIBRATION / VALIDATION ONLY. Changes nothing in production.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
cd "$W"
CAL="$W/results/ranking/diagnostics/molluscan_calibration"

python3 - "$CAL" <<'PY'
import sys
import pandas as pd

CAL = sys.argv[1]
a = pd.read_csv(f"{CAL}/molluscan_null_proteinclip3b_candidates.tsv",
                sep="\t").set_index("id")
b = pd.read_csv(f"{CAL}/molluscan_null_protrek_candidates.tsv",
                sep="\t").set_index("id")
shared = a.index.intersection(b.index)
print(f"candidates: pc3b={len(a)} protrek={len(b)} shared={len(shared)}")
if len(shared) != len(a) or len(shared) != len(b):
    raise SystemExit("FATAL: candidate id sets differ between models")

for nm, df in (("proteinclip3b", a), ("protrek", b)):
    v = df[df.inside_vert_prodbasis]
    m = df[df.pct_moll_mollusca_all <= 5]
    print(f"\n{nm}:")
    print(f"  vertebrate-envelope inside (production basis) = {len(v)}")
    print(f"    by family: {v.best_family.value_counts().to_dict()}")
    print(f"  molluscan-calibrated inside (q<=5)            = {len(m)}")
    print(f"    by family: {m.best_family.value_counts().to_dict()}")
    kept = int(v.index.isin(m.index).sum())
    print(f"  of the {len(v)} vertebrate-envelope calls, {kept} survive "
          f"molluscan calibration ({100.0 * kept / max(len(v), 1):.1f}%)")

va, vb = set(a[a.inside_vert_prodbasis].index), set(b[b.inside_vert_prodbasis].index)
ma = set(a[a.pct_moll_mollusca_all <= 5].index)
mb = set(b[b.pct_moll_mollusca_all <= 5].index)
print(f"\nCROSS-MODEL CONSENSUS")
print(f"  vertebrate-envelope inside, BOTH models : {len(va & vb):>3}  "
      f"(union {len(va | vb)})")
print(f"  molluscan-calibrated inside, BOTH models: {len(ma & mb):>3}  "
      f"(union {len(ma | mb)})")

both = sorted(ma & mb)
agree = sum(1 for i in both if a.loc[i, "best_family"] == b.loc[i, "best_family"])
print(f"\n  molluscan-calibrated CONSENSUS exclusions (n={len(both)}), "
      f"family agreement {agree}/{len(both)}:")
for i in both:
    fa, fb = a.loc[i, "best_family"], b.loc[i, "best_family"]
    print(f"    {i:<28} pc3b {fa:<22} {a.loc[i, 'pct_moll_mollusca_all']:6.2f} | "
          f"protrek {fb:<22} {b.loc[i, 'pct_moll_mollusca_all']:6.2f} | "
          f"{'AGREE' if fa == fb else 'DISAGREE'}")

out = pd.DataFrame({
    "id": both,
    "pc3b_family": [a.loc[i, "best_family"] for i in both],
    "pc3b_pct_moll": [a.loc[i, "pct_moll_mollusca_all"] for i in both],
    "protrek_family": [b.loc[i, "best_family"] for i in both],
    "protrek_pct_moll": [b.loc[i, "pct_moll_mollusca_all"] for i in both],
})
out["family_agree"] = out.pc3b_family == out.protrek_family
out.to_csv(f"{CAL}/molluscan_calibrated_consensus_exclusions.tsv",
           sep="\t", index=False)
print(f"\nwrote {CAL}/molluscan_calibrated_consensus_exclusions.tsv")
PY

echo "[consensus] DONE $(date -Is)"
