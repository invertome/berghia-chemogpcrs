#!/bin/bash
#SBATCH --job-name=moll_cal_merge
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Merge the per-shard embedding npz into one npz per model, then project the
# ESM-2 3B base through the released ProteinCLIP esm2_36 head to produce
# proteinclip3b -- exactly the derivation scripts/unity/run_proteinclip3b.sh
# uses for the production channel.
#
# Fails loudly on any id lost between the input FASTA and the merged npz: a
# silently short null would bias every percentile it feeds.
set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate esmc
set -u

W=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
PR=/scratch3/workspace/jperezmoreno_umass_edu-jorge/proteinclip
cd "$W"
CAL="$W/results/ranking/diagnostics/molluscan_calibration"
SH="$CAL/shards"

python3 - "$CAL" "$SH" <<'PY'
import glob, os, sys
import numpy as np

cal, sh = sys.argv[1], sys.argv[2]
expected = set()
with open(os.path.join(cal, "molluscan_null_all.fa")) as fh:
    for line in fh:
        if line.startswith(">"):
            expected.add(line[1:].split()[0])
print(f"[merge] expected ids: {len(expected)}")

for tag in ("esm2_3b", "protrek"):
    parts = sorted(glob.glob(os.path.join(sh, f"{tag}_shard*.npz")))
    if not parts:
        print(f"[merge] no shards for {tag}; skipping")
        continue
    merged, dims = {}, set()
    for p in parts:
        z = np.load(p)
        for k in z.files:
            if k in merged:
                raise SystemExit(f"[merge] FATAL: id {k} appears in >1 shard")
            v = np.asarray(z[k], dtype=np.float32)
            merged[k] = v
            dims.add(v.shape[-1])
    missing = expected - set(merged)
    extra = set(merged) - expected
    print(f"[merge] {tag}: {len(parts)} shards, {len(merged)} ids, dims={sorted(dims)}, "
          f"missing={len(missing)}, extra={len(extra)}")
    if len(dims) != 1:
        raise SystemExit(f"[merge] FATAL: {tag} has mixed dims {sorted(dims)}")
    if missing:
        raise SystemExit(f"[merge] FATAL: {tag} lost {len(missing)} ids "
                         f"(e.g. {sorted(missing)[:3]}); a short null biases "
                         "every percentile computed from it")
    if extra:
        raise SystemExit(f"[merge] FATAL: {tag} has {len(extra)} ids not in the input")
    bad = [k for k, v in merged.items() if not np.isfinite(v).all()]
    if bad:
        raise SystemExit(f"[merge] FATAL: {tag} has {len(bad)} non-finite vectors")
    out = os.path.join(cal, f"molluscan_null_{tag}.npz")
    np.savez(out + ".tmp.npz", **merged)
    os.replace(out + ".tmp.npz", out)
    print(f"[merge] wrote {out}")
PY

BASE="$CAL/molluscan_null_esm2_3b.npz"
PC="$CAL/molluscan_null_proteinclip3b.npz"
if [ -s "$BASE" ] && [ ! -s "$PC" ]; then
    PROTEINCLIP_REPO="$PR" PC_HEAD=36 BASE_NPZ="$BASE" OUT_NPZ="$PC" \
        python3 scratch_proteinclip_project_param.py
fi

echo "[merge] DONE $(date -Is)"
ls -la "$CAL"/molluscan_null_*.npz
