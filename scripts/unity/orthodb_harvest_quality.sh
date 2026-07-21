#!/bin/bash
#SBATCH --job-name=orthodb_qual
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Filter harvest fragments using a threshold CALIBRATED ON THE ANCHORS, not
# picked by eye.
#
# A 7tm_1 hit at E<1e-5 fires on a partial match: a sequence carrying three of
# seven transmembrane helices still hits. 19% of the harvest came in under
# 320 aa, and a 98 aa "opsin" is a fragment, not a reference sequence. Putting
# fragments into a reference envelope corrupts exactly what the envelope is
# supposed to measure - the within-family spread.
#
# So both the anchors and the harvest are scanned with --domtblout, the
# anchors' own distribution of 7tm_1 MODEL COVERAGE and length defines the
# acceptable range, and harvested sequences are required to fall inside it. The
# anchors are curated class-A receptors, so whatever coverage they achieve is
# by definition what a good class-A sequence looks like under this model.
#
# Usage:  sbatch scripts/unity/orthodb_harvest_quality.sh
set -eo pipefail

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
OUT="${OUT:-${REPO}/results/ranking/diagnostics/orthodb}"
HMM="${OUT}/hmm/PF00001.hmm"
cd "${REPO}"; mkdir -p logs

ANCHOR_FA="${OUT}/anchor_classA_working.fasta"
HARVEST_FA="${OUT}/harvest_sequences.fasta"
for f in "${HMM}" "${ANCHOR_FA}" "${HARVEST_FA}"; do
    [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

echo "[$(date +%T)] scanning anchors ($(grep -c '^>' "${ANCHOR_FA}")) and harvest ($(grep -c '^>' "${HARVEST_FA}"))"
hmmsearch --cpu 8 -E 1e-3 --domtblout "${OUT}/anchor_7tm1.domtbl"  "${HMM}" "${ANCHOR_FA}"  > /dev/null
hmmsearch --cpu 8 -E 1e-3 --domtblout "${OUT}/harvest_7tm1.domtbl" "${HMM}" "${HARVEST_FA}" > /dev/null

python3 - "${OUT}/anchor_7tm1.domtbl" "${OUT}/harvest_7tm1.domtbl" \
         "${OUT}/harvest_verified.tsv" "${OUT}/harvest_quality.tsv" <<'PY'
import csv, statistics, sys
from collections import defaultdict

adom, hdom, ver, out = sys.argv[1:5]

def coverage(path):
    """seq -> (best model coverage fraction, seq length). domtblout cols:
    0 target, 2 tlen, 3 query, 5 qlen, 15/16 hmm from/to."""
    best, tlen = defaultdict(float), {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 17:
                continue
            seq, qlen = f[0], int(f[5])
            tlen[seq] = int(f[2])
            cov = (int(f[16]) - int(f[15]) + 1) / qlen
            best[seq] = max(best[seq], cov)
    return {s: (best[s], tlen[s]) for s in best}

acov, hcov = coverage(adom), coverage(hdom)

av = sorted(c for c, _ in acov.values())
al = sorted(l for _, l in acov.values())
def pct(v, p):
    return v[min(len(v) - 1, int(p / 100 * len(v)))]

# The anchors define the acceptable region. p1 is deliberately permissive:
# the goal is to exclude fragments, not to demand the harvest look more
# canonical than curated references do.
COV_MIN, LEN_MIN = pct(av, 1), pct(al, 1)
print(f"anchors scanned: {len(acov)}")
print(f"anchor 7tm_1 coverage  p1={pct(av,1):.3f} p5={pct(av,5):.3f} "
      f"median={statistics.median(av):.3f}")
print(f"anchor length          p1={pct(al,1)} p5={pct(al,5)} "
      f"median={int(statistics.median(al))}")
print(f"=> thresholds: coverage >= {COV_MIN:.3f}, length >= {LEN_MIN}")

rows, kept = [], 0
with open(ver, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        gid = r["gene_id"]
        cov, tlen = hcov.get(gid, (0.0, int(r["seq_len"] or 0)))
        reasons = []
        if r["verdict"] != "verified_class_A":
            reasons.append(r["drop_reasons"])
        if cov < COV_MIN:
            reasons.append(f"low_7tm1_coverage({cov:.2f})")
        if tlen < LEN_MIN:
            reasons.append(f"too_short({tlen})")
        ok = not reasons
        rows.append({**{k: v for k, v in r.items() if k != "sequence"},
                     "tm7_coverage": round(cov, 4),
                     "quality_verdict": "PASS" if ok else "FAIL",
                     "quality_reasons": ";".join(x for x in reasons if x),
                     "sequence": r.get("sequence", "")})
        kept += ok

with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

from collections import Counter
passed = [r for r in rows if r["quality_verdict"] == "PASS"]
print(f"\nharvest entries: {len(rows)}   PASS: {kept}")
print("fail reasons:", dict(Counter(
    x.split("(")[0] for r in rows if r["quality_verdict"] == "FAIL"
    for x in r["quality_reasons"].split(";") if x)))
print("PASS per phylum:", dict(Counter(r["phylum"] for r in passed).most_common()))
print("PASS per family:", dict(Counter(r["family"] for r in passed).most_common()))
anchor_phyla = {"Annelida","Arthropoda","Chordata","Mollusca","Nematoda",
                "Nemertea","Platyhelminthes"}
newp = [r for r in passed if r["phylum"] and r["phylum"] not in anchor_phyla]
print(f"PASS from phyla at zero in the anchor set: {len(newp)} "
      f"across {len({r['phylum'] for r in newp})} phyla")
PY

echo "[$(date +%T)] DONE"
