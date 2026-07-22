#!/bin/bash
#SBATCH --job-name=taar_impact_rebuild
#SBATCH --partition=cpu
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# taar_impact_rebuild.sh — RIGOROUS comparison for bead 1h4c.
#
# QUESTION: does including the 36 vertebrate TAAR2-9 receptors in the
# 'aminergic' non-chemoreceptor reference change ANY Berghia class-A
# candidate classification?
#
# The local DIRECT DIAGNOSTIC already showed the divergence hypothesis is
# confirmed: a pure-TAAR profile's best hit to any of the 790 candidates is
# E~7e-76, ~29 orders of magnitude WEAKER than the aminergic operating
# threshold 2.4e-105 -> TAARs cannot INDEPENDENTLY drive an aminergic call.
# BUT one candidate (BersteEVm003482t1) sits essentially AT the threshold,
# and the LOO threshold is a whole-library property that is recomputed when
# the aminergic training set changes. So the definitive call-change answer
# needs this rebuild. Only run this if the direct diagnostic was NOT deemed
# decisive for the call-change question (it is decisive for divergence).
#
# WHY the production scorer, not a bare hmmsearch: the LOO threshold
# (2.4e-105) is calibrated under hmmscan-against-the-library (search space
# Z = number of HMMs), and production classification (classify_via_hmm.py)
# scores candidates the same way. A bare `hmmsearch aminergic.hmm vs 790`
# uses Z=790 and is on a DIFFERENT E-value scale (~79x), so it cannot be
# compared to 2.4e-105 directly. This script therefore reuses the exact
# production trio (build -> LOO validate -> classify) for BOTH arms so the
# only variable is the presence/absence of the 36 TAARs.
#
# ARMS:
#   FULL    = all 165 aminergic refs (current production state; 36 TAARs IN)
#   REDUCED = 129 aminergic refs (36 TAARs REMOVED)
# Both libraries are rebuilt FRESH on this node so the two arms are produced
# by identical code; the FULL arm should reproduce the committed 2.4e-105.
#
# Submit:  sbatch scripts/unity/taar_impact_rebuild.sh
# Do NOT submit from a login node's compute path; this is the sbatch itself.

set -eo pipefail
source ~/.miniconda3/etc/profile.d/conda.sh
conda activate berghia-gpcr
set -u

# ---- resolve repo root from this script's location ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

REF_TSV="${REPO_ROOT}/references/non_chemo_gpcr/all_references.tsv"
REF_FASTA="${REPO_ROOT}/references/non_chemo_gpcr/all_references.fasta"
CAND_FASTA="${CAND_FASTA:-${REPO_ROOT}/reports/plm_report/data/chemogpcrs_berghia_classA.fa}"
BUILD_PY="${REPO_ROOT}/scripts/build_classification_hmms.py"
VALIDATE_PY="${REPO_ROOT}/scripts/validate_classification_hmms.py"
CLASSIFY_PY="${REPO_ROOT}/scripts/classify_via_hmm.py"

OUT="${REPO_ROOT}/results/classification/taar_impact/rigorous"
mkdir -p "${OUT}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

for f in "${REF_TSV}" "${REF_FASTA}" "${CAND_FASTA}" "${BUILD_PY}" "${VALIDATE_PY}" "${CLASSIFY_PY}"; do
  [ -s "$f" ] || { echo "FATAL: missing input $f" >&2; exit 3; }
done
echo "[taar-impact] repo=${REPO_ROOT}"
echo "[taar-impact] candidates=$(grep -c '^>' "${CAND_FASTA}")  (expect 790)"

# ======================================================================
# STEP 1 — derive the 36 TAAR2-9 accessions (deterministic rule) and write
#          a REDUCED reference set (TSV + FASTA) with them removed.
#          Fails loud unless exactly 36 are selected and aminergic 165->129.
# ======================================================================
REDUCED_TSV="${OUT}/all_references_reduced.tsv"
REDUCED_FASTA="${OUT}/all_references_reduced.fasta"
EXCL_LIST="${OUT}/taar36_accessions.txt"

python3 - "$REF_TSV" "$REF_FASTA" "$REDUCED_TSV" "$REDUCED_FASTA" "$EXCL_LIST" <<'PY'
import re, sys
ref_tsv, ref_fasta, red_tsv, red_fasta, excl_out = sys.argv[1:6]

# 3 vertebrate TAAR1 orthologs + 2 genuine invertebrate tyramine receptors
EXCLUDE_KNOWN = {"Q96RJ0", "Q923Y8", "Q923Y9", "O02213", "Q19084"}
gene_re = re.compile("TAAR", re.IGNORECASE)

rows = []
with open(ref_tsv) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(header)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= len(header):
            rows.append(f)

def g(r, c): return r[idx[c]]

tyramine = [r for r in rows if g(r, "family") == "aminergic" and g(r, "subfamily") == "tyramine"]
taar_gene = [r for r in tyramine if gene_re.search(g(r, "gene"))]
taar36 = [g(r, "accession") for r in taar_gene if g(r, "accession") not in EXCLUDE_KNOWN]

# assertions on the split
present = {g(r, "accession") for r in tyramine} & EXCLUDE_KNOWN
assert present == EXCLUDE_KNOWN, f"exclusion set not fully present: missing {EXCLUDE_KNOWN - present}"
assert len(taar36) == len(set(taar36)), "duplicate TAAR accession"
if len(taar36) != 36:
    sys.stderr.write(f"FATAL: selected {len(taar36)} TAARs, expected 36:\n")
    for a in sorted(taar36):
        sys.stderr.write(f"  {a}\n")
    sys.exit(2)

taar_set = set(taar36)
with open(excl_out, "w") as out:
    for a in sorted(taar_set):
        out.write(a + "\n")

# aminergic counts before/after
amin_all = [r for r in rows if g(r, "family") == "aminergic"]
amin_reduced = [r for r in amin_all if g(r, "accession") not in taar_set]
assert len(amin_all) == 165, f"expected 165 aminergic, found {len(amin_all)}"
assert len(amin_reduced) == 129, f"expected 129 aminergic after removal, found {len(amin_reduced)}"

# write reduced TSV
with open(red_tsv, "w") as out:
    out.write("\t".join(header) + "\n")
    kept = 0
    for r in rows:
        if g(r, "accession") in taar_set:
            continue
        out.write("\t".join(r) + "\n")
        kept += 1
assert kept == len(rows) - 36, f"reduced TSV kept {kept}, expected {len(rows)-36}"

# write reduced FASTA (header '>ACC|family|subfamily|species'; acc = first '|' field)
def acc_of(hdr):
    return hdr[1:].split("|", 1)[0].strip()

removed = 0
with open(ref_fasta) as fin, open(red_fasta, "w") as fout:
    write = True
    for line in fin:
        if line.startswith(">"):
            a = acc_of(line.rstrip("\n"))
            write = a not in taar_set
            if not write:
                removed += 1
        if write:
            fout.write(line)
assert removed == 36, f"removed {removed} FASTA records, expected 36"

sys.stderr.write(f"[step1] OK: 36 TAARs selected; aminergic 165->129; "
                 f"reduced TSV rows={kept}, FASTA removed={removed}\n")
PY
echo "[taar-impact] step 1 done: reduced reference written"

# ======================================================================
# STEP 2 — build BOTH HMM libraries FRESH (reuse build_classification_hmms.py)
# ======================================================================
HMMS_FULL="${OUT}/hmms_full"
HMMS_REDUCED="${OUT}/hmms_reduced"
rm -rf "${HMMS_FULL}" "${HMMS_REDUCED}"

python3 "${BUILD_PY}" \
  --reference-fasta "${REF_FASTA}" \
  --reference-tsv   "${REF_TSV}" \
  --output-dir      "${HMMS_FULL}" \
  --threads "${THREADS}"

python3 "${BUILD_PY}" \
  --reference-fasta "${REDUCED_FASTA}" \
  --reference-tsv   "${REDUCED_TSV}" \
  --output-dir      "${HMMS_REDUCED}" \
  --threads "${THREADS}"

# verify aminergic member counts in each arm
NF=$(grep -m1 '^NSEQ' "${HMMS_FULL}/aminergic.hmm"    | awk '{print $2}')
NR=$(grep -m1 '^NSEQ' "${HMMS_REDUCED}/aminergic.hmm" | awk '{print $2}')
echo "[taar-impact] aminergic NSEQ: full=${NF} reduced=${NR}"
[ "${NF}" = "165" ] || { echo "FATAL: full aminergic NSEQ=${NF}, expected 165" >&2; exit 4; }
[ "${NR}" = "129" ] || { echo "FATAL: reduced aminergic NSEQ=${NR}, expected 129" >&2; exit 4; }

# ======================================================================
# STEP 3 — recompute LOO thresholds for BOTH libraries (reuse validator).
#          The reduced arm's aminergic threshold is what changes.
# ======================================================================
LOO_FULL="${OUT}/loo_full"
LOO_REDUCED="${OUT}/loo_reduced"

python3 "${VALIDATE_PY}" --hmm-dir "${HMMS_FULL}"    --output-dir "${LOO_FULL}"    --threads "${THREADS}"
python3 "${VALIDATE_PY}" --hmm-dir "${HMMS_REDUCED}" --output-dir "${LOO_REDUCED}" --threads "${THREADS}"

THR_FULL=$(awk -F'\t' '$1=="aminergic"{print $NF}' "${LOO_FULL}/loo_metrics.tsv")
THR_REDUCED=$(awk -F'\t' '$1=="aminergic"{print $NF}' "${LOO_REDUCED}/loo_metrics.tsv")
echo "[taar-impact] aminergic LOO threshold: full=${THR_FULL} reduced=${THR_REDUCED}"

# ======================================================================
# STEP 4 — classify the 790 candidates with EACH library at its OWN LOO
#          thresholds (production scorer: hmmscan vs consolidated library).
# ======================================================================
CALLS_FULL="${OUT}/calls_full.tsv"
CALLS_REDUCED="${OUT}/calls_reduced.tsv"

python3 "${CLASSIFY_PY}" \
  --candidate-fasta "${CAND_FASTA}" \
  --hmm-dir "${HMMS_FULL}" \
  --loo-metrics "${LOO_FULL}/loo_metrics.tsv" \
  --output-tsv "${CALLS_FULL}" \
  --threads "${THREADS}"

python3 "${CLASSIFY_PY}" \
  --candidate-fasta "${CAND_FASTA}" \
  --hmm-dir "${HMMS_REDUCED}" \
  --loo-metrics "${LOO_REDUCED}/loo_metrics.tsv" \
  --output-tsv "${CALLS_REDUCED}" \
  --threads "${THREADS}"

# ======================================================================
# STEP 5 — DIFF the per-candidate calls between the two arms. Any candidate
#          whose family / subfamily / pass-status changes is what the whole
#          test is about. Report coverage (N compared) and itemise changes.
# ======================================================================
COMPARISON="${OUT}/comparison.tsv"
python3 - "$CALLS_FULL" "$CALLS_REDUCED" "$COMPARISON" "$THR_FULL" "$THR_REDUCED" <<'PY'
import csv, sys
full_p, red_p, out_p, thr_full, thr_red = sys.argv[1:6]

def load(p):
    d = {}
    with open(p) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            d[row["candidate_id"]] = row
    return d

F, R = load(full_p), load(red_p)
ids = sorted(set(F) | set(R))
assert F.keys() == R.keys(), "candidate sets differ between arms — investigate"

n_changed = 0
changed_rows = []
with open(out_p, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["candidate_id",
                "family_full", "family_reduced",
                "subfamily_full", "subfamily_reduced",
                "evalue_full", "evalue_reduced",
                "threshold_source_full", "threshold_source_reduced",
                "changed"])
    for cid in ids:
        f, r = F[cid], R[cid]
        changed = (f["hmm_family"] != r["hmm_family"] or
                   f["hmm_subfamily"] != r["hmm_subfamily"])
        if changed:
            n_changed += 1
            changed_rows.append((cid, f["hmm_family"], r["hmm_family"]))
        w.writerow([cid,
                    f["hmm_family"], r["hmm_family"],
                    f["hmm_subfamily"], r["hmm_subfamily"],
                    f["evalue"], r["evalue"],
                    f["threshold_source"], r["threshold_source"],
                    "YES" if changed else "no"])

sys.stderr.write("\n==================== TAAR-IMPACT VERDICT ====================\n")
sys.stderr.write(f"aminergic LOO threshold: full(165)={thr_full}  reduced(129)={thr_red}\n")
sys.stderr.write(f"candidates compared: {len(ids)}\n")
sys.stderr.write(f"candidates whose classification CHANGED: {n_changed}\n")
for cid, ff, rf in changed_rows:
    sys.stderr.write(f"  {cid}: {ff}  ->  {rf}\n")
if n_changed == 0:
    sys.stderr.write("VERDICT: including the 36 TAARs changes NO candidate call "
                     "-> inclusion is inconsequential (hypothesis confirmed).\n")
else:
    sys.stderr.write("VERDICT: including the 36 TAARs DOES change the calls above "
                     "-> inclusion is consequential; see comparison.tsv.\n")
sys.stderr.write("=============================================================\n")
PY

echo "[taar-impact] DONE. Outputs in ${OUT}:"
echo "  comparison.tsv        (per-candidate diff, the primary result)"
echo "  calls_full.tsv / calls_reduced.tsv"
echo "  loo_full/loo_metrics.tsv / loo_reduced/loo_metrics.tsv (thresholds)"
