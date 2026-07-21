#!/bin/bash
#SBATCH --job-name=orthodb_verify
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Extract the harvested sequences and verify class A PER ENTRY.
#
# Class A is NOT inherited from the seeding orthogroup. An orthogroup seeded by
# five characterized class-A anchors can still contain fragments, misassemblies
# and non-rhodopsin proteins, and the whole point of the exercise is to widen a
# reference envelope without contaminating it. So every harvested sequence is
# scanned against the Pfam rhodopsin-like 7TM model itself.
#
# Also re-derives each sequence's organism from its OrthoDB gene id prefix and
# checks it against the organism recorded during selection. A gene id that
# resolves to a different organism than the selection table claims means the
# join drifted, and the entry is dropped rather than reported.
#
# Usage:  sbatch scripts/unity/orthodb_verify_harvest.sh
set -eo pipefail

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
OUT="${OUT:-${REPO}/results/ranking/diagnostics/orthodb}"
METAZOA_FA="${REPO}/species_tree_data/orthodb/odb12_metazoa.fa"
HMMDIR="${OUT}/hmm"

cd "${REPO}"
mkdir -p "${OUT}" "${HMMDIR}" logs

SEL="${OUT}/harvest_selection.tsv"
[ -s "${SEL}" ] || { echo "ERROR: ${SEL} missing" >&2; exit 1; }
echo "[$(date +%T)] harvest entries: $(( $(wc -l < "${SEL}") - 1 ))"

source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u
command -v hmmsearch >/dev/null || { echo "ERROR: hmmsearch not on PATH" >&2; exit 1; }

# ------------------------------------------------- extract the sequences --
awk -F'\t' 'NR>1{print $3}' "${SEL}" | sort -u > "${OUT}/_harvest_gene_ids.txt"
echo "[$(date +%T)] distinct gene ids to extract: $(wc -l < "${OUT}/_harvest_gene_ids.txt")"

python3 - "${METAZOA_FA}" "${OUT}/_harvest_gene_ids.txt" "${OUT}/harvest_sequences.fasta" <<'PY'
import sys
fa, idf, out = sys.argv[1], sys.argv[2], sys.argv[3]
want = {l.strip() for l in open(idf) if l.strip()}
found, name, buf = {}, None, []
def flush():
    if name is not None and name in want:
        found[name] = "".join(buf)
with open(fa) as fh:
    for line in fh:
        if line.startswith(">"):
            flush()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
flush()
with open(out, "w") as fh:
    for k, v in found.items():
        fh.write(f">{k}\n{v}\n")
missing = want - set(found)
print(f"extracted {len(found)}/{len(want)} sequences", file=sys.stderr)
if missing:
    print(f"WARNING: {len(missing)} gene ids absent from the FASTA "
          f"(e.g. {sorted(missing)[:5]})", file=sys.stderr)
PY
_n_extracted=$(grep -c '^>' "${OUT}/harvest_sequences.fasta" 2>/dev/null || true)
echo "[$(date +%T)] extracted: ${_n_extracted:-0}"

# ------------------------------------------------------- the 7TM_1 model --
if [ ! -s "${HMMDIR}/PF00001.hmm" ]; then
    echo "[$(date +%T)] fetching Pfam PF00001 (7tm_1)"
    curl -sL --retry 4 --max-time 300 \
        -o "${HMMDIR}/PF00001.hmm.gz" \
        "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00001/?annotation=hmm" \
        && gunzip -f "${HMMDIR}/PF00001.hmm.gz" || true
fi
[ -s "${HMMDIR}/PF00001.hmm" ] || { echo "ERROR: could not obtain PF00001.hmm" >&2; exit 1; }
grep -m1 '^NAME' "${HMMDIR}/PF00001.hmm"

# --------------------------------------------------------- scan per entry --
echo "[$(date +%T)] scanning harvested sequences against 7tm_1"
hmmsearch --cpu 8 -E 1e-5 --domE 1e-5 \
    --tblout "${OUT}/harvest_7tm1.tbl" \
    "${HMMDIR}/PF00001.hmm" "${OUT}/harvest_sequences.fasta" > "${OUT}/harvest_7tm1.out"
echo "[$(date +%T)] hits: $(grep -vc '^#' "${OUT}/harvest_7tm1.tbl" || echo 0)"

# ------------------------------------------------------------- verdicts --
python3 - "${SEL}" "${OUT}/harvest_7tm1.tbl" "${OUT}/harvest_sequences.fasta" \
         "${OUT}/harvest_verified.tsv" "${OUT}/odb_species.tsv" <<'PY'
import csv, re, sys
sel, tbl, fa, out, species = sys.argv[1:6]

# Resolve organism id -> taxid through the OrthoDB species table rather than
# parsing the taxid out of the organism id string. OrthoDB mints an organism id
# once and keeps it, so an organism whose NCBI taxid was later MERGED still
# carries the pre-merge id: Exaiptasia diaphana is organism 1720309_0 while its
# current taxid is 2652724. String-parsing the prefix reports that as a
# cross-organism mismatch and silently discards good cnidarian sequences.
org2taxid = {}
with open(species, newline="") as fh:
    for row in csv.reader(fh, delimiter="\t"):
        if len(row) >= 2:
            org2taxid[row[1]] = row[0]

hits = {}
with open(tbl) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) > 5:
            hits[f[0]] = {"evalue": f[4], "score": f[5]}

seqs, name = {}, None
with open(fa) as fh:
    for line in fh:
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}

GID = re.compile(r"^(\d+_\d+):")
rows, kept = [], 0
with open(sel, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        gid = r["gene_id"]
        seq = seqs.get(gid, "")
        m = GID.match(gid)
        org = m.group(1) if m else ""
        derived_taxid = org2taxid.get(org, "")
        organism_ok = (
            bool(org)
            and org == r["org_id"]          # same organism record
            and derived_taxid == r["taxid"]  # and it still resolves to the same taxid
        )
        is_a = gid in hits
        verdict = ("verified_class_A" if (is_a and organism_ok and seq)
                   else "dropped")
        reasons = []
        if not seq:
            reasons.append("no_sequence")
        if not is_a:
            reasons.append("no_7tm1_hit")
        if not organism_ok:
            reasons.append(f"organism_mismatch({derived_taxid}!={r['taxid']})")
        rows.append({**r,
                     "seq_len": len(seq),
                     "pf00001_evalue": hits.get(gid, {}).get("evalue", ""),
                     "pf00001_score": hits.get(gid, {}).get("score", ""),
                     "organism_verified": str(organism_ok),
                     "verdict": verdict,
                     "drop_reasons": ";".join(reasons),
                     "sequence": seq})
        kept += verdict == "verified_class_A"

with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

from collections import Counter
print(f"harvest entries: {len(rows)}   verified class A: {kept}")
print("drop reasons:", dict(Counter(
    x for r in rows if r["verdict"] != "verified_class_A"
    for x in r["drop_reasons"].split(";") if x)))
print("verified per phylum:", dict(Counter(
    r["phylum"] for r in rows if r["verdict"] == "verified_class_A").most_common()))
PY

rm -f "${OUT}/_harvest_gene_ids.txt"
echo "[$(date +%T)] DONE"
