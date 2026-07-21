#!/bin/bash
#SBATCH --job-name=orthodb_og2rec
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Assign the RECOVERED anchor genes to orthogroups, so the mapping rate and the
# cross-clade concordance can be recomputed on the enlarged anchor set rather
# than asserted to be unchanged.
#
# Usage:  sbatch scripts/unity/orthodb_og2genes_recovered.sh
set -eo pipefail

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
CACHE="${CACHE:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/orthodb_v12v2}"
OUT="${OUT:-${REPO}/results/ranking/diagnostics/orthodb}"
cd "${REPO}"; mkdir -p logs

MAP="${OUT}/anchor_gene_map_recovered.tsv"
[ -s "${MAP}" ] || { echo "ERROR: ${MAP} missing (run scripts/orthodb_verify_recovery.py)" >&2; exit 1; }

cut -f2 "${MAP}" | sort -u > "${OUT}/_recovered_gene_ids.txt"
echo "[$(date +%T)] gene ids (original + recovered): $(wc -l < "${OUT}/_recovered_gene_ids.txt")"

# The exact-sequence route reads gene ids out of the Greifswald BRAKER
# partition of OrthoDB, NOT out of the odb12v2 tables. Those two could be
# different releases sharing an id FORMAT while meaning different genes, which
# would make every route-C recovery silently wrong. Confirm each recovered id
# actually exists in the odb12v2 genes table before anything downstream uses it.
awk -F'\t' '$3 ~ /^RECOVERED/ {print $2}' "${MAP}" | sort -u > "${OUT}/_routeC_ids.txt"
if [ -s "${OUT}/_routeC_ids.txt" ]; then
    zcat "${CACHE}/odb12v2_genes.tab.gz" \
      | awk -F'\t' -v idf="${OUT}/_routeC_ids.txt" '
          BEGIN { while ((getline g < idf) > 0) { gsub(/[ \t\r]/,"",g); if (g!="") want[g]=1 } }
          ($1 in want) { print $1 "\t" $2 "\t" $5 }' \
      > "${OUT}/recovered_id_space_check.tsv"
    n_want=$(wc -l < "${OUT}/_routeC_ids.txt")
    n_have=$(cut -f1 "${OUT}/recovered_id_space_check.tsv" | sort -u | wc -l)
    echo "[$(date +%T)] ID-SPACE CHECK: ${n_have}/${n_want} recovered gene ids exist in odb12v2_genes"
    if [ "${n_have}" -eq 0 ]; then
        echo "ERROR: none of the recovered gene ids exist in odb12v2. The FASTA used" >&2
        echo "       by the sequence route is a different release; recoveries are void." >&2
        exit 1
    fi
fi

zcat "${CACHE}/odb12v2_OG2genes.tab.gz" \
  | awk -F'\t' -v gf="${OUT}/_recovered_gene_ids.txt" '
      BEGIN { while ((getline g < gf) > 0) { gsub(/[ \t\r]/,"",g); if (g!="") want[g]=1 } }
      ($2 in want) { print $1 "\t" $2 }' \
  > "${OUT}/anchor_og2genes_recovered.tsv"
echo "[$(date +%T)] rows: $(wc -l < "${OUT}/anchor_og2genes_recovered.tsv")"

cut -f1 "${OUT}/anchor_og2genes_recovered.tsv" | sort -u > "${OUT}/_og_ids_recovered.txt"
zcat "${CACHE}/odb12v2_OGs.tab.gz" \
  | awk -F'\t' -v of="${OUT}/_og_ids_recovered.txt" '
      BEGIN { while ((getline o < of) > 0) { gsub(/[ \t\r]/,"",o); if (o!="") want[o]=1 } }
      ($1 in want) { print }' \
  > "${OUT}/anchor_og_meta_recovered.tsv"
echo "[$(date +%T)] OG metadata rows: $(wc -l < "${OUT}/anchor_og_meta_recovered.tsv")"

rm -f "${OUT}/_recovered_gene_ids.txt" "${OUT}/_og_ids_recovered.txt"
echo "[$(date +%T)] DONE"
