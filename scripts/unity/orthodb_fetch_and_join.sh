#!/bin/bash
#SBATCH --job-name=orthodb_join
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Fetch the OrthoDB v12v2 flat-file dump tables needed to map curated GPCR
# anchors onto orthogroups, verify them against the checksums OrthoDB
# publishes, and reduce them to compact per-anchor / per-orthogroup tables.
#
# We deliberately do NOT touch species_tree_data/orthodb/odb12_metazoa.fa. That
# file is the Greifswald BRAKER partition of OrthoDB and carries no orthogroup
# assignment at all, so it cannot answer this question; and it is 21.9 GB.
#
# Downloads (~9.4 GB, one time, cached and skipped on re-run):
#   odb12v2_levels.tab.gz         16 kB   OG levels (clades)
#   odb12v2_species.tab.gz       647 kB   organisms
#   odb12v2_level2species.tab.gz 247 kB   organism -> level path
#   odb12v2_OGs.tab.gz           128 MB   OG -> level, name
#   odb12v2_gene_xrefs.tab.gz    4.7 GB   OrthoDB gene -> UniProt etc.
#   odb12v2_OG2genes.tab.gz      4.5 GB   OG -> OrthoDB gene
#
# Checksums are scraped from the OrthoDB dump page at runtime rather than
# hardcoded, so a silently-updated release fails loudly instead of being
# joined against stale expectations.
#
# Usage:  sbatch scripts/unity/orthodb_fetch_and_join.sh
set -eo pipefail

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
CACHE="${CACHE:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/orthodb_v12v2}"
OUT="${OUT:-${REPO}/results/ranking/diagnostics/orthodb}"
DUMP_PAGE="https://data.orthodb.org/current/download/odb_data_dump"
BASE="https://data.orthodb.org/current/download/odb_data_dump"

mkdir -p "${CACHE}" "${OUT}" "${REPO}/logs"

ACC="${OUT}/anchor_accessions.txt"
if [ ! -s "${ACC}" ]; then
    echo "ERROR: anchor accession list missing at ${ACC}." >&2
    echo "Run scripts/orthodb_snapshot_anchors.py and copy its outputs here first." >&2
    exit 1
fi
echo "[$(date +%T)] anchors to map: $(wc -l < "${ACC}")"

# ---------------------------------------------------------------- checksums --
echo "[$(date +%T)] scraping published checksums from ${DUMP_PAGE}"
curl -sL --max-time 120 "${DUMP_PAGE}" -o "${CACHE}/dump_index.html"
python3 - "${CACHE}/dump_index.html" "${CACHE}/published_md5.tsv" <<'PY'
import re, sys
html = open(sys.argv[1]).read()
# Each table row: <td><a href=...>FILE</a></td><td>SIZE</td><td>DESC</td><td>MD5</td>
rows = re.findall(
    r'<a href=\S+?/(odb\S+?)\s[^>]*>\1</a></td>\s*<td>(.*?)</td>\s*<td>.*?</td>\s*<td>([0-9a-f]{32})</td>',
    html, re.S)
if not rows:
    sys.exit("ERROR: could not parse any file/md5 rows from the OrthoDB dump page")
with open(sys.argv[2], 'w') as fh:
    for name, size, md5 in rows:
        fh.write(f"{name}\t{size.strip()}\t{md5}\n")
print(f"parsed {len(rows)} published checksums")
PY
cat "${CACHE}/published_md5.tsv"

FILES=(
    odb12v2_levels.tab.gz
    odb12v2_species.tab.gz
    odb12v2_level2species.tab.gz
    odb12v2_OGs.tab.gz
    odb12v2_gene_xrefs.tab.gz
    odb12v2_OG2genes.tab.gz
)

# ----------------------------------------------------------------- download --
for f in "${FILES[@]}"; do
    want=$(awk -F'\t' -v n="$f" '$1==n{print $3}' "${CACHE}/published_md5.tsv")
    if [ -z "${want}" ]; then
        echo "ERROR: ${f} has no published checksum on the dump page" >&2; exit 1
    fi
    if [ -s "${CACHE}/${f}" ]; then
        have=$(md5sum "${CACHE}/${f}" | cut -d' ' -f1)
        if [ "${have}" = "${want}" ]; then
            echo "[$(date +%T)] ${f} cached and checksum-verified"; continue
        fi
        echo "[$(date +%T)] ${f} cached but checksum mismatch -> re-downloading"
    fi
    echo "[$(date +%T)] downloading ${f}"
    curl -sL --retry 5 --retry-delay 10 --max-time 7200 -o "${CACHE}/${f}" "${BASE}/${f}"
    have=$(md5sum "${CACHE}/${f}" | cut -d' ' -f1)
    if [ "${have}" != "${want}" ]; then
        echo "ERROR: checksum mismatch for ${f}: got ${have}, published ${want}" >&2; exit 1
    fi
    echo "[$(date +%T)] ${f} downloaded and checksum-verified"
done

# ------------------------------------------- step 1: UniProt -> OrthoDB gene --
# gene_xrefs columns: 1 OrthoDB gene id, 2 external id, 3 external DB name.
echo "[$(date +%T)] step 1: mapping anchor UniProt accessions to OrthoDB gene ids"
zcat "${CACHE}/odb12v2_gene_xrefs.tab.gz" \
  | awk -F'\t' -v accfile="${ACC}" '
      BEGIN { while ((getline a < accfile) > 0) { gsub(/[ \t\r]/,"",a); if (a!="") want[a]=1 } }
      $3=="UniProt" {
          acc=$2; sub(/[-.][0-9]+$/,"",acc)      # strip isoform/version suffix
          if (acc in want) print acc "\t" $1 "\t" $2
      }' \
  > "${OUT}/anchor_gene_map.tsv"
echo "[$(date +%T)] step 1 rows: $(wc -l < "${OUT}/anchor_gene_map.tsv")"
echo "[$(date +%T)] step 1 distinct anchors mapped: $(cut -f1 "${OUT}/anchor_gene_map.tsv" | sort -u | wc -l)"

if [ ! -s "${OUT}/anchor_gene_map.tsv" ]; then
    echo "ERROR: zero anchors mapped to OrthoDB gene ids. The join key is wrong -- aborting" >&2
    echo "       rather than emitting statistics computed on an empty intersection." >&2
    exit 1
fi

cut -f2 "${OUT}/anchor_gene_map.tsv" | sort -u > "${OUT}/anchor_gene_ids.txt"

# ------------------------------------------------ step 2: anchor gene -> OG --
# OG2genes columns: 1 OG id, 2 OrthoDB gene id.
echo "[$(date +%T)] step 2: assigning anchor genes to orthogroups at every level"
zcat "${CACHE}/odb12v2_OG2genes.tab.gz" \
  | awk -F'\t' -v gf="${OUT}/anchor_gene_ids.txt" '
      BEGIN { while ((getline g < gf) > 0) { gsub(/[ \t\r]/,"",g); if (g!="") want[g]=1 } }
      ($2 in want) { print $1 "\t" $2 }' \
  > "${OUT}/anchor_og2genes.tsv"
echo "[$(date +%T)] step 2 rows: $(wc -l < "${OUT}/anchor_og2genes.tsv")"

cut -f1 "${OUT}/anchor_og2genes.tsv" | sort -u > "${OUT}/anchor_og_ids.txt"
echo "[$(date +%T)] distinct anchor-containing OGs: $(wc -l < "${OUT}/anchor_og_ids.txt")"

# ------------------------------------------------- step 3: OG metadata rows --
# OGs columns: 1 OG id, 2 level tax id, 3 OG name.
echo "[$(date +%T)] step 3: extracting OG metadata"
zcat "${CACHE}/odb12v2_OGs.tab.gz" \
  | awk -F'\t' -v of="${OUT}/anchor_og_ids.txt" '
      BEGIN { while ((getline o < of) > 0) { gsub(/[ \t\r]/,"",o); if (o!="") want[o]=1 } }
      ($1 in want) { print }' \
  > "${OUT}/anchor_og_meta.tsv"
echo "[$(date +%T)] step 3 rows: $(wc -l < "${OUT}/anchor_og_meta.tsv")"

# ------------------------------- step 4: full membership of anchor-bearing OGs --
# Needed for the expansion ceiling: every gene in every OG that contains a
# characterized anchor, so members can be counted per clade.
echo "[$(date +%T)] step 4: full membership of anchor-bearing orthogroups"
zcat "${CACHE}/odb12v2_OG2genes.tab.gz" \
  | awk -F'\t' -v of="${OUT}/anchor_og_ids.txt" '
      BEGIN { while ((getline o < of) > 0) { gsub(/[ \t\r]/,"",o); if (o!="") want[o]=1 } }
      ($1 in want) { print $1 "\t" $2 }' \
  > "${OUT}/anchor_og_membership.tsv"
echo "[$(date +%T)] step 4 rows: $(wc -l < "${OUT}/anchor_og_membership.tsv")"

# ------------------------------------------------ step 5: reference metadata --
echo "[$(date +%T)] step 5: copying levels / species / level2species"
zcat "${CACHE}/odb12v2_levels.tab.gz"        > "${OUT}/odb_levels.tsv"
zcat "${CACHE}/odb12v2_species.tab.gz"       > "${OUT}/odb_species.tsv"
zcat "${CACHE}/odb12v2_level2species.tab.gz" > "${OUT}/odb_level2species.tsv"

# Definitive presence/absence check for the target species, against the
# authoritative organism table rather than a search endpoint.
{
  echo -e "check\tvalue"
  echo -e "orthodb_release\todb12v2"
  echo -e "organisms_total\t$(wc -l < "${OUT}/odb_species.tsv")"
  echo -e "berghia_taxid_1287507_rows\t$(awk -F'\t' '$1==1287507' "${OUT}/odb_species.tsv" | wc -l)"
  echo -e "berghia_name_rows\t$(grep -ci berghia "${OUT}/odb_species.tsv" || true)"
  echo -e "mollusca_6447_organisms\t$(awk -F'\t' '$1==6447{print $5}' "${OUT}/odb_levels.tsv")"
} > "${OUT}/berghia_presence_check.tsv"
cat "${OUT}/berghia_presence_check.tsv"

echo "[$(date +%T)] DONE. Outputs in ${OUT}"
ls -la "${OUT}"
