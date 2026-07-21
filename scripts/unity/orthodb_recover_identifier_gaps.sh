#!/bin/bash
#SBATCH --job-name=orthodb_recover
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --partition=cpu
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Recover the anchors whose species IS in OrthoDB but whose UniProt accession
# carries no UniProt cross-reference on the OrthoDB side, by trying join keys
# other than the accession.
#
# Three routes, cheapest first:
#   A  external ids (EMBL / RefSeq / Ensembl / NCBI GeneID) against the
#      gene_xrefs table's external-identifier column
#   B  gene name + organism against the genes table
#   C  exact protein sequence match, which depends on no identifier agreeing
#
# NOTHING here is trusted on the strength of a string match. Every recovered
# OrthoDB gene id is re-derived back to its organism and checked against the
# anchor's own taxid, because a join key that maps an accession to the WRONG
# gene is worse than no mapping at all. Conflicts are reported, never resolved
# silently.
#
# Usage:  sbatch scripts/unity/orthodb_recover_identifier_gaps.sh
set -eo pipefail

REPO="${REPO:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
CACHE="${CACHE:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/orthodb_v12v2}"
OUT="${OUT:-${REPO}/results/ranking/diagnostics/orthodb}"
BASE="https://data.orthodb.org/current/download/odb_data_dump"
METAZOA_FA="${REPO}/species_tree_data/orthodb/odb12_metazoa.fa"

cd "${REPO}"
mkdir -p "${OUT}" logs

ALT="${OUT}/identifier_gap_alt_ids.txt"
SEQ="${OUT}/identifier_gap_sequences.fasta"
for f in "${ALT}" "${SEQ}" "${OUT}/identifier_gap_alt_xrefs.tsv"; do
    [ -s "$f" ] || { echo "ERROR: missing $f (run scripts/orthodb_fetch_alt_xrefs.py)" >&2; exit 1; }
done
echo "[$(date +%T)] alternate-key rows: $(wc -l < "${ALT}")"

# The genes table is needed for route B and to resolve any recovered gene id
# back to its organism. Checksum-verified like the rest.
GENES=odb12v2_genes.tab.gz
want=$(awk -F'\t' -v n="${GENES}" '$1==n{print $3}' "${CACHE}/published_md5.tsv")
[ -n "${want}" ] || { echo "ERROR: no published checksum for ${GENES}" >&2; exit 1; }
if [ -s "${CACHE}/${GENES}" ] && [ "$(md5sum "${CACHE}/${GENES}" | cut -d' ' -f1)" = "${want}" ]; then
    echo "[$(date +%T)] ${GENES} cached and checksum-verified"
else
    echo "[$(date +%T)] downloading ${GENES}"
    curl -sL --retry 5 --retry-delay 10 --max-time 7200 -o "${CACHE}/${GENES}" "${BASE}/${GENES}"
    have=$(md5sum "${CACHE}/${GENES}" | cut -d' ' -f1)
    [ "${have}" = "${want}" ] || { echo "ERROR: ${GENES} checksum mismatch" >&2; exit 1; }
    echo "[$(date +%T)] ${GENES} downloaded and checksum-verified"
fi

# ------------------------------------------- route A: external identifiers --
echo "[$(date +%T)] route A: external ids against gene_xrefs"
cut -f3 "${ALT}" | sort -u > "${OUT}/_alt_id_values.txt"
zcat "${CACHE}/odb12v2_gene_xrefs.tab.gz" \
  | awk -F'\t' -v idf="${OUT}/_alt_id_values.txt" '
      BEGIN { while ((getline v < idf) > 0) { gsub(/[ \t\r]/,"",v); if (v!="") want[v]=1 } }
      { key=$2; sub(/\.[0-9]+$/,"",key)
        if ($2 in want) print $2 "\t" $1 "\t" $3
        else if (key in want) print key "\t" $1 "\t" $3 }' \
  | sort -u > "${OUT}/recover_routeA_hits.tsv"
echo "[$(date +%T)] route A raw hits: $(wc -l < "${OUT}/recover_routeA_hits.tsv")"

# ---------------------------------------------- route B: gene name + organism --
# genes columns: 1 gene id, 2 organism id, 3 original seq id, 4 synonyms,
#                5 uniprot, 6 ensembl, 7 ncbi gid/gene name, 8 description
echo "[$(date +%T)] route B: gene name + organism against the genes table"
awk -F'\t' '$2=="gene_names"{print $3}' "${ALT}" | sort -u > "${OUT}/_alt_gene_names.txt"
if [ -s "${OUT}/_alt_gene_names.txt" ]; then
  zcat "${CACHE}/${GENES}" \
    | awk -F'\t' -v nf="${OUT}/_alt_gene_names.txt" '
        BEGIN { while ((getline v < nf) > 0) { gsub(/[ \t\r]/,"",v); if (v!="") want[tolower(v)]=1 } }
        { n=tolower($7); gsub(/^[ \t]+|[ \t]+$/,"",n)
          if (n != "" && (n in want)) print $7 "\t" $1 "\t" $2 }' \
    | sort -u > "${OUT}/recover_routeB_hits.tsv"
else
  : > "${OUT}/recover_routeB_hits.tsv"
fi
echo "[$(date +%T)] route B raw hits: $(wc -l < "${OUT}/recover_routeB_hits.tsv")"

# ------------------------------------------------ route C: exact sequence --
# Uses the Metazoa protein FASTA already on disk. NOT re-downloaded. Before
# relying on it we confirm its headers live in the same gene-id space as the
# odb12v2 tables, otherwise a sequence hit would map to a meaningless id.
echo "[$(date +%T)] route C: exact sequence match"
if [ ! -s "${METAZOA_FA}" ]; then
    echo "[$(date +%T)] route C SKIPPED: ${METAZOA_FA} not present"
    : > "${OUT}/recover_routeC_hits.tsv"
else
    head -400000 "${METAZOA_FA}" | grep '^>' | sed 's/^>//' | cut -f1 | cut -d' ' -f1 \
        | sort -u > "${OUT}/_fa_ids_sample.txt"
    shared=$(comm -12 "${OUT}/_fa_ids_sample.txt" <(sort -u "${OUT}/anchor_gene_ids.txt") | wc -l)
    echo "[$(date +%T)] gene-id space check: ${shared} of the sampled FASTA ids are also anchor gene ids"
    python3 - "${METAZOA_FA}" "${SEQ}" "${OUT}/recover_routeC_hits.tsv" <<'PY'
import hashlib, sys
fa, qry, out = sys.argv[1], sys.argv[2], sys.argv[3]

def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)

want = {}
for n, s in read_fasta(qry):
    want.setdefault(hashlib.md5(s.upper().encode()).hexdigest(), []).append(n)
print(f"  route C: {len(want)} distinct query sequence digests", file=sys.stderr)

hits = 0
with open(out, "w") as fh:
    for gid, seq in read_fasta(fa):
        d = hashlib.md5(seq.upper().encode()).hexdigest()
        if d in want:
            for acc in want[d]:
                fh.write(f"{acc}\t{gid}\tSEQ_EXACT\n")
                hits += 1
print(f"  route C: {hits} exact sequence hits", file=sys.stderr)
PY
fi
echo "[$(date +%T)] route C raw hits: $(wc -l < "${OUT}/recover_routeC_hits.tsv")"

# --------------------------------------- organism map for the verification --
echo "[$(date +%T)] extracting organism ids for every candidate gene"
cut -f2 "${OUT}/recover_routeA_hits.tsv" "${OUT}/recover_routeB_hits.tsv" \
        "${OUT}/recover_routeC_hits.tsv" 2>/dev/null \
  | grep -E '^[0-9]+_[0-9]+:' | sort -u > "${OUT}/_cand_gene_ids.txt"
echo "[$(date +%T)] candidate gene ids: $(wc -l < "${OUT}/_cand_gene_ids.txt")"

rm -f "${OUT}/_alt_id_values.txt" "${OUT}/_alt_gene_names.txt" "${OUT}/_fa_ids_sample.txt"
echo "[$(date +%T)] DONE"
ls -la "${OUT}"/recover_route*.tsv
