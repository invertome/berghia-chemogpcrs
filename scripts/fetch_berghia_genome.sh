#!/bin/bash
# fetch_berghia_genome.sh — Download Berghia stephanieae genome + RefSeq
# annotations from NCBI and link them at the canonical pipeline paths.
#
# Bead -4nu. The annotation is on RefSeq (GCF_034508935.2), not GenBank
# (GCA_034508935.4), released 2026-01-29 by the NCBI annotation pipeline.
# Source: Goodheart et al. 2024 BMC Biology 22:9.
#
# Usage:
#   bash scripts/fetch_berghia_genome.sh
#
# Requires: NCBI Datasets CLI (`datasets`) — typically `pip install ncbi-datasets-cli`
#           or available via conda-forge (already in environment.yml prerequisites).
#
# Idempotent: if the canonical files already exist, exits with success.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Allow override via environment
ACCESSION="${BERGHIA_GENOME_ACCESSION:-GCF_034508935.2}"
VERSION="${BERGHIA_GENOME_VERSION:-UCSD_Bste_1.3}"
GENOME_DIR="${GENOME_DIR:-$BASE_DIR/genomes}"

# Canonical paths use the real NCBI taxid + genus + species
# (e.g. "1287507_berghia_stephanieae"). Must match BERGHIA_FILE_PREFIX in
# config.sh. We respect any override via the environment.
PREFIX="${BERGHIA_FILE_PREFIX:-${BERGHIA_TAXID:-1287507}_${BERGHIA_GENUS:-berghia}_${BERGHIA_SPECIES:-stephanieae}}"
GENOME_FA="$GENOME_DIR/${PREFIX}.fasta"
GENOME_GFF="$GENOME_DIR/${PREFIX}.gff3"
GENOME_PROT="$GENOME_DIR/${PREFIX}.proteins.fa"
GENOME_CDS="$GENOME_DIR/${PREFIX}.cds.fna"

# Legacy aliases (point at the canonical files; kept so any pre-2026-05 code
# that hardcoded "taxid_berghia_berghia.*" still works during transition).
LEGACY_FA="$GENOME_DIR/taxid_berghia_berghia.fasta"
LEGACY_GFF="$GENOME_DIR/taxid_berghia_berghia.gff3"
LEGACY_PROT="$GENOME_DIR/taxid_berghia_berghia.proteins.fa"
LEGACY_CDS="$GENOME_DIR/taxid_berghia_berghia.cds.fna"

mkdir -p "$GENOME_DIR"

if [ -e "$GENOME_FA" ] && [ -e "$GENOME_GFF" ] && [ -e "$GENOME_PROT" ] && [ -e "$GENOME_CDS" ]; then
    echo "[fetch_berghia_genome] All canonical genome files already present." >&2
    ls -lh "$GENOME_FA" "$GENOME_GFF" "$GENOME_PROT" "$GENOME_CDS" >&2
    exit 0
fi

if ! command -v datasets &>/dev/null; then
    echo "ERROR: NCBI Datasets CLI ('datasets') not found in PATH." >&2
    echo "       Install via 'conda install -c conda-forge ncbi-datasets-cli'" >&2
    echo "       or 'pip install ncbi-datasets-cli'." >&2
    exit 2
fi

ZIP_PATH="$GENOME_DIR/berghia_refseq_${ACCESSION}.zip"
if [ ! -e "$ZIP_PATH" ]; then
    echo "[fetch_berghia_genome] Downloading $ACCESSION (~370MB) ..." >&2
    datasets download genome accession "$ACCESSION" \
        --include genome,gff3,protein,cds \
        --filename "$ZIP_PATH"
fi

echo "[fetch_berghia_genome] Unzipping ..." >&2
(cd "$GENOME_DIR" && unzip -q -o "$ZIP_PATH")

DATA_DIR="$GENOME_DIR/ncbi_dataset/data/$ACCESSION"
if [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: Expected $DATA_DIR after unzip but it is missing." >&2
    exit 3
fi

# Pattern restricted to GCA/GCF accession prefix to avoid matching
# 'cds_from_genomic.fna'.
SRC_FA=$(ls "$DATA_DIR"/GC[AF]_*_genomic.fna 2>/dev/null | head -1)
SRC_GFF="$DATA_DIR/genomic.gff"
SRC_PROT="$DATA_DIR/protein.faa"
SRC_CDS="$DATA_DIR/cds_from_genomic.fna"

for src in "$SRC_FA" "$SRC_GFF" "$SRC_PROT" "$SRC_CDS"; do
    if [ ! -e "$src" ]; then
        echo "ERROR: Missing expected file from RefSeq archive: $src" >&2
        exit 4
    fi
done

ln -sf "$SRC_FA" "$GENOME_FA"
ln -sf "$SRC_GFF" "$GENOME_GFF"
ln -sf "$SRC_PROT" "$GENOME_PROT"
ln -sf "$SRC_CDS" "$GENOME_CDS"

# Legacy aliases for backward compat with pre-2026-05 hardcoded paths.
ln -sf "$GENOME_FA"   "$LEGACY_FA"
ln -sf "$GENOME_GFF"  "$LEGACY_GFF"
ln -sf "$GENOME_PROT" "$LEGACY_PROT"
ln -sf "$GENOME_CDS"  "$LEGACY_CDS"

echo "[fetch_berghia_genome] OK." >&2
echo "  Genome:   $GENOME_FA -> $(readlink "$GENOME_FA")" >&2
echo "  GFF3:     $GENOME_GFF -> $(readlink "$GENOME_GFF")" >&2
echo "  Proteins: $GENOME_PROT -> $(readlink "$GENOME_PROT")" >&2
echo "  CDS:      $GENOME_CDS -> $(readlink "$GENOME_CDS")" >&2
echo "  Legacy aliases at $LEGACY_FA etc. for backward compat." >&2

N_SCAFF=$(awk '/^>/ {n++} END {print n}' "$GENOME_FA")
N_GENES=$(grep -c $'\tgene\t' "$GENOME_GFF")
N_PROT=$(grep -c '^>' "$GENOME_PROT")
N_CDS=$(grep -c '^>' "$GENOME_CDS")
echo "  Stats: ${N_SCAFF} scaffolds, ${N_GENES} genes, ${N_PROT} proteins, ${N_CDS} CDS" >&2
