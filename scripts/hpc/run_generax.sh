#!/bin/bash
# run_generax.sh — GeneRax replaces Notung as default reconciliation.
#
# Bead -30g. GeneRax (Morel 2020 MBE 37:2763) estimates joint phylogenetic +
# reconciliation likelihoods. Notung is parsimony-style and inferior
# (Morel 2023 GBE 15:evad134). HPC-only — runtime ~hours per family.

#SBATCH --job-name=generax
#SBATCH --output=logs/generax_%j.out
#SBATCH --error=logs/generax_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$PROJECT_DIR/config.sh"
source "$PROJECT_DIR/functions.sh"

GENERAX_BIN="${GENERAX:-generax}"
SPECIES_TREE="${SPECIES_TREE:-${RESULTS_DIR}/busco/busco_species_tree.tre}"
GENE_FAMILIES="${1:?GeneRax families.txt config}"
OUTPUT_DIR="${RESULTS_DIR}/reconciliation/generax"
mkdir -p "$OUTPUT_DIR"

if ! command -v "$GENERAX_BIN" &>/dev/null; then
    echo "ERROR: GeneRax binary '$GENERAX_BIN' not found. " \
         "See docs/installation_hpc.md for build instructions." >&2
    exit 2
fi
if [ ! -f "$SPECIES_TREE" ]; then
    echo "ERROR: species tree not found at $SPECIES_TREE" >&2
    exit 1
fi

log "GeneRax: $GENE_FAMILIES against $SPECIES_TREE -> $OUTPUT_DIR"

mpiexec -n "${SLURM_CPUS_PER_TASK:-8}" "$GENERAX_BIN" \
    --species-tree "$SPECIES_TREE" \
    --families "$GENE_FAMILIES" \
    --rec-model UndatedDTL \
    --strategy SPR \
    --prefix "$OUTPUT_DIR" \
    > "$OUTPUT_DIR/generax.log" 2>&1

log "GeneRax complete -> $OUTPUT_DIR"
log "Update lse_refine.py to consume GeneRax NHX-tagged trees from " \
    "$OUTPUT_DIR/results/*/geneTree.newick for cleaner LSE classification."
