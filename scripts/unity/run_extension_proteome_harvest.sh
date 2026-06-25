#!/bin/bash
# run_extension_proteome_harvest.sh — harvest Phase-1d extension proteomes.
#
# Reinstates the Phase-1a download step for the extension set: selects
# NCBI-annotated extension species, downloads their protein FASTAs via
# NCBI Datasets CLI, appends them to proteome_manifest.tsv, then re-derives
# the extension + genome inventories and produces a consolidation report.
#
# STOPS after writing updated manifests + the consolidation report.
# NO git commit, NO Phase-1f release. The updated manifests are for human
# review before the next pipeline stage.
#
# Bead: berghia-chemogpcrs-w2x.
#
# Usage:
#   sbatch scripts/unity/run_extension_proteome_harvest.sh
#
# Optional env overrides (--export=ALL,VAR=val):
#   EXTENSION_TSV   path to extension_inventory.tsv (default: auto-derived from REPO_ROOT)
#   PROTEOME_MANIFEST  path to proteome_manifest.tsv
#   NCBI_CACHE_DIR  where downloaded FASTAs land
#
#SBATCH --job-name=ext_prot_harvest
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/ext_prot_harvest-%j.out
#SBATCH --error=logs/ext_prot_harvest-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"

# ---------------------------------------------------------------------------
# Paths (all relative to REPO_ROOT; override via env if needed)
# ---------------------------------------------------------------------------
REFS="${REPO_ROOT}/references/species_tree"
EXTENSION_TSV="${EXTENSION_TSV:-${REFS}/extension_inventory.tsv}"
PROTEOME_MANIFEST="${PROTEOME_MANIFEST:-${REFS}/proteome_manifest.tsv}"
GENOME_INVENTORY="${REFS}/genome_inventory.tsv"
HARVEST_WORKDIR="${REFS}/harvest_extension_work"
NCBI_CACHE_DIR="${NCBI_CACHE_DIR:-${REPO_ROOT}/species_tree_data/ncbi_proteomes}"
DATASETS_BIN="${DATASETS_BIN:-datasets}"

DOWNLOAD_MANIFEST="${HARVEST_WORKDIR}/download_manifest.tsv"
STAGED_MANIFEST="${HARVEST_WORKDIR}/staged_manifest.tsv"
ACCESSIONS_FILE="${HARVEST_WORKDIR}/extension_accessions.txt"
DATASETS_JSONL="${HARVEST_WORKDIR}/extension_datasets_summary.jsonl"
DOWNLOAD_REPORT="${HARVEST_WORKDIR}/download_report.tsv"
CONSOLIDATE_OUTDIR="${REPO_ROOT}/species_tree_data/genome_wide_orthofinder_input"

mkdir -p "${HARVEST_WORKDIR}" "${NCBI_CACHE_DIR}" logs

echo "[$(date '+%H:%M:%S')] run_extension_proteome_harvest: starting"
echo "  REPO_ROOT         = ${REPO_ROOT}"
echo "  EXTENSION_TSV     = ${EXTENSION_TSV}"
echo "  PROTEOME_MANIFEST = ${PROTEOME_MANIFEST}"
echo "  NCBI_CACHE_DIR    = ${NCBI_CACHE_DIR}"

# ---------------------------------------------------------------------------
# Step (a): extract extension accessions for NCBI datasets summary
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Step (a): extracting extension accessions"
[ -f "${EXTENSION_TSV}" ] || { echo "ERROR: EXTENSION_TSV not found: ${EXTENSION_TSV}" >&2; exit 1; }

awk -F'\t' 'NR>1 && $6 != "" {print $6}' "${EXTENSION_TSV}" > "${ACCESSIONS_FILE}"
N_ACC=$(wc -l < "${ACCESSIONS_FILE}")
echo "  ${N_ACC} accessions extracted to ${ACCESSIONS_FILE}"

# ---------------------------------------------------------------------------
# Step (b): NCBI datasets summary genome accession
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Step (b): querying NCBI datasets summary"
"${DATASETS_BIN}" summary genome accession \
    --inputfile "${ACCESSIONS_FILE}" \
    --as-json-lines \
    > "${DATASETS_JSONL}"
echo "  datasets summary written to ${DATASETS_JSONL}"

# ---------------------------------------------------------------------------
# Step (c): harvest select — pick annotated, write download + staged manifests
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Step (c): harvest select"
python3 "${REPO_ROOT}/scripts/harvest_extension_proteomes.py" select \
    --extension-tsv "${EXTENSION_TSV}" \
    --datasets-jsonl "${DATASETS_JSONL}" \
    --download-manifest-out "${DOWNLOAD_MANIFEST}" \
    --staged-out "${STAGED_MANIFEST}"

N_SEL=$(awk 'NR>1' "${DOWNLOAD_MANIFEST}" | wc -l)
echo "  ${N_SEL} annotated extension species selected"
if [ "${N_SEL}" -eq 0 ]; then
    echo "[$(date '+%H:%M:%S')] WARNING: no annotated extension species found; nothing to download."
    echo "  Check that EXTENSION_TSV has annotated entries with matching NCBI records."
fi

# ---------------------------------------------------------------------------
# Step (d): download proteomes via download_species_tree_phase1a.py
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Step (d): downloading proteomes"
python3 "${REPO_ROOT}/scripts/download_species_tree_phase1a.py" \
    --manifest "${DOWNLOAD_MANIFEST}" \
    --cache-dir "${NCBI_CACHE_DIR}" \
    --report "${DOWNLOAD_REPORT}" \
    --datasets-bin "${DATASETS_BIN}" || {
    # exit 2 = some downloads failed (not fatal; append filters by status)
    echo "  WARNING: download step returned non-zero (some species may have failed; see ${DOWNLOAD_REPORT})" >&2
    true
}
echo "  download report written to ${DOWNLOAD_REPORT}"

# ---------------------------------------------------------------------------
# Step (e): harvest append — add successfully-downloaded species to proteome_manifest
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Step (e): harvest append"
python3 "${REPO_ROOT}/scripts/harvest_extension_proteomes.py" append \
    --staged "${STAGED_MANIFEST}" \
    --download-results "${DOWNLOAD_REPORT}" \
    --proteome-manifest "${PROTEOME_MANIFEST}"

# ---------------------------------------------------------------------------
# Step (f): re-derive inventories + consolidation report
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Step (f): re-derive extension inventory"
python3 "${REPO_ROOT}/scripts/build_species_tree_phase1d_extension_inventory.py" \
    --phase1a-manifest "${PROTEOME_MANIFEST}" \
    --out "${REFS}/extension_inventory_rederived.tsv"

echo "[$(date '+%H:%M:%S')] Step (f): re-derive unified genome inventory"
python3 "${REPO_ROOT}/scripts/build_genome_inventory.py" \
    --manifest "${GENOME_INVENTORY}"

echo "[$(date '+%H:%M:%S')] Step (f): consolidation report"
mkdir -p "${CONSOLIDATE_OUTDIR}"
python3 "${REPO_ROOT}/scripts/consolidate_proteomes_for_genome_wide_og.py" \
    --out-dir "${CONSOLIDATE_OUTDIR}" \
    --phase1a-manifest "${PROTEOME_MANIFEST}" \
    --manifest "${GENOME_INVENTORY}"

# ---------------------------------------------------------------------------
# Summary — stop here; no git commit, no Phase-1f release
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] DONE. Updated manifests + consolidation report written."
echo "  proteome_manifest : ${PROTEOME_MANIFEST}"
echo "  genome_inventory  : ${GENOME_INVENTORY}"
echo "  consolidation dir : ${CONSOLIDATE_OUTDIR}"
echo "  download report   : ${DOWNLOAD_REPORT}"
echo ""
echo "  Review the above before releasing Phase-1f or committing changes."
