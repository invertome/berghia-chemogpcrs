#!/bin/bash
# diag_nautilus_es.sh — quick quality read on the ES (evidence-free) BRAKER4
# annotation of Nautilus macromphalus (bead -u25), to decide whether ES is good
# enough to harvest as-is or whether the evidence-poor lineage needs an
# enriched-evidence re-annotation.
#
#   1. proteome BUSCO (metazoa_odb12, proteins mode) on braker.longest.aa
#      -> compare proteome completeness to the 98.4% genome-level BUSCO ceiling.
#   2. GPCR proxy: hmmsearch the TIAMMAT Mollusca GPCR HMM set vs the proteome
#      -> are our chemoreceptor targets plausibly annotated?
#
# Read-only QC: writes only into a new diag_es_assessment/ subdir; touches no
# production proteome. busco_env = BUSCO 5.8.0 (genome run used container 6.0.0;
# same metazoa_odb12 lineage + cutoffs, so the % is directly comparable).
#
# Submit:  sbatch scripts/unity/diag_nautilus_es.sh
#
#SBATCH --job-name=naut_es_diag
#SBATCH --partition=cpu
#SBATCH --qos=long
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/naut_es_diag-%j.out
#SBATCH --error=logs/naut_es_diag-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate busco_env
set -u

REPO=/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05
SP=34576_Nautilus_macromphalus
ESDIR="${REPO}/species_tree_data/braker4_run/array_runs/nautilus_es/output/${SP}"
PROT="${ESDIR}/braker.longest.aa"
LINDL="${REPO}/external/braker4/shared_data/busco_downloads"
TIAMMAT="${REPO}/references/tiammat_mollusca_gpcr.hmm"
OUT="${ESDIR}/diag_es_assessment"
CPUS="${SLURM_CPUS_PER_TASK:-8}"

[ -s "$PROT" ]    || { echo "[diag] ERROR: missing $PROT"; exit 1; }
[ -s "$TIAMMAT" ] || { echo "[diag] ERROR: missing $TIAMMAT"; exit 1; }
mkdir -p "${OUT}"; cd "${OUT}"
echo "[diag] proteome: ${PROT} ($(grep -c '^>' "$PROT") genes); busco $(busco --version 2>/dev/null)"

# --- 1. proteome-level BUSCO ------------------------------------------------
busco -i "${PROT}" -o proteins_busco -m proteins -l metazoa_odb12 \
      --offline --download_path "${LINDL}" -c "${CPUS}" -f --out_path "${OUT}" \
    || echo "[diag] WARN: busco exited nonzero (see log)"
echo "===== PROTEOME BUSCO (metazoa_odb12) ====="
grep -E "C:|Complete|Single|Duplicated|Fragmented|Missing|n:" \
     "${OUT}/proteins_busco/short_summary."*.txt 2>/dev/null || echo "(no summary produced)"

# --- 2. GPCR presence proxy (TIAMMAT Mollusca GPCR HMMs) --------------------
hmmsearch --cpu "${CPUS}" -E 1e-5 --tblout "${OUT}/tiammat_gpcr.tblout" \
          "${TIAMMAT}" "${PROT}" > "${OUT}/tiammat_gpcr.hmmsearch.log" 2>&1 \
    || echo "[diag] WARN: hmmsearch exited nonzero (see log)"
n_gpcr=$(awk '!/^#/{print $1}' "${OUT}/tiammat_gpcr.tblout" 2>/dev/null | sort -u | wc -l)
echo "===== GPCR PROXY (TIAMMAT, E<1e-5) ====="
echo "[diag] ${n_gpcr} unique proteins with >=1 7TM/chemoreceptor-family hit"
echo "[diag] DONE -> ${OUT}"
