#!/bin/bash
# alphafast_setup_databases.sh — download AlphaFast's pre-built AF3 MMseqs2 DBs
# (HuggingFace, protein-only) for the stage-08 AlphaFast backend (bead q0o.1).
#
# WHY protein-only: we fold protein monomers (chemoreceptor GPCRs), so the RNA
# MMseqs2 DBs (~1.48 TB of the 2.22 TB HF repo) are never used. --protein-only
# fetches just mmseqs/* (padded protein DBs, ~584 GB) + mmcif_files (~57 GB).
#
# WHY a local wrapper (not AlphaFast's scripts/setup_databases.sbatch): that one
# hard-codes Duke's `--partition=common` and activates no conda env. We override
# to Unity's `cpu` partition and activate the `alphafast_dl` env (hf CLI + zstd).
#
# Submit FROM the workspace dir so the relative logs/ path resolves there:
#   cd /scratch3/workspace/$USER-jorge && sbatch <repo>/scripts/unity/alphafast_setup_databases.sh
#
#SBATCH --job-name=alphafast_db
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/alphafast_db-%j.out
#SBATCH --error=logs/alphafast_db-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate alphafast_dl          # hf CLI + zstd; activate BEFORE set -u
set -u

AFR="${ALPHAFAST_REPO:-/scratch3/workspace/$USER-jorge/alphafast}"
DBDIR="${ALPHAFAST_DB_DIR:-/scratch3/workspace/$USER-jorge/alphafast_db}"
mkdir -p "$DBDIR"

echo "[$(date)] AlphaFast DB setup (protein-only, HuggingFace pre-built)"
echo "  repo:   $AFR"
echo "  db_dir: $DBDIR"
command -v hf zstd tar

"$AFR/scripts/setup_databases.sh" "$DBDIR" --protein-only

echo "[$(date)] DONE. db_dir contents:"
ls -la "$DBDIR" "$DBDIR/mmseqs" 2>/dev/null | head -40
echo "[$(date)] db_dir usage:"; du -sh "$DBDIR" 2>/dev/null || true
