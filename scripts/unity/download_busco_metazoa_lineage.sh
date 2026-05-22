#!/bin/bash
#SBATCH --job-name=busco_lineage_dl
#SBATCH --partition=cpu
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/busco_lineage_dl-%j.out
#SBATCH --error=logs/busco_lineage_dl-%j.err
#SBATCH --account=pi_pkatz_umass_edu

# One-time download of the metazoa_odb12 BUSCO lineage database for
# BRAKER4's compleasm rule (bead -p49). The `run_compleasm` rule expects
# a flat layout at:
#   external/braker4/shared_data/busco_downloads/lineages/metazoa_odb12/
# When the lineage is missing, compleasm tries to download it mid-rule
# (busco-data.ezlab.org), which is fragile under concurrent sbatch jobs.
# Prefetching it once is the recommended pattern.
#
# Uses compleasm.py from inside the BRAKER3 container (where it's at
# /opt/compleasm_kit/compleasm.py — not on $PATH, so we invoke explicitly).
#
# Idempotent: if metazoa_odb12 already exists at the target, exits 0.

set -eo pipefail
mkdir -p logs

cd /scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05

LINEAGE_DIR="external/braker4/shared_data/busco_downloads/lineages"
SIF="species_tree_data/braker4_run/.singularity_cache/teambraker_braker3_v3.0.10.sif"

if [ -d "${LINEAGE_DIR}/metazoa_odb12" ] && [ "$(ls -A "${LINEAGE_DIR}/metazoa_odb12" 2>/dev/null)" ]; then
    echo "metazoa_odb12 already present at ${LINEAGE_DIR}/metazoa_odb12 — skipping"
    ls -lh "${LINEAGE_DIR}/metazoa_odb12" | head -10
    exit 0
fi

if [ ! -f "$SIF" ]; then
    echo "ERROR: BRAKER3 container not found at $SIF" >&2
    exit 2
fi

mkdir -p "$LINEAGE_DIR"

module load apptainer/latest 2>/dev/null || true

echo "Downloading metazoa_odb12 lineage into ${LINEAGE_DIR}/ ..."
apptainer exec -B "$(realpath "$LINEAGE_DIR"):$(realpath "$LINEAGE_DIR")" "$SIF" \
    python3 /opt/compleasm_kit/compleasm.py download metazoa \
        -L "$(realpath "$LINEAGE_DIR")" \
        --odb odb12

echo
echo "=== Verifying download ==="
if [ -d "${LINEAGE_DIR}/metazoa_odb12" ]; then
    echo "metazoa_odb12 contents:"
    ls -lh "${LINEAGE_DIR}/metazoa_odb12" | head -10
    echo "Total size:"
    du -sh "${LINEAGE_DIR}/metazoa_odb12"
else
    echo "ERROR: metazoa_odb12 directory not created" >&2
    echo "Library contents:"
    ls -la "$LINEAGE_DIR" >&2 || true
    exit 3
fi
