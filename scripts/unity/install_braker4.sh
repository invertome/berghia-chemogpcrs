#!/bin/bash
# install_braker4.sh — one-shot Unity setup for Phase 1f BRAKER4 (bead -p49).
#
# Idempotent: re-runs skip steps whose artifacts already exist. Designed
# to run on a Unity compute node via srun (NOT the login node). The clone
# + container pull + OrthoDB download all need network; Unity compute
# nodes have outbound HTTPS.
#
# Five steps:
#   1. Clone Gaius-Augustus/BRAKER4 to external/braker4/ (gitignored)
#   2. Create conda env from scripts/unity/envs/braker4_snakemake.yml
#   3. Pull teambraker/braker3:v3.0.10 into Apptainer image cache
#   4. Download OrthoDB metazoa partition (~2-3 GB) and gunzip
#   5. Copy AUGUSTUS config dir from container to writable scratch
#
# Submit via:
#   srun --partition=cpu --qos=long --time=2:00:00 \
#        --cpus-per-task=4 --mem=16G --account=pi_pkatz_umass_edu \
#        bash scripts/unity/install_braker4.sh

set -eo pipefail

# Conda activation BEFORE set -u (per the global trap rule)
source ~/.miniconda3/etc/profile.d/conda.sh
set +u
set -u 2>/dev/null || true

# Resolve script + project paths (assumes cwd is project root on Unity)
PROJECT_DIR="$(pwd)"
SCRIPT_DIR="$PROJECT_DIR/scripts/unity"

# Scratch root is the workspace where genomes / orthodb / .singularity_cache live.
# Allow override via env var; default to the standard workspace layout.
: "${SCRATCH_ROOT:=$PROJECT_DIR/species_tree_data}"
mkdir -p "$SCRATCH_ROOT"

EXTERNAL_DIR="$PROJECT_DIR/external/braker4"
ORTHODB_DIR="$SCRATCH_ROOT/orthodb"
ORTHODB_FA="$ORTHODB_DIR/odb12_metazoa.fa"
BRAKER4_RUN_DIR="$SCRATCH_ROOT/braker4_run"
SIF_CACHE="$BRAKER4_RUN_DIR/.singularity_cache"
AUGUSTUS_CONFIG="$BRAKER4_RUN_DIR/augustus_config"

echo "=== Phase 1f BRAKER4 install starting at $(date -Iseconds) ==="
echo "  PROJECT_DIR:        $PROJECT_DIR"
echo "  SCRATCH_ROOT:       $SCRATCH_ROOT"
echo "  EXTERNAL_DIR:       $EXTERNAL_DIR"
echo "  ORTHODB_DIR:        $ORTHODB_DIR"
echo "  BRAKER4_RUN_DIR:    $BRAKER4_RUN_DIR"

# -------------------- Step 1: clone BRAKER4 upstream --------------------

if [ -d "$EXTERNAL_DIR/.git" ]; then
    echo "[step 1] $EXTERNAL_DIR exists; skipping clone."
else
    echo "[step 1] Cloning Gaius-Augustus/BRAKER4 into $EXTERNAL_DIR ..."
    mkdir -p "$(dirname "$EXTERNAL_DIR")"
    git clone --depth 1 https://github.com/Gaius-Augustus/BRAKER4.git "$EXTERNAL_DIR"
fi

# -------------------- Step 2: conda env --------------------

ENV_YML="$SCRIPT_DIR/envs/braker4_snakemake.yml"
if ! command -v mamba &>/dev/null; then
    echo "[step 2] mamba not found in base; installing once via conda..."
    conda install -y -n base -c conda-forge mamba
fi

if conda env list 2>/dev/null | awk '{print $1}' | grep -qx braker4_runner; then
    echo "[step 2] conda env 'braker4_runner' already exists; skipping create."
else
    echo "[step 2] Creating conda env 'braker4_runner' from $ENV_YML ..."
    mamba env create -f "$ENV_YML"
fi

# Activate for subsequent steps (mostly to make `snakemake` available
# for the AUGUSTUS-config-copy step if needed).
set +u
conda activate braker4_runner
set -u 2>/dev/null || true

# -------------------- Step 3: pull Apptainer container --------------------

mkdir -p "$SIF_CACHE"
BRAKER3_SIF="$SIF_CACHE/teambraker_braker3_v3.0.10.sif"
if [ -f "$BRAKER3_SIF" ]; then
    echo "[step 3] Container already cached at $BRAKER3_SIF; skipping pull."
else
    echo "[step 3] Pulling teambraker/braker3:v3.0.10 to $BRAKER3_SIF ..."
    apptainer pull "$BRAKER3_SIF" docker://teambraker/braker3:v3.0.10
fi

# -------------------- Step 4: download OrthoDB metazoa --------------------

mkdir -p "$ORTHODB_DIR"
if [ -f "$ORTHODB_FA" ] && [ -s "$ORTHODB_FA" ]; then
    echo "[step 4] OrthoDB metazoa already at $ORTHODB_FA ($(du -h "$ORTHODB_FA" | cut -f1)); skipping."
else
    # The Bremen partitioned OrthoDB12 download is the canonical EP source per
    # the BRAKER4 docs. Fall back to Greifswald mirror if the primary is down.
    echo "[step 4] Downloading OrthoDB metazoa (~2-3 GB) ..."
    URL_PRIMARY="https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Metazoa.fa.gz"
    cd "$ORTHODB_DIR"
    if ! wget -q --show-progress -O Metazoa.fa.gz "$URL_PRIMARY"; then
        echo "[step 4] Primary URL failed; halting (no fallback configured)."
        exit 4
    fi
    gunzip -f Metazoa.fa.gz
    mv Metazoa.fa "$(basename "$ORTHODB_FA")"
    cd "$PROJECT_DIR"
    echo "[step 4] OrthoDB metazoa unpacked: $(du -h "$ORTHODB_FA" | cut -f1)"
fi

# -------------------- Step 5: writable AUGUSTUS config --------------------

if [ -d "$AUGUSTUS_CONFIG" ] && [ -d "$AUGUSTUS_CONFIG/species" ]; then
    echo "[step 5] $AUGUSTUS_CONFIG already populated; skipping copy."
else
    echo "[step 5] Copying AUGUSTUS config from container to $AUGUSTUS_CONFIG ..."
    mkdir -p "$AUGUSTUS_CONFIG"
    apptainer exec "$BRAKER3_SIF" bash -lc \
        'cp -r ${AUGUSTUS_CONFIG_PATH:-/opt/Augustus/config}/* '"$AUGUSTUS_CONFIG"'/'
fi

# -------------------- Summary --------------------

echo
echo "=== Install complete at $(date -Iseconds) ==="
echo "  BRAKER4 repo:        $EXTERNAL_DIR ($(cd "$EXTERNAL_DIR" && git rev-parse --short HEAD 2>/dev/null || echo '?'))"
echo "  Conda env:           braker4_runner ($(snakemake --version 2>/dev/null || echo 'snakemake not on PATH'))"
echo "  Container:           $BRAKER3_SIF ($(du -h "$BRAKER3_SIF" | cut -f1))"
echo "  OrthoDB metazoa:     $ORTHODB_FA ($(du -h "$ORTHODB_FA" | cut -f1))"
echo "  AUGUSTUS config:     $AUGUSTUS_CONFIG ($(find "$AUGUSTUS_CONFIG/species" -maxdepth 1 -type d 2>/dev/null | wc -l) species dirs)"
echo
echo "Next: download genomes via scripts/download_species_tree_phase1f_genomes.py"
