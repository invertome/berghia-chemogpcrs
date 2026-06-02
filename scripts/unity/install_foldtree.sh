#!/bin/bash
#SBATCH --job-name=install_foldtree
#SBATCH --partition=cpu
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/install_foldtree-%j.out
#SBATCH --error=logs/install_foldtree-%j.err
# install_foldtree.sh — install FoldTree on Unity (bead -68w).
#
# FoldTree (Dessimoz Lab) builds a phylogenetic tree from protein structures via
# Foldseek structural distances. Repo (verified): https://github.com/DessimozLab/fold_tree
# It is a Snakemake pipeline exposing a `foldtree` click CLI, so it gets its OWN
# conda env (python 3.10 + snakemake-minimal) — NOT berghia-gpcr: the tool needs
# its own interpreter and spawns per-rule conda envs at runtime. A fresh env
# leaves berghia-gpcr (and its tmbed/transformers pin) untouched.
#
# Real usage (per upstream README):
#   foldtree --cores 4 --folder <dir> --custom-structs -p
#   input layout: <dir>/structs/*.pdb   ->   outputs Foldtree/LDDT/TM trees
# NOTE: stage 08 currently invokes ${FOLDTREE} --input_dir/--output/--method,
# which does NOT match this CLI. Reconciling stage 08 (FOLDTREE path + a thin
# CLI adapter, or a stage-08 update) is a separate follow-up.
#
# Idempotent: re-clone is a git pull; existing env/package installs are no-ops.
# Submit: sbatch scripts/unity/install_foldtree.sh

set -eo pipefail
mkdir -p logs

WORKDIR="/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05"
TOOLS_DIR="${WORKDIR}/tools"
REPO_DIR="${TOOLS_DIR}/fold_tree"
ENV_NAME="${FOLDTREE_ENV:-foldtree}"
mkdir -p "${TOOLS_DIR}"

# shellcheck disable=SC1091
source ~/.miniconda3/etc/profile.d/conda.sh
SOLVER=mamba
command -v "$SOLVER" >/dev/null 2>&1 || SOLVER=conda
echo "[foldtree] solver: $SOLVER"

# 1. Clone or update the repo.
if [ -d "${REPO_DIR}/.git" ]; then
    echo "[foldtree] repo present at ${REPO_DIR}; updating"
    git -C "${REPO_DIR}" pull --ff-only || echo "[foldtree] WARN: pull skipped"
else
    echo "[foldtree] cloning DessimozLab/fold_tree"
    git clone https://github.com/DessimozLab/fold_tree.git "${REPO_DIR}"
fi

# 2. Dedicated env (idempotent).
if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
    echo "[foldtree] env '${ENV_NAME}' already exists"
else
    echo "[foldtree] creating env '${ENV_NAME}' (python=3.10)"
    "$SOLVER" create -y -n "${ENV_NAME}" python=3.10
fi
conda activate "${ENV_NAME}"
echo "[foldtree] env active: $(which python)"

# 3. Pipeline deps + the package itself.
# gemmi: used by scripts/run_foldtree.sh to convert AF3 mmCIF -> .pdb, since
# fold_tree's structs2fasta.py globs structs/*.pdb only.
"$SOLVER" install -y -c bioconda -c conda-forge click snakemake-minimal gemmi
( cd "${REPO_DIR}" && python -m pip install . --no-deps --no-build-isolation --no-cache-dir )

# 4. Verify.
FOLDTREE_BIN="$(command -v foldtree || true)"
echo "[foldtree] binary: ${FOLDTREE_BIN:-MISSING}"
if [ -z "${FOLDTREE_BIN}" ]; then
    echo "[foldtree] FATAL: foldtree not on PATH after install." >&2
    exit 1
fi
foldtree --help 2>&1 | head -20 || true

echo "===================================================================="
echo "[foldtree] installed OK."
echo "  FOLDTREE path (for config.sh): ${FOLDTREE_BIN}"
echo "  Follow-up: stage 08's --input_dir/--output/--method call needs a thin"
echo "  adapter (real CLI: --folder <dir w/ structs/> --custom-structs -p)."
echo "===================================================================="
