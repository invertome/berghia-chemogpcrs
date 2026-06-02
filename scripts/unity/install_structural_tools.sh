#!/bin/bash
#SBATCH --job-name=install_struct_tools
#SBATCH --partition=cpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/install_struct_tools-%j.out
#SBATCH --error=logs/install_struct_tools-%j.err
# install_structural_tools.sh — install foldseek + TMalign into the
# berghia-gpcr conda env for Gate 7 (stage 08) structural-phylogeny work
# (bead -05o).
#
# Background (Gate 7 audit 2026-05-18):
#   - Unity already has AlphaFold3 v3.0.1 as a system Apptainer module
#     (loaded with `module load alphafold3/latest`) and its 2 TB of
#     databases at /datasets/bio/alphafold3 — no local install needed
#     for AlphaFold itself; that wiring is a separate task.
#   - foldseek is installed in conda BASE, not berghia-gpcr — so the
#     stage 08 sbatch (which activates berghia-gpcr) cannot find it.
#   - TMalign is not installed anywhere.
#   - config.sh references ${FOLDTREE}=foldtree (Reading Lab repo,
#     separate clone — out of scope for this task).
#
# This script makes foldseek + TMalign available inside berghia-gpcr
# via bioconda. Mamba is the preferred solver on Unity for dev-side env ops.
#
# Idempotent: re-running is a no-op when both binaries are present.
#
# Submit: sbatch scripts/unity/install_structural_tools.sh

set -eo pipefail
mkdir -p logs

ENV_NAME="${ENV_NAME:-berghia-gpcr}"

# shellcheck disable=SC1091
source ~/.miniconda3/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

echo "[install-struct] env: $ENV_NAME ($(which python))"

# Probe each binary inside the active env (not base PATH).
need_install=()
if ! command -v foldseek >/dev/null 2>&1; then
    need_install+=("foldseek")
else
    echo "[install-struct] foldseek: $(which foldseek)"
fi
if ! command -v TMalign >/dev/null 2>&1; then
    need_install+=("tmalign")
else
    echo "[install-struct] TMalign: $(which TMalign)"
fi

if [ "${#need_install[@]}" -eq 0 ]; then
    echo "[install-struct] All tools already installed — nothing to do."
    exit 0
fi

# bioconda's tmalign package provides /bin/TMalign; foldseek provides /bin/foldseek.
# Use mamba (faster solver) per user preference for dev-side env ops.
SOLVER=mamba
command -v "$SOLVER" >/dev/null 2>&1 || SOLVER=conda

echo "[install-struct] Installing ${need_install[*]} via $SOLVER + bioconda..."
# --freeze-installed: never mutate already-installed packages (protects the
# tmbed/transformers<5 pin in berghia-gpcr); a conflicting tool fails loudly.
"$SOLVER" install -y --freeze-installed -n "$ENV_NAME" -c bioconda -c conda-forge "${need_install[@]}"

# Reactivate so PATH picks up the new binaries deterministically.
conda deactivate
conda activate "$ENV_NAME"

# Verify.
FOLDSEEK_BIN="$(command -v foldseek || true)"
TMALIGN_BIN="$(command -v TMalign || true)"
echo "[install-struct] foldseek: ${FOLDSEEK_BIN:-MISSING}"
echo "[install-struct] TMalign:  ${TMALIGN_BIN:-MISSING}"

if [ -z "$FOLDSEEK_BIN" ] || [ -z "$TMALIGN_BIN" ]; then
    echo "[install-struct] FATAL: one or more binaries missing after install." >&2
    exit 1
fi

# Smoke test: foldseek prints version banner; TMalign prints usage on bare invocation.
echo "[install-struct] foldseek version:"
foldseek version 2>&1 | head -3 || true
echo "[install-struct] TMalign banner:"
TMalign 2>&1 | head -3 || true

echo "[install-struct] Done. Stage 08 can now use foldseek + TMalign from the berghia-gpcr env."
