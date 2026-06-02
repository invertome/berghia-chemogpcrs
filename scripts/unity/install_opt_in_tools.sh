#!/bin/bash
#SBATCH --job-name=install_opt_in_tools
#SBATCH --partition=cpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/install_opt_in_tools-%j.out
#SBATCH --error=logs/install_opt_in_tools-%j.err
# install_opt_in_tools.sh — install opt-in / fallback pipeline tools into the
# berghia-gpcr env (bead -68w).
#
# Uses `--freeze-installed` so already-installed packages (notably the
# transformers<5 pin tmbed needs) can NEVER be mutated — a tool that would
# require changing them fails loudly and is skipped, instead of silently
# breaking the env. Safe to run while other jobs use berghia-gpcr.
#
# Tools (bioconda), each invoked only by an opt-in / non-default code path:
#   fastml   — stage 05 ASR fallback         (ASR_BACKEND=fastml; iqtree3 is default)
#   cafe     — stage 03c CAFE5 expansion      (also needs a species tree, bd-dnk)
#   mrbayes  — stage 04 USE_MRBAYES=true       (binary 'mb'; IQ-TREE is canonical)
#   mcscanx  — stage 06 SYNTENY_BACKEND=mcscanx (jcvi mcscan is default)
#
# NOT handled here:
#   foldseek / TMalign — install_structural_tools.sh (bead -05o).
#   foldtree           — DEFERRED: config only references a bare `foldtree`
#                        binary with a --input_dir/--output/--method CLI that
#                        does not match the Dessimoz fold_tree Snakemake tool;
#                        the correct repo/install must be verified before
#                        installing (Gate 7 / stage 08, not urgent).
#
# Each tool installs independently; one failure does not abort the others.
# Idempotent: present tools are skipped. End-of-run summary reports each status.
#
# Submit: sbatch scripts/unity/install_opt_in_tools.sh
# (Do NOT run concurrently with install_structural_tools.sh — two solvers
#  modifying one env at once can corrupt it; chain with --dependency=afterok.)

set -eo pipefail
mkdir -p logs
ENV_NAME="${ENV_NAME:-berghia-gpcr}"

# shellcheck disable=SC1091
source ~/.miniconda3/etc/profile.d/conda.sh
conda activate "$ENV_NAME"
echo "[install-opt-in] env: $ENV_NAME ($(which python))"

SOLVER=mamba
command -v "$SOLVER" >/dev/null 2>&1 || SOLVER=conda
echo "[install-opt-in] solver: $SOLVER"

# install_one <probe-binary> <bioconda-pkg> <label>
# Returns 0 if present-or-installed, 1 on install failure.
install_one() {
    local probe="$1" pkg="$2" label="$3"
    if command -v "$probe" >/dev/null 2>&1; then
        echo "[install-opt-in] $label already present ($(command -v "$probe")) — skip"
        return 0
    fi
    echo "[install-opt-in] installing $label ($pkg) with --freeze-installed ..."
    if "$SOLVER" install -y --freeze-installed -n "$ENV_NAME" \
            -c bioconda -c conda-forge "$pkg"; then
        echo "[install-opt-in] $label: install OK"
        return 0
    fi
    echo "[install-opt-in] WARN: $label ($pkg) install FAILED — likely a" \
         "--freeze-installed conflict; skipped (env left untouched)." >&2
    return 1
}

declare -A STATUS
install_one fastml  fastml  fastml   && STATUS[fastml]=ok  || STATUS[fastml]=fail
install_one cafe5   cafe    cafe5    && STATUS[cafe]=ok    || STATUS[cafe]=fail
install_one mb      mrbayes mrbayes  && STATUS[mrbayes]=ok || STATUS[mrbayes]=fail
install_one MCScanX mcscanx mcscanx  && STATUS[mcscanx]=ok || STATUS[mcscanx]=fail

# Reactivate so PATH picks up any new binaries.
conda deactivate
conda activate "$ENV_NAME"

echo "==================== install summary ===================="
for t in fastml cafe mrbayes mcscanx; do
    printf '  %-10s %s\n' "$t" "${STATUS[$t]:-unknown}"
done
echo "  foldseek/TMalign  -> install_structural_tools.sh (-05o)"
echo "  foldtree          -> DEFERRED (verify repo/CLI first)"
echo "========================================================="

# Best-effort version banners for what installed.
command -v fastml  >/dev/null 2>&1 && { echo "[ver] fastml:";  fastml  2>&1 | head -2 || true; }
command -v cafe5   >/dev/null 2>&1 && { echo "[ver] cafe5:";   cafe5 --version 2>&1 | head -1 || true; }
command -v mb      >/dev/null 2>&1 && echo "[ver] mb present: $(command -v mb)"
command -v MCScanX >/dev/null 2>&1 && echo "[ver] MCScanX present: $(command -v MCScanX)"

echo "[install-opt-in] done."
