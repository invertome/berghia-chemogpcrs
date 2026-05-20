#!/bin/bash
#SBATCH --job-name=install_pytest
#SBATCH --partition=cpu
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/install_pytest-%j.out
#SBATCH --error=logs/install_pytest-%j.err
# install_pytest_berghia_gpcr.sh — install pytest stack into the
# berghia-gpcr conda env so the Python unit-test suite actually runs
# under sbatch_run_tests.sh (bead -gli).
#
# Background: sbatch_run_tests.sh imports the env and calls `pytest
# tests/unit/`, but pytest isn't in the env — so every Python unit test
# is silently skipped. Only the bash unit tests run, giving false-green
# confidence. This script fixes that:
#
#   1. activate berghia-gpcr
#   2. pip install pytest pytest-cov pytest-mock
#   3. verify pytest >= 7
#   4. run the full Python suite to confirm collection works
#
# Idempotent: re-running is a no-op if pytest already installed.
#
# Submit: sbatch scripts/unity/install_pytest_berghia_gpcr.sh

set -eo pipefail
mkdir -p logs

ENV_NAME="${ENV_NAME:-berghia-gpcr}"

# shellcheck disable=SC1091
source ~/.miniconda3/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

echo "[install-pytest] env: $ENV_NAME ($(which python))"

# Skip the install if already present at a recent enough version.
if python3 -c "import pytest, sys; sys.exit(0 if pytest.__version__.split('.')[0] >= '7' else 1)" 2>/dev/null; then
    echo "[install-pytest] pytest already installed: $(python3 -c 'import pytest; print(pytest.__version__)')"
else
    echo "[install-pytest] Installing pytest pytest-cov pytest-mock..."
    pip install --quiet 'pytest>=7' pytest-cov pytest-mock
fi

PY_VER="$(python3 -c 'import pytest; print(pytest.__version__)')"
echo "[install-pytest] pytest $PY_VER ready"

# Confirm pytest can at least COLLECT the suite without import errors.
# We don't enforce all-green here; the goal is "Python suite is now actually
# runnable under sbatch_run_tests.sh". A separate run validates correctness.
cd "$(dirname "$0")/../.."
echo "[install-pytest] Collecting tests/unit/..."
if pytest --collect-only -q tests/unit/ 2>&1 | tail -10; then
    echo "[install-pytest] Collection OK."
else
    echo "[install-pytest] WARN: pytest collection had issues — see above." >&2
fi

echo "[install-pytest] Done. Re-run sbatch_run_tests.sh to execute the full suite."
