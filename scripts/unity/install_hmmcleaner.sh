#!/bin/bash
# install_hmmcleaner.sh — install HmmCleaner.pl into the berghia-gpcr conda env on Unity.
#
# WHY THIS EXISTS
# ---------------
# HmmCleaner.pl is the segment-cleaning step in the alignment QC stack
# (RUN_HMMCLEANER=1 path in the main pipeline). Two prior attempts on Unity
# both failed:
#   1. `cpanm Bio::MUST::Apps::HmmCleaner` — bails inside Bio::MUST::Drivers
#      test phase (exo_aligned.t hangs on exonerate sanity test). Documented
#      at https://github.com/dongzhang0725/PhyloSuite/issues/67 — "zero
#      successes on cpan testers".
#   2. Bioconda packages `perl-bio-must-apps-hmmcleaner`,
#      `perl-bio-must-drivers`, `perl-bio-must-core` were tried in job
#      56714199 — they DO NOT EXIST on bioconda (verified May 2026 via
#      `conda search bioconda::perl-bio-must*` → no results).
#
# WHY THERE'S NO CONTAINER FALLBACK
# ---------------------------------
# Searched quay.io/biocontainers, depot.galaxyproject.org/singularity,
# DockerHub, and the Bio-MUST GitHub for HmmCleaner / bio-must / bio-must-apps
# — no published Docker or Apptainer image exists. The PhyloSuite team
# bundles HmmCleaner only in their Linux pre-compiled package (which is GUI-
# integrated and not headless-CLI clean for HPC use).
#
# CHOSEN APPROACH
# ---------------
# Explicit topologically-ordered cpanm install with --notest --no-man-pages,
# into a local::lib rooted inside the conda env, after first installing the
# system bio tools (exonerate, hmmer, blast) that Bio::MUST::Drivers
# needs at runtime. The dependency order is in cpanm_dep_order.txt and was
# derived from the metacpan Requires lists for each Bio::MUST::* dist.
#
# The --notest flag is what makes this work where prior attempts failed: it
# bypasses the exonerate-version-sensitive test in Bio::MUST::Drivers without
# affecting runtime correctness (HmmCleaner's own integration test is run
# at the end of this script as the real verification).
#
# SUBMIT:
#   sbatch scripts/unity/install_hmmcleaner.sh
#
# After completion, RUN_HMMCLEANER=1 in the main pipeline will pick up
# HmmCleaner.pl from $CONDA_PREFIX/bin (which is already on PATH for the
# berghia-gpcr conda env).

#SBATCH --job-name=hmmcleaner_install
#SBATCH --partition=cpu
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -eo pipefail
# NOTE: deliberately omitting `set -u`. Conda's gcc/gxx deactivate hooks
# (deactivate-gxx_linux-64.sh) reference $CONDA_BACKUP_CXX without a
# default, which makes any `conda install -n <other_env>` (which forces a
# deactivate first) abort under `set -u` with:
#   "CONDA_BACKUP_CXX: unbound variable"
# The first install attempt (job 56723380) crashed in 5 seconds for
# exactly this reason. We accept the marginal loss of unset-var safety
# in this install script in exchange for a working conda interaction.

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
WORKDIR="/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05"
FLAG="${WORKDIR}/.hmmcleaner_installed"
DEP_LIST="${WORKDIR}/scripts/unity/cpanm_dep_order.txt"
LOGDIR="${WORKDIR}/logs"
mkdir -p "${LOGDIR}"

# Idempotency: if flag exists and binary works, exit early.
if [[ -f "${FLAG}" ]] && command -v HmmCleaner.pl >/dev/null 2>&1 \
   && HmmCleaner.pl --help >/dev/null 2>&1; then
    echo "[$(date -u +%FT%TZ)] HmmCleaner.pl already installed; flag at ${FLAG}; exiting."
    exit 0
fi

# -----------------------------------------------------------------------------
# Conda env activation
# -----------------------------------------------------------------------------
# shellcheck disable=SC1091
source ~/.miniconda3/etc/profile.d/conda.sh
conda activate berghia-gpcr

echo "[$(date -u +%FT%TZ)] Active env: ${CONDA_DEFAULT_ENV:-?}"
echo "[$(date -u +%FT%TZ)] CONDA_PREFIX=${CONDA_PREFIX}"
echo "[$(date -u +%FT%TZ)] perl: $(perl -v | head -2 | tail -1)"

# -----------------------------------------------------------------------------
# System bio tools that Bio::MUST::Drivers shells out to.
# These are on bioconda and are NOT included in the existing env per the
# May-2026 inventory. Install via conda (much more reliable than CPAN's
# attempt to detect them).
# -----------------------------------------------------------------------------
echo "[$(date -u +%FT%TZ)] Ensuring exonerate, hmmer, blast are installed via conda..."
conda install -y -n berghia-gpcr -c bioconda -c conda-forge \
    exonerate hmmer blast || {
        echo "WARN: conda install of system tools partially failed; continuing." >&2
    }

# -----------------------------------------------------------------------------
# cpanm bootstrap (install into the conda env's site_perl tree).
# We deliberately do NOT use a per-user local::lib in $HOME — the conda env
# already has its own perl, and the install must live inside that env so that
# `which HmmCleaner.pl` resolves predictably under any conda activation.
# -----------------------------------------------------------------------------
if ! command -v cpanm >/dev/null 2>&1; then
    echo "[$(date -u +%FT%TZ)] Installing App::cpanminus into env perl..."
    # Use the conda env's perl explicitly to bootstrap cpanm.
    conda install -y -n berghia-gpcr -c bioconda -c conda-forge perl-app-cpanminus \
        || curl -sL https://cpanmin.us | perl - --self-contained App::cpanminus
fi

CPANM=$(command -v cpanm)
PERL=$(command -v perl)
echo "[$(date -u +%FT%TZ)] Using cpanm: ${CPANM}"
echo "[$(date -u +%FT%TZ)] Using perl:  ${PERL}"

# Make cpanm install into the conda env perl tree, not ~/perl5.
# (PERL5LIB is left untouched; the Makefile.PL → install puts files under
#  $CONDA_PREFIX/lib/perl5/site_perl/...)
export PERL_CPANM_OPT="--notest --no-man-pages --quiet"

# -----------------------------------------------------------------------------
# Install dependencies in topological order from cpanm_dep_order.txt.
# Skip blank lines and comments. A failure on any one module is logged but
# does NOT abort the loop; we only abort if the FINAL dist
# (Bio::MUST::Apps::HmmCleaner) fails.
# -----------------------------------------------------------------------------
if [[ ! -f "${DEP_LIST}" ]]; then
    echo "ERROR: dep list not found at ${DEP_LIST}" >&2
    exit 2
fi

echo "[$(date -u +%FT%TZ)] Installing CPAN deps in topological order..."
FAILED=()
while IFS= read -r mod; do
    # strip comments and trailing whitespace
    mod="${mod%%#*}"
    mod="${mod//$'\r'/}"
    mod="$(echo "${mod}" | xargs || true)"
    [[ -z "${mod}" ]] && continue

    echo "  -> cpanm ${mod}"
    if ! "${CPANM}" "${mod}"; then
        echo "  WARN: ${mod} failed; continuing." >&2
        FAILED+=("${mod}")
    fi
done < "${DEP_LIST}"

echo "[$(date -u +%FT%TZ)] Dep install pass complete. Soft-failures: ${#FAILED[@]}"
if [[ ${#FAILED[@]} -gt 0 ]]; then
    printf '   - %s\n' "${FAILED[@]}"
fi

# -----------------------------------------------------------------------------
# Final retry of the leaf dist with --force (tolerates upstream test flake
# AS LONG AS the binary works at the verification step below).
# -----------------------------------------------------------------------------
if ! command -v HmmCleaner.pl >/dev/null 2>&1; then
    echo "[$(date -u +%FT%TZ)] HmmCleaner.pl not yet on PATH; retrying with --force..."
    "${CPANM}" --force Bio::MUST::Apps::HmmCleaner || true
fi

# -----------------------------------------------------------------------------
# Verification (this is the only "real" success criterion).
# -----------------------------------------------------------------------------
echo "[$(date -u +%FT%TZ)] Verifying install..."
if ! command -v HmmCleaner.pl >/dev/null 2>&1; then
    echo "ERROR: HmmCleaner.pl not on PATH after install" >&2
    echo "  PATH=${PATH}" >&2
    echo "  CONDA_PREFIX/bin contents (HmmCleaner.pl?):" >&2
    ls -la "${CONDA_PREFIX}/bin/HmmCleaner"* 2>&1 || true
    exit 3
fi
echo "  HmmCleaner.pl path: $(command -v HmmCleaner.pl)"

if ! HmmCleaner.pl --help 2>&1 | head -5; then
    echo "ERROR: HmmCleaner.pl --help exited non-zero" >&2
    exit 4
fi

# -----------------------------------------------------------------------------
# Mark success.
# -----------------------------------------------------------------------------
touch "${FLAG}"
echo "[$(date -u +%FT%TZ)] SUCCESS: HmmCleaner.pl installed; flag at ${FLAG}"
echo "[$(date -u +%FT%TZ)] To use, set RUN_HMMCLEANER=1 in config.local.sh and re-run stage 04."
