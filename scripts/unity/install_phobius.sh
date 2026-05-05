#!/bin/bash
# install_phobius.sh — install Phobius (TM + signal-peptide HMM predictor)
# into the berghia-gpcr workdir on Unity, as the tertiary fallback in the
# pipeline's TM-prediction chain.
#
# WHY THIS EXISTS
# ---------------
# The pipeline's TM-prediction chain is:
#   1. TMbed   (primary)   — installed via git clone of BernhoferM/TMbed; works.
#   2. DeepTMHMM (secondary, three sub-strategies: Apptainer / venv / system)
#                          — none currently installed; not blocking.
#   3. Phobius (tertiary)  — HMM-based, classical, very robust.
#                          — Bioconda install attempted in job 56714199 with
#                            `perl-bio-must-...` (the WRONG package name) and
#                            failed. There is no real Phobius package on the
#                            main bioconda channel — see "license note" below.
#
# Phobius (Käll, Krogh, Sonnhammer 2004 J Mol Biol 338:1027-36,
# doi:10.1016/j.jmb.2004.03.016; web service Käll et al. 2007 NAR
# 35:W429-W432, doi:10.1093/nar/gkm256) ships under an academic-use license
# from Stockholm University (sbc.su.se). The license requires explicit human
# acceptance via a web form before the tarball is released:
#   https://phobius.sbc.su.se/data.html
#   http://software.sbc.su.se/cgi-bin/request.cgi?project=phobius
# This is why no automated pull can fetch Phobius — bioconda's main channel
# does not host it, and the only third-party mirror (the `predector` channel
# on anaconda.org) ships an "empty placeholder" recipe whose `phobius-register`
# script still requires the user to provide the manually-downloaded tarball.
# We therefore CANNOT auto-download. The script below requires the user to
# place the tarball in a known location after accepting the license; see
# README_phobius_install.md for the step-by-step.
#
# WHY NO CONTAINER FALLBACK
# -------------------------
# Same root cause: redistributing the Phobius binary in a public container
# would violate the academic-use license, so no public quay.io / DockerHub /
# Galaxy depot image exists. (InterProScan ships Phobius support but only
# *activates* it after the user drops the same tarball into its data/ tree —
# see ebi-pf-team/interproscan-docs ActivatingLicensedAnalyses.rst.)
#
# CHOSEN APPROACH
# ---------------
# 1. User accepts license + downloads `phobius101_linux.tar.gz` manually,
#    scp's it to `${WORKDIR}/tools/phobius/phobius_dist.tar.gz`.
# 2. This script untars to `${WORKDIR}/tools/phobius/dist/`, normalizes the
#    decodeanhmm binary (the tarball ships a 32-bit `decodeanhmm` and a
#    `decodeanhmm.64bit`; on x86_64 we promote .64bit to be the canonical
#    one — this is the documented Phobius fix for the well-known
#    "could not read provided fasta sequence" bug on 64-bit Linux,
#    biostars.org/p/238642 and predector#82).
# 3. Writes a shim wrapper at `~/.local/bin/phobius.pl` that exec's the
#    unpacked phobius.pl with PATH adjusted so its internal call to
#    decodeanhmm resolves. `~/.local/bin` is on the conda env's PATH per the
#    existing berghia-gpcr setup; if it isn't, the script prints a
#    diagnostic and exits 5.
# 4. Smoke-tests the wrapper on a single GPCR-like sequence shipped beside
#    this script (`phobius_smoke_test.fasta`).
#
# SUBMIT (after the manual download step in the README):
#   sbatch scripts/unity/install_phobius.sh
#
# After completion, the pipeline's TM-prediction chain will find phobius.pl
# on PATH and use it as the third-tier fallback when both TMbed and DeepTMHMM
# are unavailable or fail.

#SBATCH --job-name=phobius_install
#SBATCH --partition=cpu
#SBATCH --time=30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
WORKDIR="/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05"
PHOBIUS_ROOT="${WORKDIR}/tools/phobius"
PHOBIUS_TARBALL="${PHOBIUS_ROOT}/phobius_dist.tar.gz"
PHOBIUS_DIST="${PHOBIUS_ROOT}/dist"
WRAPPER_DIR="${HOME}/.local/bin"
WRAPPER="${WRAPPER_DIR}/phobius.pl"
FLAG="${WORKDIR}/.phobius_installed"
LOGDIR="${WORKDIR}/logs"
SMOKE_FASTA_SRC="${WORKDIR}/scripts/unity/phobius_smoke_test.fasta"
SMOKE_FASTA_TMP="${LOGDIR}/phobius_smoke_test.fasta"

mkdir -p "${LOGDIR}" "${PHOBIUS_ROOT}" "${WRAPPER_DIR}"

# -----------------------------------------------------------------------------
# Idempotency: if flag exists and wrapper works on the smoke test, exit early.
# -----------------------------------------------------------------------------
if [[ -f "${FLAG}" ]] \
   && command -v phobius.pl >/dev/null 2>&1 \
   && phobius.pl -h >/dev/null 2>&1; then
    echo "[$(date -u +%FT%TZ)] Phobius already installed; flag at ${FLAG}; exiting."
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
echo "[$(date -u +%FT%TZ)] perl: $(perl -v 2>&1 | head -2 | tail -1)"

# -----------------------------------------------------------------------------
# Tarball presence check. We MUST NOT auto-download — see header note.
# -----------------------------------------------------------------------------
if [[ ! -f "${PHOBIUS_TARBALL}" ]]; then
    cat >&2 <<EOF
ERROR: Phobius distribution tarball not found at:
       ${PHOBIUS_TARBALL}

Phobius is academic-licensed software (Käll et al. 2004, Stockholm
University) and the license requires manual acceptance before download.
This script cannot automate that step. To proceed:

  1. On a workstation with a browser, visit:
       https://phobius.sbc.su.se/data.html
     and follow the link to the request form
       http://software.sbc.su.se/cgi-bin/request.cgi?project=phobius
     Accept the academic-use license. The site emails (or directly serves)
     the tarball "phobius101_linux.tar.gz" (~2 MB).

  2. Rename and scp it to Unity at the expected path:
       scp phobius101_linux.tar.gz \\
         unity:${PHOBIUS_TARBALL}

  3. Re-submit this job:
       sbatch scripts/unity/install_phobius.sh

See README_phobius_install.md beside this script for full instructions.
EOF
    exit 2
fi

# -----------------------------------------------------------------------------
# Untar.
# The tarball lays out as: phobius/phobius.pl, phobius/decodeanhmm,
# phobius/decodeanhmm.64bit, phobius/phobius.model, phobius/phobius.options,
# phobius/example.fasta (in some versions).
# We unpack into ${PHOBIUS_DIST}/ so that ${PHOBIUS_DIST}/phobius/phobius.pl
# is the canonical entrypoint.
# -----------------------------------------------------------------------------
echo "[$(date -u +%FT%TZ)] Unpacking ${PHOBIUS_TARBALL} into ${PHOBIUS_DIST}..."
rm -rf "${PHOBIUS_DIST}"
mkdir -p "${PHOBIUS_DIST}"
tar -xzf "${PHOBIUS_TARBALL}" -C "${PHOBIUS_DIST}"

# Locate the unpacked phobius dir robustly: most releases use "phobius/", but
# some intermediate releases used "phobius101_linux/". Glob and pick the first
# directory that contains phobius.pl.
PHOBIUS_BIN_DIR=""
for cand in "${PHOBIUS_DIST}"/*/; do
    if [[ -f "${cand}phobius.pl" ]]; then
        PHOBIUS_BIN_DIR="${cand%/}"
        break
    fi
done
if [[ -z "${PHOBIUS_BIN_DIR}" ]]; then
    echo "ERROR: tarball unpacked but no */phobius.pl found under ${PHOBIUS_DIST}" >&2
    echo "Contents:" >&2
    ls -la "${PHOBIUS_DIST}" >&2 || true
    exit 3
fi
echo "[$(date -u +%FT%TZ)] Phobius binary dir: ${PHOBIUS_BIN_DIR}"

# -----------------------------------------------------------------------------
# 64-bit normalisation.
# The phobius.pl bundled in the tarball calls "decodeanhmm" (no extension).
# On x86_64 hosts, the shipped 32-bit binary segfaults / silently produces
# the "could not read provided fasta sequence at phobius.pl line 408" error
# (biostars.org/p/238642, predector#82, funannotate#696). The fix is to
# use decodeanhmm.64bit — promote it to be the canonical decodeanhmm.
# -----------------------------------------------------------------------------
if [[ -f "${PHOBIUS_BIN_DIR}/decodeanhmm.64bit" ]]; then
    echo "[$(date -u +%FT%TZ)] Promoting decodeanhmm.64bit -> decodeanhmm (x86_64 fix)"
    cp -f "${PHOBIUS_BIN_DIR}/decodeanhmm.64bit" "${PHOBIUS_BIN_DIR}/decodeanhmm"
    chmod +x "${PHOBIUS_BIN_DIR}/decodeanhmm"
elif [[ -f "${PHOBIUS_BIN_DIR}/decodeanhmm" ]]; then
    echo "[$(date -u +%FT%TZ)] No decodeanhmm.64bit shipped; using decodeanhmm as-is."
    chmod +x "${PHOBIUS_BIN_DIR}/decodeanhmm"
else
    echo "ERROR: no decodeanhmm binary in ${PHOBIUS_BIN_DIR}" >&2
    ls -la "${PHOBIUS_BIN_DIR}" >&2 || true
    exit 4
fi
chmod +x "${PHOBIUS_BIN_DIR}/phobius.pl"

# Sanity check: model + options files are present (phobius.pl bails without).
for f in phobius.model phobius.options; do
    if [[ ! -f "${PHOBIUS_BIN_DIR}/${f}" ]]; then
        echo "ERROR: missing required Phobius data file: ${PHOBIUS_BIN_DIR}/${f}" >&2
        exit 4
    fi
done

# -----------------------------------------------------------------------------
# Wrapper: phobius.pl on PATH.
# phobius.pl resolves decodeanhmm and the model files via $0's directory
# (it does `use FindBin; my $PHOBIUS_DIR = $FindBin::Bin;`). So the wrapper
# must `exec` the unpacked phobius.pl directly, preserving $0's dirname —
# a symlink works here, but we prefer a tiny shell shim so the path is
# explicit in `command -v` output and easier to grep for in logs.
# -----------------------------------------------------------------------------
cat > "${WRAPPER}" <<EOF
#!/bin/bash
# Auto-generated by scripts/unity/install_phobius.sh on $(date -u +%FT%TZ).
# Do not edit by hand; rerun the install script to regenerate.
exec "${PHOBIUS_BIN_DIR}/phobius.pl" "\$@"
EOF
chmod +x "${WRAPPER}"
echo "[$(date -u +%FT%TZ)] Wrapper installed at ${WRAPPER}"

# Confirm ~/.local/bin is on PATH inside this conda env.
if ! command -v phobius.pl >/dev/null 2>&1; then
    echo "ERROR: phobius.pl wrapper not on PATH after install." >&2
    echo "  WRAPPER=${WRAPPER}" >&2
    echo "  PATH=${PATH}" >&2
    echo "Add ${WRAPPER_DIR} to PATH (e.g. via ~/.bashrc or conda env_vars) and rerun." >&2
    exit 5
fi

# -----------------------------------------------------------------------------
# Smoke test.
# Use the bundled smoke FASTA if present; otherwise fall back to whatever
# example file the tarball shipped (some Phobius releases include
# example.fasta, others don't).
# -----------------------------------------------------------------------------
if [[ -f "${SMOKE_FASTA_SRC}" ]]; then
    cp -f "${SMOKE_FASTA_SRC}" "${SMOKE_FASTA_TMP}"
elif [[ -f "${PHOBIUS_BIN_DIR}/example.fasta" ]]; then
    cp -f "${PHOBIUS_BIN_DIR}/example.fasta" "${SMOKE_FASTA_TMP}"
else
    # Inline fallback: bovine rhodopsin (UniProt P02699) residues 1-100,
    # a canonical 7TM GPCR — Phobius should report >=1 TM helix in the
    # first ~100 aa.
    cat > "${SMOKE_FASTA_TMP}" <<'EOF'
>P02699_RHO_BOVIN_1_100 Rhodopsin Bos taurus (residues 1-100)
MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYV
TVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLH
EOF
fi

echo "[$(date -u +%FT%TZ)] Smoke test on ${SMOKE_FASTA_TMP}..."
SMOKE_OUT="${LOGDIR}/phobius_smoke_test.out"
if ! phobius.pl -short "${SMOKE_FASTA_TMP}" > "${SMOKE_OUT}" 2>&1; then
    echo "ERROR: phobius.pl smoke test exited non-zero." >&2
    echo "--- output (head) ---" >&2
    head -50 "${SMOKE_OUT}" >&2 || true
    exit 6
fi

# In `-short` mode Phobius emits one line per protein; the second column
# is the TM count. The bovine rhodopsin fragment has 1+ TMs in residues
# 1-100, so we require TM>=1. Bail with a clear error otherwise.
TM_COUNT=$(awk 'NR>1 && NF>=2 {print $2; exit}' "${SMOKE_OUT}" || echo "0")
if [[ -z "${TM_COUNT}" || "${TM_COUNT}" -lt 1 ]]; then
    echo "ERROR: Phobius smoke test reported TM=${TM_COUNT}; expected >=1." >&2
    echo "--- output ---" >&2
    cat "${SMOKE_OUT}" >&2 || true
    exit 7
fi
echo "[$(date -u +%FT%TZ)] Smoke test OK (TM count = ${TM_COUNT})"

# -----------------------------------------------------------------------------
# Mark success.
# -----------------------------------------------------------------------------
touch "${FLAG}"
echo "[$(date -u +%FT%TZ)] SUCCESS: Phobius installed."
echo "[$(date -u +%FT%TZ)]   binary dir: ${PHOBIUS_BIN_DIR}"
echo "[$(date -u +%FT%TZ)]   wrapper:    ${WRAPPER}"
echo "[$(date -u +%FT%TZ)]   flag:       ${FLAG}"
echo "[$(date -u +%FT%TZ)] Pipeline TM-fallback chain will pick up phobius.pl"
echo "[$(date -u +%FT%TZ)] when TMbed and DeepTMHMM are both unavailable."
