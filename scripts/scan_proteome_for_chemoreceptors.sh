#!/usr/bin/env bash
# scan_proteome_for_chemoreceptors.sh
#
# Per-proteome 6TM chemoreceptor scanner.
#
# Wraps the stage-02 inner-loop flow (HMM-first GPCR filter → length filter →
# TMbed/DeepTMHMM TM prediction → ≥6 TM filter) as a standalone command that
# can be invoked on ANY input proteome FASTA, not just the Berghia focal sample.
#
# Used by P1 (GPCR class classifier) and P2 (per-class pool builder).
#
# Usage — explicit mode:
#   bash scan_proteome_for_chemoreceptors.sh \
#       --proteome <fasta> --out-dir <dir> --taxid <int> --binomial <str> \
#       [--min-tm 6] [--min-confidence 0.5] [--threads N] [--force]
#
# Usage — SLURM array mode (SLURM_ARRAY_TASK_ID must be set):
#   bash scan_proteome_for_chemoreceptors.sh \
#       --manifest <tsv> --out-dir <dir> \
#       [--min-tm 6] [--min-confidence 0.5] [--threads N] [--force]
#
# Manifest TSV columns (header required): taxid, binomial, proteome_path
#
# Outputs (written to <out_dir>/):
#   <taxid>_<sanitized_binomial>.chemo_candidates.fa   — FASTA of candidates
#   <taxid>_<sanitized_binomial>.scan_record.tsv       — per-sequence record
#   <taxid>_<sanitized_binomial>.zero_candidates       — flag if 0 candidates
#
# Idempotent: skips if both outputs already exist (with candidates) unless --force.

set -eo pipefail

# ---------------------------------------------------------------------------
# Locate config.sh + functions.sh relative to this script
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

if [[ -f "${BASE_DIR}/config.sh" ]]; then
    # shellcheck source=../config.sh
    source "${BASE_DIR}/config.sh"
fi
if [[ -f "${BASE_DIR}/functions.sh" ]]; then
    # shellcheck source=../functions.sh
    source "${BASE_DIR}/functions.sh"
fi

# ---------------------------------------------------------------------------
# Default values (overridden by config.sh env vars, then by CLI flags)
# ---------------------------------------------------------------------------

ARG_PROTEOME=""
ARG_MANIFEST=""
ARG_OUT_DIR=""
ARG_TAXID=""
ARG_BINOMIAL=""
ARG_MIN_TM="${MIN_TM_REGIONS:-6}"
ARG_MIN_CONF="${DEEPTMHMM_MIN_CONFIDENCE:-0.5}"
ARG_THREADS="${SLURM_CPUS_PER_TASK:-${CPUS:-4}}"
ARG_FORCE=0

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

usage() {
    cat >&2 <<EOF
Usage: $(basename "$0") [OPTIONS]

Explicit mode:
  --proteome FILE     Input proteome FASTA (required unless --manifest used)
  --out-dir DIR       Output directory (required)
  --taxid INT         NCBI taxid for this proteome (required in explicit mode)
  --binomial STR      Binomial species name (required in explicit mode)

SLURM array mode (reads row SLURM_ARRAY_TASK_ID from manifest):
  --manifest FILE     TSV with columns: taxid, binomial, proteome_path
  --out-dir DIR       Output directory (required)

Common options:
  --min-tm N          Minimum TM regions [default: ${ARG_MIN_TM}]
  --min-confidence F  Minimum TMbed confidence [default: ${ARG_MIN_CONF}]
  --threads N         CPU threads [default: ${ARG_THREADS}]
  --force             Re-run even if outputs already exist
  -h, --help          Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --proteome)       ARG_PROTEOME="$2";  shift 2 ;;
        --manifest)       ARG_MANIFEST="$2";  shift 2 ;;
        --out-dir)        ARG_OUT_DIR="$2";   shift 2 ;;
        --taxid)          ARG_TAXID="$2";     shift 2 ;;
        --binomial)       ARG_BINOMIAL="$2";  shift 2 ;;
        --min-tm)         ARG_MIN_TM="$2";    shift 2 ;;
        --min-confidence) ARG_MIN_CONF="$2";  shift 2 ;;
        --threads)        ARG_THREADS="$2";   shift 2 ;;
        --force)          ARG_FORCE=1;        shift ;;
        -h|--help)        usage ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done

# ---------------------------------------------------------------------------
# SLURM array manifest resolution
# ---------------------------------------------------------------------------

if [[ -n "${ARG_MANIFEST}" ]]; then
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        echo "[scan] ERROR: --manifest requires SLURM_ARRAY_TASK_ID to be set" >&2
        exit 1
    fi
    if [[ ! -f "${ARG_MANIFEST}" ]]; then
        echo "[scan] ERROR: manifest not found: ${ARG_MANIFEST}" >&2
        exit 1
    fi
    # Row 1 (after header) = SLURM_ARRAY_TASK_ID 1, etc.
    TASK_ROW=$(awk -F'\t' -v n="${SLURM_ARRAY_TASK_ID}" \
        'NR == 1 { for (i=1;i<=NF;i++) col[$i]=i; next }
         NR == n+1 { print $col["taxid"] "\t" $col["binomial"] "\t" $col["proteome_path"]; exit }
        ' "${ARG_MANIFEST}")
    if [[ -z "${TASK_ROW}" ]]; then
        echo "[scan] ERROR: no row ${SLURM_ARRAY_TASK_ID} in manifest (1-indexed: row 1 = task 1)" >&2
        exit 1
    fi
    ARG_TAXID=$(echo "${TASK_ROW}"    | cut -f1)
    ARG_BINOMIAL=$(echo "${TASK_ROW}" | cut -f2)
    ARG_PROTEOME=$(echo "${TASK_ROW}" | cut -f3)
fi

# ---------------------------------------------------------------------------
# Validate required args
# ---------------------------------------------------------------------------

if [[ -z "${ARG_OUT_DIR}" ]]; then
    echo "[scan] ERROR: --out-dir is required" >&2; usage
fi
if [[ -z "${ARG_PROTEOME}" ]]; then
    echo "[scan] ERROR: --proteome (or --manifest) is required" >&2; usage
fi
if [[ -z "${ARG_TAXID}" ]]; then
    echo "[scan] ERROR: --taxid is required" >&2; usage
fi
if [[ -z "${ARG_BINOMIAL}" ]]; then
    echo "[scan] ERROR: --binomial is required" >&2; usage
fi
if [[ ! -f "${ARG_PROTEOME}" ]]; then
    echo "[scan] ERROR: proteome file not found: ${ARG_PROTEOME}" >&2; exit 1
fi

# ---------------------------------------------------------------------------
# Derive the sanitized sample name (matches _scan_proteome_filter_helper.py)
# Bash-native equivalent of Python: re.sub(r'[^A-Za-z0-9_]+', '_', s).strip('_')
# with run deduplication.  Used only to build output file paths.
# ---------------------------------------------------------------------------

_sanitize() {
    local s="$1"
    # Replace non-[A-Za-z0-9_] runs with '_', then collapse '__+', strip edges
    echo "$s" | sed 's/[^A-Za-z0-9_]\+/_/g; s/_\+/_/g; s/^_//; s/_$//'
}

SAMPLE_NAME="$(_sanitize "${ARG_TAXID}_${ARG_BINOMIAL}")"

mkdir -p "${ARG_OUT_DIR}"

OUT_FA="${ARG_OUT_DIR}/${SAMPLE_NAME}.chemo_candidates.fa"
OUT_TSV="${ARG_OUT_DIR}/${SAMPLE_NAME}.scan_record.tsv"
ZERO_FLAG="${ARG_OUT_DIR}/${SAMPLE_NAME}.zero_candidates"
WORK_DIR="${ARG_OUT_DIR}/_work_${SAMPLE_NAME}"

echo "[scan] Sample: ${SAMPLE_NAME}" >&2
echo "[scan] Proteome: ${ARG_PROTEOME}" >&2
echo "[scan] Output dir: ${ARG_OUT_DIR}" >&2

# ---------------------------------------------------------------------------
# Idempotency check
# ---------------------------------------------------------------------------

if [[ "${ARG_FORCE}" -eq 0 ]]; then
    HAS_CANDS=0
    if [[ -f "${OUT_FA}" ]] && grep -q '^>' "${OUT_FA}" 2>/dev/null; then
        HAS_CANDS=1
    fi
    HAS_ZERO=0
    [[ -f "${ZERO_FLAG}" ]] && HAS_ZERO=1
    if [[ -f "${OUT_FA}" && -f "${OUT_TSV}" && ( "${HAS_CANDS}" -eq 1 || "${HAS_ZERO}" -eq 1 ) ]]; then
        echo "[scan] skipped: outputs already exist for ${SAMPLE_NAME}" >&2
        exit 0
    fi
fi

mkdir -p "${WORK_DIR}"

GPCR_FA="${WORK_DIR}/gpcrs.fa"
TMBED_INPUT_FA="${WORK_DIR}/tmbed_input.aa"
TMBED_OUT_DIR="${WORK_DIR}/tmbed_out"
HMM_IDS_FILE="${WORK_DIR}/hmm_positive_ids.txt"

# ---------------------------------------------------------------------------
# Step 1: HMM-first GPCR filter → GPCR-positive subset
#
# Calls identify_gpcr_candidates from functions.sh (same function used by
# stage 02 for non-Berghia reference samples).  The function runs
# classify_via_hmm.py (curated Swiss-Prot HMMs + Pfam fallback) +
# hmmsearch against TIAMMAT_GPCR_HMM + curated lophotrochozoan HMMs,
# merges via identify_gpcrs.py, and writes a GPCR-positive FASTA.
# ---------------------------------------------------------------------------

echo "[scan] Step 1: HMM-first GPCR filter" >&2

identify_gpcr_candidates "${ARG_PROTEOME}" "${GPCR_FA}"

if [[ ! -s "${GPCR_FA}" ]]; then
    echo "[scan] WARN: no GPCR-positive sequences found; writing empty outputs" >&2
    : > "${OUT_FA}"
    printf 'seq_id\tgpcr_positive\ttm_count\ttm_confidence\tpassed_gate\tnotes\n' > "${OUT_TSV}"
    touch "${ZERO_FLAG}"
    exit 0
fi

# Collect the HMM-positive IDs for the Python helper
grep '^>' "${GPCR_FA}" | sed 's/^>//' | awk '{print $1}' > "${HMM_IDS_FILE}"
N_GPCR=$(wc -l < "${HMM_IDS_FILE}")
echo "[scan] Step 1: ${N_GPCR} GPCR-positive sequences" >&2

# ---------------------------------------------------------------------------
# Step 2: Length pre-filter — drop sequences > MAX_AA_LENGTH before TMbed
# (mirrors stage-02 filter_fasta_by_length)
# ---------------------------------------------------------------------------

echo "[scan] Step 2: length filter (max ${MAX_AA_LENGTH:-1500} aa)" >&2

filter_fasta_by_length "${GPCR_FA}" "${TMBED_INPUT_FA}" "${MAX_AA_LENGTH:-1500}"

if [[ ! -s "${TMBED_INPUT_FA}" ]]; then
    echo "[scan] WARN: no sequences survived length filter; writing empty outputs" >&2
    : > "${OUT_FA}"
    printf 'seq_id\tgpcr_positive\ttm_count\ttm_confidence\tpassed_gate\tnotes\n' > "${OUT_TSV}"
    touch "${ZERO_FLAG}"
    exit 0
fi

# ---------------------------------------------------------------------------
# Step 3: TMbed / DeepTMHMM on the GPCR-positive, length-OK subset
# ---------------------------------------------------------------------------

echo "[scan] Step 3: TMbed/DeepTMHMM" >&2

TMBED_PREDICTION="${TMBED_OUT_DIR}/prediction"

mkdir -p "${TMBED_OUT_DIR}"
${DEEPTMHMM:-bash "${SCRIPT_DIR}/run_deeptmhmm.sh"} \
    -f "${TMBED_INPUT_FA}" \
    -o "${TMBED_OUT_DIR}" \
    2> "${WORK_DIR}/tmbed.err" || {
        echo "[scan] ERROR: TMbed/DeepTMHMM failed; see ${WORK_DIR}/tmbed.err" >&2
        exit 1
    }

if [[ ! -f "${TMBED_PREDICTION}" ]]; then
    echo "[scan] ERROR: TMbed prediction file not found: ${TMBED_PREDICTION}" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 4: Apply ≥min_tm + ≥min_confidence filter + emit outputs via Python helper
#
# The helper (scripts/_scan_proteome_filter_helper.py):
#   - Reproduces stage-02 awk filter exactly
#     (NF >= 5 && $5+0 >= min_tm && $3+0 >= min_conf)
#   - Intersects TMbed-positive set with HMM-positive set
#   - Writes the 6-column scan_record.tsv
#   - Writes the candidate FASTA
#   - Writes a zero_candidates flag if 0 candidates pass
# ---------------------------------------------------------------------------

echo "[scan] Step 4: filter + emit outputs" >&2

FORCE_FLAG=""
[[ "${ARG_FORCE}" -eq 1 ]] && FORCE_FLAG="--force"

python3 "${SCRIPT_DIR}/_scan_proteome_filter_helper.py" \
    --tmbed-prediction "${TMBED_PREDICTION}" \
    --hmm-positive-ids "${HMM_IDS_FILE}" \
    --input-fa "${GPCR_FA}" \
    --out-fa "${OUT_FA}" \
    --out-tsv "${OUT_TSV}" \
    --min-tm "${ARG_MIN_TM}" \
    --min-confidence "${ARG_MIN_CONF}" \
    ${FORCE_FLAG}

N_CANDS=$(grep -c '^>' "${OUT_FA}" 2>/dev/null || echo 0)
echo "[scan] Done: ${N_CANDS} chemoreceptor candidates for ${SAMPLE_NAME}" >&2
echo "[scan] FASTA: ${OUT_FA}" >&2
echo "[scan] TSV:   ${OUT_TSV}" >&2
