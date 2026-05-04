#!/bin/bash
# run_cds_preprocess.sh — One-shot preprocessing for reference CDS recovery.
#
# Bead -325. Downloads the 55 reference-species genome assemblies needed
# for miniprot-based CDS recovery, runs the recovery, and writes the
# merged reference CDS file + per-sequence provenance manifest that
# stage 05's find_nucleotide_sequences expects to read.
#
# Run this ONCE before the main pipeline (stages 01-09). Idempotent:
# - Genomes are cached under genomes/ref_assemblies/ and reused across runs
#   (use --redownload to force re-fetch; or delete genomes/ref_assemblies/).
# - Per-species CDS files are skipped if already present.
# - Safe to re-run after adding new species to SPECIES_MAP.
#
# Usage:
#   bash scripts/run_cds_preprocess.sh                    # all species
#   bash scripts/run_cds_preprocess.sh --species cobe,dipe  # subset
#   bash scripts/run_cds_preprocess.sh --threads 16       # parallelism
#   bash scripts/run_cds_preprocess.sh --redownload       # force re-fetch
#
# Expected runtime: 2-8 hours depending on network + CPU. Output sizes:
# ~50GB of cached genomes (kept), ~50MB of recovered CDS + manifest.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck disable=SC1091
source "$BASE_DIR/config.sh"

THREADS="${CPUS:-8}"
SPECIES=""
REDOWNLOAD=0
KEEP_GENOMES=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads=*|--threads)
            if [[ "$1" == --threads=* ]]; then THREADS="${1#*=}"; shift
            else THREADS="$2"; shift 2; fi ;;
        --species=*|--species)
            if [[ "$1" == --species=* ]]; then SPECIES="${1#*=}"; shift
            else SPECIES="$2"; shift 2; fi ;;
        --redownload)        REDOWNLOAD=1; shift ;;
        --no-keep-genomes)   KEEP_GENOMES=0; shift ;;
        -h|--help)
            sed -n '2,30p' "$0" | sed 's/^# //;s/^#//' >&2
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

GENOME_CACHE="${BASE_DIR}/genomes/ref_assemblies"
RECOVERED_DIR="${RESULTS_DIR}/reference_sequences/cds_recovered"
NATIVE_CDS_DIR="${RESULTS_DIR}/reference_sequences/cds"
MANIFEST="${RESULTS_DIR}/reference_sequences/cds_provenance.csv"
mkdir -p "$GENOME_CACHE" "$RECOVERED_DIR" "$NATIVE_CDS_DIR"

if [ "$REDOWNLOAD" = "1" ]; then
    echo "[run_cds_preprocess] --redownload set; clearing $GENOME_CACHE" >&2
    rm -rf "$GENOME_CACHE"
    mkdir -p "$GENOME_CACHE"
fi

# --- Step 1: efetch-based fast path (fetch_reference_cds.py) ---
# This catches ~4,000 references whose protein record has a linked NUCCORE
# CDS — much faster than miniprot. Idempotent (skips if already done).
NATH_DIR="${NATH_ET_AL_DIR:-${REFERENCE_DIR}/nath_et_al}"
if [ -d "$NATH_DIR" ] && [ -f "$SCRIPT_DIR/fetch_reference_cds.py" ]; then
    if ! find "$NATIVE_CDS_DIR" -name "*_cds.fna" -type f 2>/dev/null | grep -q .; then
        echo "[run_cds_preprocess] Step 1/3: efetch CDS for reference proteins..." >&2
        python3 "$SCRIPT_DIR/fetch_reference_cds.py" "$NATH_DIR" \
            -o "$NATIVE_CDS_DIR" \
            --log-failed "${RESULTS_DIR}/reference_sequences/failed_cds_accessions.csv"
    else
        echo "[run_cds_preprocess] Step 1/3: efetch CDS already present (skipping)" >&2
    fi
else
    echo "[run_cds_preprocess] Step 1/3: skipped (Nath dir or fetch script not found)" >&2
fi

# --- Step 2: miniprot recovery for the remaining ~22k references ---
echo "[run_cds_preprocess] Step 2/3: miniprot CDS recovery (cached at $GENOME_CACHE)..." >&2
EXISTING_CDS=""
if [ -f "${NATIVE_CDS_DIR}/all_references_cds.fna" ]; then
    EXISTING_CDS="${NATIVE_CDS_DIR}/all_references_cds.fna"
else
    # Pre-concatenate per-species efetch results so the recovery script
    # can identify which proteins are still missing.
    EXISTING_CDS="${NATIVE_CDS_DIR}/all_references_cds.fna"
    : > "$EXISTING_CDS"
    find "$NATIVE_CDS_DIR" -name "*_cds.fna" -type f \
        ! -name "all_references_cds.fna" -exec cat {} + >> "$EXISTING_CDS" 2>/dev/null || true
fi

species_arg=""
[ -n "$SPECIES" ] && species_arg="--species $SPECIES"
keep_arg=""
[ "$KEEP_GENOMES" = "1" ] && keep_arg="--keep-genomes"

python3 "$SCRIPT_DIR/recover_cds_from_assemblies.py" \
    --og-dir "$NATH_DIR" \
    --cds-file "$EXISTING_CDS" \
    --output-dir "$RECOVERED_DIR" \
    --genome-dir "$GENOME_CACHE" \
    --threads "$THREADS" \
    --manifest "$MANIFEST" \
    $keep_arg \
    $species_arg

# --- Step 3: merge native + recovered into the single CDS file stage 05 reads ---
echo "[run_cds_preprocess] Step 3/3: merging native + recovered CDS..." >&2
ALL_REF_CDS="${NATIVE_CDS_DIR}/all_references_cds.fna"
: > "${ALL_REF_CDS}.tmp"
find "$NATIVE_CDS_DIR" -name "*.fna" -type f \
    ! -name "all_references_cds.fna*" -exec cat {} + >> "${ALL_REF_CDS}.tmp" 2>/dev/null || true
find "$RECOVERED_DIR" -name "*.fna" -type f \
    ! -name "all_recovered_cds.fna*" -exec cat {} + >> "${ALL_REF_CDS}.tmp" 2>/dev/null || true
mv "${ALL_REF_CDS}.tmp" "$ALL_REF_CDS"

n_seq=$(grep -c '^>' "$ALL_REF_CDS" 2>/dev/null || echo 0)
n_native=$(awk -F, 'NR>1 && $2=="native"' "$MANIFEST" 2>/dev/null | wc -l || echo 0)
n_miniprot=$(awk -F, 'NR>1 && $2=="miniprot"' "$MANIFEST" 2>/dev/null | wc -l || echo 0)

echo "" >&2
echo "[run_cds_preprocess] DONE. Summary:" >&2
echo "  Merged CDS file:         $ALL_REF_CDS" >&2
echo "  Total CDS sequences:     $n_seq" >&2
echo "  Provenance — miniprot:   $n_miniprot   (recovered via miniprot)" >&2
echo "  Provenance — native:     $n_native     (efetch / pre-existing)" >&2
echo "  Provenance manifest:     $MANIFEST" >&2
echo "  Genome cache:            $GENOME_CACHE" >&2
echo "" >&2
echo "Now ready to start the main pipeline (sbatch 01_reference_processing.sh)." >&2
