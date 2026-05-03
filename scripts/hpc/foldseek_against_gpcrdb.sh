#!/bin/bash
# foldseek_against_gpcrdb.sh — Structural augmentation stub (bead -s6v, P3).
#
# Compares AlphaFold-predicted Berghia structures against GPCRdb 2025
# (Pándy-Szekeres 2025 NAR 53:D425), which now hosts AlphaFold-Multistate
# models for 814 human ORs. Foldseek (Kempen 2023, Nat Biotechnol 41:243)
# clusters by structural similarity in seconds.
#
# This is a SCAFFOLD — install Foldseek + download GPCRdb structures
# yourself; the script orchestrates the comparison.

#SBATCH --job-name=foldseek_gpcrdb
#SBATCH --output=logs/foldseek_gpcrdb_%j.out
#SBATCH --error=logs/foldseek_gpcrdb_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$PROJECT_DIR/config.sh"

QUERY_DIR="${1:?Query AlphaFold PDB dir}"
GPCRDB_DB="${GPCRDB_FOLDSEEK_DB:-${PROJECT_DIR}/references/foldseek/gpcrdb_2025}"
OUT_DIR="${RESULTS_DIR}/foldseek"
mkdir -p "$OUT_DIR"

if ! command -v foldseek &>/dev/null; then
    echo "ERROR: foldseek not found. Install: " \
         "https://github.com/steineggerlab/foldseek/releases" >&2
    exit 2
fi
if [ ! -d "$GPCRDB_DB" ]; then
    echo "ERROR: GPCRdb foldseek DB not found at $GPCRDB_DB." >&2
    echo "       Build via: foldseek easy-search <gpcrdb_pdbs/> $GPCRDB_DB tmp/" >&2
    exit 2
fi

# Build query DB from candidate AlphaFold PDBs
QUERY_DB="$OUT_DIR/query_db"
foldseek createdb "$QUERY_DIR" "$QUERY_DB"

# Search against GPCRdb
foldseek search "$QUERY_DB" "$GPCRDB_DB" \
    "$OUT_DIR/result" "$OUT_DIR/tmp" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

# Convert to readable hits TSV
foldseek convertalis "$QUERY_DB" "$GPCRDB_DB" \
    "$OUT_DIR/result" "$OUT_DIR/hits.tsv" \
    --format-output "query,target,fident,alnlen,evalue,bits,prob"

echo "Foldseek vs GPCRdb complete: $OUT_DIR/hits.tsv" >&2
