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

# --- Refuse a zero-structure query set ---------------------------------------
# `foldseek createdb` on an empty (or simply wrong) directory does NOT fail: it
# builds an empty DB, the search returns nothing, convertalis writes a zero-row
# hits.tsv, and the job exits 0. Downstream that is indistinguishable from
# "searched everything, found no structural evidence" -- the silent dormancy
# bead 5ubd fixed in scripts/unity/run_foldseek_candidates.sh, which guards its
# own staging the same way (N_MODELS == 0 -> exit 3).
#
# Counted at depth 1 by the extensions foldseek actually ingests, and with -L so
# the SYMLINKS run_foldseek_candidates.sh stages (<cand_id>.cif -> the nested AF3
# model) are counted as the structures they point at. Requiring a structure
# extension -- not merely "some file" -- also catches the wrong-directory case:
# a FASTA dir, or stage 08's alphafold/ root whose depth-1 entries are all
# per-candidate DIRECTORIES.
if [ ! -d "$QUERY_DIR" ]; then
    echo "ERROR: query structure dir not found: $QUERY_DIR" >&2
    echo "       Expected a FLAT directory of <candidate_id>.{cif,pdb} models" >&2
    echo "       (stage 08 output, staged by scripts/unity/run_foldseek_candidates.sh)." >&2
    exit 3
fi
N_QUERIES="$(find -L "$QUERY_DIR" -maxdepth 1 -type f \
    \( -iname '*.cif' -o -iname '*.mmcif' -o -iname '*.bcif' \
       -o -iname '*.pdb' -o -iname '*.ent' \
       -o -iname '*.cif.gz' -o -iname '*.mmcif.gz' -o -iname '*.bcif.gz' \
       -o -iname '*.pdb.gz' -o -iname '*.ent.gz' \) 2>/dev/null | wc -l)"
if [ "$N_QUERIES" -eq 0 ]; then
    echo "ERROR: no query structures in $QUERY_DIR (0 .cif/.pdb/.mmcif/.ent files)." >&2
    echo "       Refusing to run foldseek over zero structures -- createdb would" >&2
    echo "       succeed on an empty DB and the whole search would exit 0 with an" >&2
    echo "       empty hits.tsv, which reads downstream as 'no structural evidence'." >&2
    echo "       Point \$1 at the flat <candidate_id>.<ext> query dir (run stage 08" >&2
    echo "       first, or use scripts/unity/run_foldseek_candidates.sh which stages it)." >&2
    exit 3
fi
echo "[foldseek-gpcrdb] query structures: ${N_QUERIES} in ${QUERY_DIR}" >&2

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
