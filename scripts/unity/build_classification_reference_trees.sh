#!/bin/bash
#SBATCH --job-name=class_ref_trees
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# build_classification_reference_trees.sh — Phase 3 Tasks 3.1 + 3.2 + 3.3
# of the non-chemoreceptor classification pipeline.
#
# Builds three IQ-TREE 3 reference trees from the curated non-chemoreceptor
# GPCR reference set (built in Phase 1):
#
#   1. backbone        — 54 reps (6 per coarse family × 9 families) for
#                        first-pass family-level placement of candidates.
#   2. aminergic       — all 165 reviewed aminergic refs for medium-
#                        granularity placement (5HT/dopamine/NE/HA/OA/Tyr).
#   3. peptide         — all 189 reviewed peptide refs for medium-
#                        granularity placement.
#
# Per-tree pipeline: PRODUCTION FILTER STACK
# (run_alignment_filter_stack from functions.sh) — PREQUAL + MAFFT
# canonical+4 variants+FAMSA ensemble + CLOAK consensus mask + TAPER
# residue-outlier mask, then ClipKit kpic-smart-gap, then IQ-TREE 3
# ModelFinder (LG/VT/WAG/JTT/Dayhoff/mtREV/cpREV) + UFBoot 1000 +
# SH-aLRT 1000 + TBE with deterministic seed.
#
# Earlier version of this script used a lighter MAFFT-only pipeline,
# arguing the Swiss-Prot reviewed inputs were already clean. Reverted
# 2026-05-07 for METHODOLOGICAL CONSISTENCY with stage 04 (the rest of
# the project's phylogenies all use the full stack) and because the
# subtrees span closely-related receptors (5HT vs dopamine etc.) where
# alignment uncertainty in ECL/ICL loops genuinely matters even on
# reviewed inputs. Defends against reviewer comments on inconsistency.
# Compute cost: ~60 extra min on Unity (one-time, cached).
#
# Submit:  sbatch scripts/unity/build_classification_reference_trees.sh

set -eo pipefail
mkdir -p logs

WORKDIR="/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05"
cd "$WORKDIR"

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate berghia-gpcr

# shellcheck disable=SC1091
source config.sh
# shellcheck disable=SC1091
source functions.sh

THREADS="${SLURM_CPUS_PER_TASK:-8}"

REF_FASTA="${WORKDIR}/references/non_chemo_gpcr/all_references.fasta"
REF_TSV="${WORKDIR}/references/non_chemo_gpcr/all_references.tsv"
OUT_DIR="${WORKDIR}/results/classification/trees"
mkdir -p "$OUT_DIR"

[[ -f "$REF_FASTA" ]] || { echo "ERROR: reference FASTA not found: $REF_FASTA" >&2; exit 1; }
[[ -f "$REF_TSV"   ]] || { echo "ERROR: reference TSV not found: $REF_TSV"   >&2; exit 1; }

echo "===================================================================="
echo "[$(date -u +%FT%TZ)] Building classification reference trees"
echo "===================================================================="
echo "Threads: $THREADS"
echo "Reference: $REF_FASTA ($(grep -c '^>' "$REF_FASTA") sequences)"
echo ""

# IQ-TREE 3 default arguments (consistent with main pipeline stage 04
# IQTREE_BOOT_FLAGS in config.sh / 04_phylogenetic_analysis.sh).
IQTREE_ARGS=(
    -m MFP
    -mset LG,VT,WAG,JTT,Dayhoff,mtREV,cpREV
    -B 1000
    -alrt 1000
    --tbe
    -seed 12345
    -T "$THREADS"
)

# ---- helper: build one tree end-to-end ---------------------------------
# args: <name> <input_fasta>
build_tree() {
    local name="$1"
    local input="$2"
    local n=$(grep -c '^>' "$input")

    echo "==== [tree:$name] $n sequences ===="
    local aln="$OUT_DIR/${name}.aln"
    local trimmed="$OUT_DIR/${name}.trimmed.aln"
    local stack_workdir="$OUT_DIR/_filter_stack_${name}"
    local prefix="$OUT_DIR/${name}"

    if [[ ! -s "$aln" ]]; then
        echo "  [tree:$name] Filter stack (PREQUAL → MAFFT ensemble → CLOAK → TAPER)..."
        run_alignment_filter_stack \
            "$input" "$aln" "$stack_workdir" "$name" "$THREADS" \
            || { echo "  [tree:$name] filter stack FAILED" >&2; return 1; }
    else
        echo "  [tree:$name] MSA already present: $aln"
    fi

    if [[ ! -s "$trimmed" ]]; then
        echo "  [tree:$name] ClipKit kpic-smart-gap..."
        clipkit "$aln" -m kpic-smart-gap -o "$trimmed" \
            > "$OUT_DIR/${name}.clipkit.log" 2>&1
    else
        echo "  [tree:$name] Trimmed alignment already present"
    fi

    if [[ ! -s "${prefix}.treefile" ]]; then
        echo "  [tree:$name] IQ-TREE 3..."
        iqtree3 -s "$trimmed" --prefix "$prefix" "${IQTREE_ARGS[@]}" \
            > "$OUT_DIR/${name}.iqtree.log" 2>&1
    else
        echo "  [tree:$name] Tree already present: ${prefix}.treefile"
    fi

    if [[ -s "${prefix}.treefile" ]]; then
        echo "  [tree:$name] OK -> ${prefix}.treefile"
    else
        echo "  [tree:$name] FAILED"
        return 1
    fi
}

# ---- 1. Backbone tree ---------------------------------------------------
echo ""
echo "[1/3] Backbone tree (~54 reps across 9 coarse families)"
BB_FASTA="$OUT_DIR/backbone.fasta"
BB_TSV="$OUT_DIR/backbone.tsv"
if [[ ! -s "$BB_FASTA" ]]; then
    python3 scripts/select_backbone_reps.py \
        --reference-fasta "$REF_FASTA" \
        --reference-tsv "$REF_TSV" \
        --output-fasta "$BB_FASTA" \
        --output-tsv "$BB_TSV" \
        --quota-per-family 6
fi
build_tree "backbone" "$BB_FASTA"

# ---- 2. Aminergic subtree -----------------------------------------------
echo ""
echo "[2/3] Aminergic subtree (all reviewed aminergic refs)"
AM_FASTA="$OUT_DIR/aminergic_subtree.fasta"
AM_TSV="$OUT_DIR/aminergic_subtree.tsv"
if [[ ! -s "$AM_FASTA" ]]; then
    # Extract aminergic refs from the full set, write per-leaf TSV.
    python3 -c "
import csv, sys
keep = set()
rows = []
with open('${REF_TSV}') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        if row['family'] == 'aminergic':
            keep.add(row['accession'])
            rows.append(row)
with open('${AM_TSV}', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['accession', 'family', 'subfamily', 'species', 'gene'])
    for row in rows:
        w.writerow([row['accession'], row['family'], row['subfamily'],
                    row['species'], row['gene']])
with open('${REF_FASTA}') as src, open('${AM_FASTA}', 'w') as dst:
    keep_block = False
    for line in src:
        if line.startswith('>'):
            acc = line[1:].split('|', 1)[0].strip()
            keep_block = acc in keep
        if keep_block:
            dst.write(line)
print(f'aminergic subtree: {len(keep)} refs', file=sys.stderr)
"
fi
build_tree "aminergic_subtree" "$AM_FASTA"

# ---- 3. Peptide subtree -------------------------------------------------
echo ""
echo "[3/3] Peptide subtree (all reviewed peptide refs)"
PEP_FASTA="$OUT_DIR/peptide_subtree.fasta"
PEP_TSV="$OUT_DIR/peptide_subtree.tsv"
if [[ ! -s "$PEP_FASTA" ]]; then
    python3 -c "
import csv, sys
keep = set()
rows = []
with open('${REF_TSV}') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        if row['family'] == 'peptide':
            keep.add(row['accession'])
            rows.append(row)
with open('${PEP_TSV}', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['accession', 'family', 'subfamily', 'species', 'gene'])
    for row in rows:
        w.writerow([row['accession'], row['family'], row['subfamily'],
                    row['species'], row['gene']])
with open('${REF_FASTA}') as src, open('${PEP_FASTA}', 'w') as dst:
    keep_block = False
    for line in src:
        if line.startswith('>'):
            acc = line[1:].split('|', 1)[0].strip()
            keep_block = acc in keep
        if keep_block:
            dst.write(line)
print(f'peptide subtree: {len(keep)} refs', file=sys.stderr)
"
fi
build_tree "peptide_subtree" "$PEP_FASTA"

# ---- Manifest -----------------------------------------------------------
echo ""
echo "===================================================================="
echo "[$(date -u +%FT%TZ)] Tree manifest"
echo "===================================================================="
MANIFEST="$OUT_DIR/manifest.tsv"
{
    echo -e "tree_name\tfasta\ttreefile\tleaf_tsv\tn_leaves\tgit_sha"
    for name in backbone aminergic_subtree peptide_subtree; do
        fasta="$OUT_DIR/${name}.fasta"
        tre="$OUT_DIR/${name}.treefile"
        tsv="$OUT_DIR/${name}.tsv"
        n=$(grep -c '^>' "$fasta" 2>/dev/null || echo 0)
        sha=$(git -C "$WORKDIR" rev-parse --short HEAD 2>/dev/null || echo unknown)
        echo -e "${name}\t${fasta}\t${tre}\t${tsv}\t${n}\t${sha}"
    done
} > "$MANIFEST"
cat "$MANIFEST"

echo ""
echo "[$(date -u +%FT%TZ)] DONE — 3 reference trees in $OUT_DIR"
