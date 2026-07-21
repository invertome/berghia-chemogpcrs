#!/bin/bash
# 03d_notung_reconciliation.sh
# Purpose: Reconcile gene trees with species tree using NOTUNG to infer duplication/loss events.
# Inputs: Gene trees from phylogenetic analysis, species tree from BUSCO
# Outputs: Reconciled trees with duplication/loss annotations
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=notung_reconcile
#SBATCH --output=${LOGS_DIR}/03d_notung_reconciliation_%j.out
#SBATCH --error=${LOGS_DIR}/03d_notung_reconciliation_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Initialize pipeline
init_pipeline

# Create output directories
NOTUNG_DIR="${RESULTS_DIR}/notung"
mkdir -p "${NOTUNG_DIR}" "${LOGS_DIR}" || { log --level=ERROR "Cannot create directories"; exit 1; }

# Check for --force flag
FORCE_RERUN=false
if [ "$1" = "--force" ]; then
    FORCE_RERUN=true
    shift
fi

if [ "$FORCE_RERUN" = false ]; then
    skip_if_completed "notung_reconciliation"
fi

log "Starting gene tree reconciliation"

# Bead -30g: RECONCILIATION_BACKEND=generax (default) uses GeneRax
# (Morel 2020 MBE 37:2763) which estimates joint phylogenetic +
# reconciliation likelihoods. RECONCILIATION_BACKEND=notung is the legacy
# parsimony-style fallback. We auto-fall-back to Notung when GeneRax
# is not installed.
RECONCILIATION_BACKEND="${RECONCILIATION_BACKEND:-generax}"
GENERAX_BIN="${GENERAX:-generax}"
if [ "$RECONCILIATION_BACKEND" = "generax" ]; then
    if ! command -v "$GENERAX_BIN" &>/dev/null; then
        log --level=WARN "GeneRax binary '$GENERAX_BIN' not found; falling back to Notung."
        log --level=WARN "Install GeneRax: see docs/installation_hpc.md"
        RECONCILIATION_BACKEND="notung"
    fi
fi
log "Reconciliation backend: $RECONCILIATION_BACKEND"

# --- Check Dependencies ---
# Step 04 creates step_completed_04.txt
check_file "${RESULTS_DIR}/step_completed_04.txt"
check_file "${RESULTS_DIR}/step_completed_busco_species_tree.txt"

# Find species tree
SPECIES_TREE_FILE=""
for tree_file in \
    "${SPECIES_TREE}" \
    "${RESULTS_DIR}/busco/busco_species_tree.tre" \
    "${RESULTS_DIR}/cafe/species_tree_ultrametric.tre"; do
    if [ -f "$tree_file" ]; then
        SPECIES_TREE_FILE="$tree_file"
        break
    fi
done

if [ -z "$SPECIES_TREE_FILE" ]; then
    log --level=ERROR "Species tree not found"
    exit 1
fi

log "Using species tree: ${SPECIES_TREE_FILE}"

# Find gene trees
PHYLO_DIR="${RESULTS_DIR}/phylogenies"
# Bead F6: `-name A -o -name B` is now explicitly grouped, and TreeShrink's
# rollback copies are excluded. run_treeshrink.sh copies every tree it edits to
# `<tree>.original.treefile`, which matches `*.treefile`; counting it as a
# separate gene tree reconciles the same orthogroup twice and double-counts its
# duplications. The Python fallback below applies the identical definition so
# the two backends agree on what a gene tree is.
GENE_TREES=$(find "$PHYLO_DIR" \( -name "*.treefile" -o -name "*_tree.tre" \) \
    ! -name "*.original.treefile" 2>/dev/null)

if [ -z "$GENE_TREES" ]; then
    log --level=WARN "No gene trees found in ${PHYLO_DIR}"
    log "Creating checkpoint and exiting"
    create_checkpoint "notung_reconciliation"
    exit 0
fi

# Discovery and the processing cap are kept SEPARATE. They used to be one
# statement (`find ... | head -100`), so TREE_COUNT counted the truncated list
# and the stage logged "Found 100 gene trees to reconcile" whether 100 trees
# existed or 4,000 did -- a run covering 2.5% of the orthogroups read exactly
# like a complete one, in the log and in events_summary.json alike.
#
# The cap is deliberately UNCHANGED: raising or removing it is a resource
# decision for the user, not one this stage should make silently. What changes
# is that the denominator survives, so that decision can be made on evidence.
TREE_TOTAL=$(printf '%s\n' "$GENE_TREES" | grep -c .)
RECONCILE_TREE_CAP=100
GENE_TREES=$(printf '%s\n' "$GENE_TREES" | head -n "$RECONCILE_TREE_CAP")
TREE_COUNT=$(printf '%s\n' "$GENE_TREES" | grep -c .)

if [ "$TREE_COUNT" -lt "$TREE_TOTAL" ]; then
    log --level=WARN "TRUNCATED: reconciling ${TREE_COUNT} of ${TREE_TOTAL} gene trees found (cap RECONCILE_TREE_CAP=${RECONCILE_TREE_CAP}); $((TREE_TOTAL - TREE_COUNT)) will NOT be reconciled."
    log --level=WARN "Duplication/loss counts from this run are LOWER BOUNDS over a ${TREE_COUNT}/${TREE_TOTAL} subset, not totals. Raise the cap in ${BASH_SOURCE[0]##*/} to cover the full set."
else
    log "Found ${TREE_COUNT} gene trees to reconcile (complete set; cap ${RECONCILE_TREE_CAP} not reached)"
fi

# --- GeneRax branch (preferred when binary available) ---
if [ "$RECONCILIATION_BACKEND" = "generax" ]; then
    GENERAX_OUT="${RESULTS_DIR}/reconciliation/generax"
    mkdir -p "$GENERAX_OUT"
    # Build a minimal families.txt: one entry per gene tree.
    FAMILIES_TXT="$GENERAX_OUT/families.txt"
    {
        echo "[FAMILIES]"
        for tree in $GENE_TREES; do
            base=$(basename "$tree" .treefile)
            base="${base%.tre}"
            # Bead F6: locate the alignment relative to the TREE'S OWN
            # directory. The previous lookup used the pre-refactor FLAT paths
            # (phylogenies/protein/<base>_trimmed.fa), which nothing writes
            # since stage 04's per-class refactor -- so `[ -z "$aln" ] &&
            # continue` dropped EVERY family and families.txt was left as a
            # bare [FAMILIES] header, i.e. GeneRax reconciled nothing.
            #
            # Anchoring on dirname covers all three layouts stage 04 emits,
            # without stage 03d having to re-derive the orthogroup's class:
            #   per-OG        class_<C>/<base>/<base>.treefile + <base>_trimmed.fa
            #   per-class     class_<C>/class_<C>.treefile     + class_<C>_trimmed.fa
            #   back-compat   protein/all_berghia_refs.treefile + ..._trimmed.fa
            # The <base>/trimmed.fa candidates cover the LSE trees, whose
            # IQ-TREE --prefix puts the treefile one level ABOVE its alignment
            # directory (04:436, 04:452).
            tree_dir=$(dirname "$tree")
            aln=""
            for cand in "${tree_dir}/${base}_trimmed.fa" \
                        "${tree_dir}/${base}_aligned.fa" \
                        "${tree_dir}/${base}/trimmed.fa" \
                        "${tree_dir}/${base}/aligned.fa"; do
                [ -f "$cand" ] && { aln="$cand"; break; }
            done
            [ -z "$aln" ] && continue
            echo "- $base"
            echo "starting_gene_tree = $tree"
            echo "alignment = $aln"
            echo "subst_model = LG+G"
        done
    } > "$FAMILIES_TXT"
    # Bead F6: report the number of families ACTUALLY written, not the number
    # of trees found. The old message reported $TREE_COUNT unconditionally, so
    # it cheerfully announced "Submitting GeneRax with 100 families" while
    # families.txt held nothing but its header -- which is how a reconciliation
    # step that reconciled zero families went unnoticed.
    FAMILY_COUNT=$(grep -c '^- ' "$FAMILIES_TXT" 2>/dev/null || true)
    FAMILY_COUNT="${FAMILY_COUNT:-0}"
    if [ "$FAMILY_COUNT" -eq 0 ]; then
        log --level=ERROR "families.txt is empty: none of the ${TREE_COUNT} gene trees had a matching alignment beside them. GeneRax would reconcile nothing; falling back to Notung."
        RECONCILIATION_BACKEND="notung"
    else
        log "Submitting GeneRax with ${FAMILY_COUNT} families (from ${TREE_COUNT} gene trees)"
    fi
fi
if [ "$RECONCILIATION_BACKEND" = "generax" ]; then
    bash "${SCRIPTS_DIR}/hpc/run_generax.sh" "$FAMILIES_TXT" \
        2>> "${LOGS_DIR}/generax.err" \
        || { log --level=WARN "GeneRax failed; falling back to Notung."; \
             RECONCILIATION_BACKEND="notung"; }
    if [ "$RECONCILIATION_BACKEND" = "generax" ]; then
        # Success — emit completion flag and exit
        touch "${RESULTS_DIR}/step_completed_notung_reconciliation.txt"
        create_checkpoint "notung_reconciliation"
        log "GeneRax reconciliation complete -> $GENERAX_OUT"
        exit 0
    fi
fi

# --- Check for NOTUNG (legacy fallback) ---
NOTUNG_JAR=""
if [ -f "${BASE_DIR}/tools/Notung-2.9.jar" ]; then
    NOTUNG_JAR="${BASE_DIR}/tools/Notung-2.9.jar"
elif [ -f "${BASE_DIR}/tools/Notung.jar" ]; then
    NOTUNG_JAR="${BASE_DIR}/tools/Notung.jar"
fi

if [ -z "$NOTUNG_JAR" ] || ! command -v java &> /dev/null; then
    log --level=WARN "NOTUNG not found or Java not available"
    log "To install: download from http://www.cs.cmu.edu/~durand/Notung/"
    log "Place Notung-2.9.jar in ${BASE_DIR}/tools/"
    log "Falling back to simple reconciliation analysis"

    # Fall back to Python-based reconciliation
    USE_PYTHON_FALLBACK=true
else
    USE_PYTHON_FALLBACK=false
    log "Using NOTUNG: ${NOTUNG_JAR}"
fi

# --- Step 1: Prepare Species Tree for NOTUNG ---
log "Step 1: Preparing species tree"

# NOTUNG requires specific format - convert if needed
NOTUNG_SPECIES_TREE="${NOTUNG_DIR}/species_tree_notung.tre"

python3 << PYTHON_SCRIPT
from ete3 import Tree
import os

species_tree_file = "${SPECIES_TREE_FILE}"
output_file = "${NOTUNG_SPECIES_TREE}"

# Read tree
tree = Tree(species_tree_file)

# NOTUNG requires rooted tree with internal node labels
# Label internal nodes if not already labeled
node_counter = 1
for node in tree.traverse():
    if not node.is_leaf() and not node.name:
        node.name = f"N{node_counter}"
        node_counter += 1

# Write in Newick format
tree.write(outfile=output_file, format=1)
print(f"Prepared species tree: {output_file}")
PYTHON_SCRIPT

check_file "$NOTUNG_SPECIES_TREE"

# --- Step 2: Run Reconciliation ---
log "Step 2: Running gene tree reconciliation"

export NOTUNG_DIR NOTUNG_SPECIES_TREE PHYLO_DIR RESULTS_DIR

if [ "$USE_PYTHON_FALLBACK" = true ]; then
    # Python-based reconciliation (simpler but functional)
    log "Using Python-based reconciliation analysis"

    python3 << 'PYTHON_SCRIPT'
import os
import sys
import json
from pathlib import Path
from collections import defaultdict
from ete3 import Tree

notung_dir = os.environ.get('NOTUNG_DIR', 'results/notung')
species_tree_file = os.environ.get('NOTUNG_SPECIES_TREE', '')
phylo_dir = os.environ.get('PHYLO_DIR', 'results/phylogenies')

# Load species tree
species_tree = Tree(species_tree_file)
species_names = set(leaf.name for leaf in species_tree)

print(f"Species in tree: {len(species_names)}", file=sys.stderr)

# Find gene trees.
# Bead F6: this globbed NON-recursively while the bash `find` above recurses,
# so after stage 04's per-class refactor moved every per-OG tree two levels
# deeper the two reconciliation backends disagreed about what a gene tree is --
# and this one always saw zero. rglob + the same .original.treefile exclusion
# (TreeShrink rollback copies) keeps the definitions identical.
gene_tree_files = [p for p in list(Path(phylo_dir).rglob("*.treefile")) + list(Path(phylo_dir).rglob("*_tree.tre")) if not p.name.endswith(".original.treefile")]

results = []
duplication_events = defaultdict(list)
loss_events = defaultdict(list)

# This fallback applies its OWN cap, independent of the bash branch's. Keep the
# denominator so events_summary.json can say what fraction was covered: without
# it, 'total_trees_analyzed' reads as a population when it is a capped sample.
RECONCILE_TREE_CAP = 100
gene_trees_discovered = len(gene_tree_files)
selected_trees = gene_tree_files[:RECONCILE_TREE_CAP]
reconciliation_truncated = gene_trees_discovered > len(selected_trees)
if reconciliation_truncated:
    print(f"WARN: TRUNCATED -- reconciling {len(selected_trees)} of "
          f"{gene_trees_discovered} gene trees (cap {RECONCILE_TREE_CAP}). "
          f"Duplication/loss counts below are LOWER BOUNDS, not totals.",
          file=sys.stderr)

for tree_file in selected_trees:
    try:
        gene_tree = Tree(str(tree_file), format=1)
        tree_name = tree_file.stem

        # Get gene names and map to species
        gene_to_species = {}
        for leaf in gene_tree:
            # Try to extract species from gene name (format: species_geneid or taxid_species_geneid)
            parts = leaf.name.split('_')
            species = None
            for sp in species_names:
                if sp in leaf.name or any(sp in p for p in parts):
                    species = sp
                    break
            gene_to_species[leaf.name] = species

        # Simple duplication detection: find nodes where children have same species
        duplications = 0
        losses = 0

        for node in gene_tree.traverse():
            if node.is_leaf():
                continue

            # Get species in each child clade
            child_species = []
            for child in node.children:
                child_sp = set()
                for leaf in child:
                    sp = gene_to_species.get(leaf.name)
                    if sp:
                        child_sp.add(sp)
                child_species.append(child_sp)

            # If same species appears in multiple children, likely duplication
            if len(child_species) >= 2:
                overlap = child_species[0] & child_species[1]
                if overlap:
                    duplications += 1
                    for sp in overlap:
                        duplication_events[sp].append(tree_name)

        # Estimate losses (species present in species tree but not gene tree)
        gene_species = set(gene_to_species.values()) - {None}
        missing_species = species_names - gene_species
        losses = len(missing_species)
        for sp in missing_species:
            loss_events[sp].append(tree_name)

        results.append({
            'tree': tree_name,
            'n_genes': len(gene_tree),
            'n_species': len(gene_species),
            'duplications': duplications,
            'losses': losses,
            'missing_species': ','.join(missing_species) if missing_species else ''
        })

    except Exception as e:
        print(f"Warning: Could not process {tree_file}: {e}", file=sys.stderr)
        continue

# Write results
import pandas as pd

df = pd.DataFrame(results)
df.to_csv(f"{notung_dir}/reconciliation_summary.tsv", sep='\t', index=False)

# Write duplication/loss events per species
events_summary = {
    'duplication_events': {sp: len(trees) for sp, trees in duplication_events.items()},
    'loss_events': {sp: len(trees) for sp, trees in loss_events.items()},
    'total_trees_analyzed': len(results),
    'total_duplications': int(df['duplications'].sum()),
    'total_losses': int(df['losses'].sum()),
    # Coverage denominator. The three 'total_*' keys above are sums over the
    # trees ACTUALLY reconciled; when reconciliation_truncated is true they are
    # lower bounds over a capped subset, and a consumer reading the JSON alone
    # has no other way to tell that apart from a complete census.
    'gene_trees_discovered': gene_trees_discovered,
    'gene_trees_reconciled': len(selected_trees),
    'reconciliation_tree_cap': RECONCILE_TREE_CAP,
    'reconciliation_truncated': reconciliation_truncated,
}

with open(f"{notung_dir}/events_summary.json", 'w') as f:
    json.dump(events_summary, f, indent=2)

# Print summary
print(f"\n=== Reconciliation Summary ===", file=sys.stderr)
print(f"Trees analyzed: {len(results)}", file=sys.stderr)
print(f"Total duplications detected: {df['duplications'].sum()}", file=sys.stderr)
print(f"Total losses estimated: {df['losses'].sum()}", file=sys.stderr)

# Top species by duplications
print("\nTop species by duplication events:", file=sys.stderr)
for sp, count in sorted(events_summary['duplication_events'].items(), key=lambda x: -x[1])[:5]:
    print(f"  {sp}: {count}", file=sys.stderr)
PYTHON_SCRIPT

else
    # Use NOTUNG for reconciliation
    log "Running NOTUNG reconciliation"

    RECONCILED_DIR="${NOTUNG_DIR}/reconciled"
    mkdir -p "$RECONCILED_DIR"

    # Process each gene tree
    for gene_tree in $GENE_TREES; do
        tree_name=$(basename "$gene_tree" .treefile)
        tree_name=${tree_name%_tree.tre}

        log "Processing: ${tree_name}"

        # Run NOTUNG
        java -jar "$NOTUNG_JAR" \
            -g "$gene_tree" \
            -s "$NOTUNG_SPECIES_TREE" \
            --reconcile \
            --rearrange \
            --threshold 90 \
            --savepng \
            --outputdir "$RECONCILED_DIR" \
            2>> "${LOGS_DIR}/notung_${tree_name}.err" || {
                log --level=WARN "NOTUNG failed for ${tree_name}"
                continue
            }
    done

    # Parse NOTUNG outputs
    python3 << 'PYTHON_SCRIPT'
import os
import json
from pathlib import Path
from collections import defaultdict

notung_dir = os.environ.get('NOTUNG_DIR', 'results/notung')
reconciled_dir = f"{notung_dir}/reconciled"

results = []
duplication_events = defaultdict(int)
loss_events = defaultdict(int)

# Parse NOTUNG output files
for output_file in Path(reconciled_dir).glob("*.reconciled"):
    try:
        with open(output_file) as f:
            content = f.read()

        tree_name = output_file.stem

        # Count duplication and loss events from NOTUNG output
        duplications = content.count('[D]') + content.count('[Dup]')
        losses = content.count('[L]') + content.count('[Loss]')

        results.append({
            'tree': tree_name,
            'duplications': duplications,
            'losses': losses
        })

    except Exception as e:
        print(f"Warning: Could not parse {output_file}: {e}")

# Write summary
import pandas as pd

if results:
    df = pd.DataFrame(results)
    df.to_csv(f"{notung_dir}/reconciliation_summary.tsv", sep='\t', index=False)

    summary = {
        'total_trees_analyzed': len(results),
        'total_duplications': df['duplications'].sum(),
        'total_losses': df['losses'].sum(),
        'mean_duplications_per_tree': df['duplications'].mean(),
        'mean_losses_per_tree': df['losses'].mean()
    }

    with open(f"{notung_dir}/events_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"Reconciliation complete: {len(results)} trees processed")
    print(f"Total duplications: {summary['total_duplications']}")
    print(f"Total losses: {summary['total_losses']}")
PYTHON_SCRIPT
fi

# Export variables for Python
export NOTUNG_DIR NOTUNG_SPECIES_TREE PHYLO_DIR

# --- Step 3: Generate Candidate Annotations ---
log "Step 3: Generating candidate duplication annotations"

python3 << 'PYTHON_SCRIPT'
import os
import json
import pandas as pd
from pathlib import Path

notung_dir = os.environ.get('NOTUNG_DIR', 'results/notung')
results_dir = os.environ.get('RESULTS_DIR', 'results')

# Load reconciliation summary
summary_file = f"{notung_dir}/reconciliation_summary.tsv"
if os.path.exists(summary_file):
    df = pd.read_csv(summary_file, sep='\t')

    # Create annotation file for candidates
    # Map trees back to orthogroups/candidates
    annotations = []

    for _, row in df.iterrows():
        tree_name = row['tree']
        duplications = row.get('duplications', 0)
        losses = row.get('losses', 0)

        annotations.append({
            'orthogroup': tree_name,
            'duplications': duplications,
            'losses': losses,
            'dup_loss_ratio': duplications / max(losses, 1),
            'expansion_signal': duplications > losses
        })

    ann_df = pd.DataFrame(annotations)
    ann_df.to_csv(f"{notung_dir}/candidate_annotations.tsv", sep='\t', index=False)

    # Identify orthogroups with strong expansion signal
    expanded = ann_df[ann_df['expansion_signal'] & (ann_df['duplications'] > 2)]
    expanded['orthogroup'].to_csv(
        f"{notung_dir}/expanded_orthogroups.txt",
        index=False, header=False
    )

    print(f"Candidate annotations written: {len(ann_df)} orthogroups")
    print(f"Expanded orthogroups identified: {len(expanded)}")
PYTHON_SCRIPT

# Mark step as completed
create_checkpoint "notung_reconciliation"

log "NOTUNG reconciliation completed"
log "Output files:"
log "  Reconciliation summary: ${NOTUNG_DIR}/reconciliation_summary.tsv"
log "  Events summary: ${NOTUNG_DIR}/events_summary.json"
log "  Candidate annotations: ${NOTUNG_DIR}/candidate_annotations.tsv"
log "  Expanded orthogroups: ${NOTUNG_DIR}/expanded_orthogroups.txt"
