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

log "Starting NOTUNG gene tree reconciliation"

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
GENE_TREES=$(find "$PHYLO_DIR" -name "*.treefile" -o -name "*_tree.tre" 2>/dev/null | head -100)

if [ -z "$GENE_TREES" ]; then
    log --level=WARN "No gene trees found in ${PHYLO_DIR}"
    log "Creating checkpoint and exiting"
    create_checkpoint "notung_reconciliation"
    exit 0
fi

TREE_COUNT=$(echo "$GENE_TREES" | wc -l)
log "Found ${TREE_COUNT} gene trees to reconcile"

# --- Check for NOTUNG ---
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

# Find gene trees
gene_tree_files = list(Path(phylo_dir).glob("*.treefile")) + list(Path(phylo_dir).glob("*_tree.tre"))

results = []
duplication_events = defaultdict(list)
loss_events = defaultdict(list)

for tree_file in gene_tree_files[:100]:  # Limit to 100 trees
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
    'total_duplications': df['duplications'].sum(),
    'total_losses': df['losses'].sum()
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
