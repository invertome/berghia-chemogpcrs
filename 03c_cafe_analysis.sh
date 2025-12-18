#!/bin/bash
# 03c_cafe_analysis.sh
# Purpose: Run CAFE5 analysis to detect significant gene family expansions/contractions.
# Inputs: OrthoFinder results, BUSCO species tree
# Outputs: CAFE5 results with significant LSE candidates
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=cafe_analysis
#SBATCH --output=${LOGS_DIR}/03c_cafe_analysis_%j.out
#SBATCH --error=${LOGS_DIR}/03c_cafe_analysis_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Initialize pipeline
init_pipeline

# Create output directories
CAFE_DIR="${RESULTS_DIR}/cafe"
mkdir -p "${CAFE_DIR}" "${LOGS_DIR}" || { log --level=ERROR "Cannot create directories"; exit 1; }

# Check for --force flag
FORCE_RERUN=false
if [ "$1" = "--force" ]; then
    FORCE_RERUN=true
    shift
fi

if [ "$FORCE_RERUN" = false ]; then
    skip_if_completed "cafe_analysis"
fi

log "Starting CAFE5 gene family evolution analysis"

# --- Check Dependencies ---
check_file "${RESULTS_DIR}/step_completed_orthology_clustering.txt"
check_file "${RESULTS_DIR}/step_completed_busco_species_tree.txt"

# Find OrthoFinder results
ORTHOFINDER_DIR="${RESULTS_DIR}/orthofinder"
ORTHOGROUPS_FILE=""

for og_file in \
    "${ORTHOFINDER_DIR}/Results_*/Orthogroups/Orthogroups.GeneCount.tsv" \
    "${ORTHOFINDER_DIR}/Orthogroups.GeneCount.tsv" \
    "${ORTHOFINDER_DIR}/Results_*/Orthogroups.GeneCount.tsv"; do
    if [ -f "$og_file" ]; then
        ORTHOGROUPS_FILE="$og_file"
        break
    fi
done

if [ -z "$ORTHOGROUPS_FILE" ]; then
    log --level=ERROR "OrthoFinder gene count file not found"
    exit 1
fi

# Find species tree
BUSCO_TREE=""
for tree_file in \
    "${SPECIES_TREE}" \
    "${RESULTS_DIR}/busco/busco_species_tree.tre" \
    "${RESULTS_DIR}/busco/species_tree.tre"; do
    if [ -f "$tree_file" ]; then
        BUSCO_TREE="$tree_file"
        break
    fi
done

if [ -z "$BUSCO_TREE" ]; then
    log --level=ERROR "BUSCO species tree not found"
    exit 1
fi

check_file "$ORTHOGROUPS_FILE" "$BUSCO_TREE"

log "Using OrthoFinder results: ${ORTHOGROUPS_FILE}"
log "Using species tree: ${BUSCO_TREE}"

# --- Step 1: Generate Ultrametric Tree ---
log "Step 1: Converting species tree to ultrametric"

ULTRAMETRIC_TREE="${CAFE_DIR}/species_tree_ultrametric.tre"

run_command "generate_ultrametric" \
    --inputs="${BUSCO_TREE}" \
    --outputs="${ULTRAMETRIC_TREE}" \
    ${RSCRIPT} "${SCRIPTS_DIR}/generate_ultrametric.R" \
    "${BUSCO_TREE}" "${ULTRAMETRIC_TREE}" "correlated" "1"

check_file "$ULTRAMETRIC_TREE"

# --- Step 2: Prepare CAFE Input Files ---
log "Step 2: Preparing CAFE5 input files"

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import os
import re
from ete3 import Tree

# Get paths
cafe_dir = os.environ.get('CAFE_DIR', 'results/cafe')
og_file = os.environ.get('ORTHOGROUPS_FILE', '')
tree_file = f"{cafe_dir}/species_tree_ultrametric.tre"

# Read OrthoFinder gene counts
df = pd.read_csv(og_file, sep='\t')

# Get species names from tree
tree = Tree(tree_file)
tree_species = set(leaf.name for leaf in tree)

# Find matching columns in gene count file
# OrthoFinder columns are typically: Orthogroup, Species1, Species2, ..., Total
species_cols = []
for col in df.columns:
    # Clean column name and check if it's in the tree
    clean_col = col.strip()
    if clean_col in tree_species:
        species_cols.append(col)
    else:
        # Try partial matching (e.g., taxid_berghia matches berghia)
        for sp in tree_species:
            if sp in clean_col or clean_col in sp:
                species_cols.append(col)
                break

print(f"Species in tree: {len(tree_species)}")
print(f"Matching columns: {len(species_cols)}")

if len(species_cols) < 2:
    print("WARNING: Not enough matching species found")
    print(f"Tree species: {tree_species}")
    print(f"DataFrame columns: {list(df.columns)}")

# Prepare CAFE input format
# CAFE format: Description\tID\tSpecies1\tSpecies2\t...
cafe_df = pd.DataFrame()
cafe_df['Description'] = df['Orthogroup'].apply(lambda x: f"OG_{x}")
cafe_df['ID'] = df['Orthogroup']

for col in species_cols:
    # Use the species name as it appears in the tree
    clean_name = col.strip()
    cafe_df[clean_name] = df[col].fillna(0).astype(int)

# Filter out families with too many genes (can cause CAFE to fail)
MAX_FAMILY_SIZE = 100
row_sums = cafe_df.iloc[:, 2:].sum(axis=1)
cafe_df = cafe_df[row_sums <= MAX_FAMILY_SIZE]

# Filter out families with all zeros
cafe_df = cafe_df[row_sums > 0]

# Write CAFE input file
output_file = f"{cafe_dir}/gene_families.tsv"
cafe_df.to_csv(output_file, sep='\t', index=False)

print(f"\nCAFE input file written: {output_file}")
print(f"  Gene families: {len(cafe_df)}")
print(f"  Species: {len(species_cols)}")
print(f"  (Families with >{MAX_FAMILY_SIZE} genes excluded)")

# Write species mapping
mapping_file = f"{cafe_dir}/species_mapping.txt"
with open(mapping_file, 'w') as f:
    for col in species_cols:
        f.write(f"{col}\n")

print(f"Species mapping written: {mapping_file}")
PYTHON_SCRIPT

check_file "${CAFE_DIR}/gene_families.tsv"

# --- Step 3: Run CAFE5 ---
log "Step 3: Running CAFE5 analysis"

# Check if CAFE5 is available
if ! command -v ${CAFE} &> /dev/null; then
    log --level=WARN "CAFE5 not found in PATH. Checking for local installation..."
    if [ -f "${BASE_DIR}/tools/cafe5" ]; then
        CAFE="${BASE_DIR}/tools/cafe5"
    else
        log --level=ERROR "CAFE5 not found. Please install: conda install -c bioconda cafe"
        exit 1
    fi
fi

# Prepare CAFE command
CAFE_OUTPUT="${CAFE_DIR}/cafe_results"
mkdir -p "${CAFE_OUTPUT}"

# Run CAFE5 with lambda search
if [ "${CAFE_LAMBDA_SEARCH}" = "true" ]; then
    log "Running CAFE5 with lambda search..."
    run_command "cafe5_analysis" \
        --inputs="${CAFE_DIR}/gene_families.tsv,${ULTRAMETRIC_TREE}" \
        --outputs="${CAFE_OUTPUT}" \
        --allow-fail \
        ${CAFE} \
        -i "${CAFE_DIR}/gene_families.tsv" \
        -t "${ULTRAMETRIC_TREE}" \
        -o "${CAFE_OUTPUT}" \
        -p
else
    log "Running CAFE5 with fixed lambda..."
    run_command "cafe5_analysis" \
        --inputs="${CAFE_DIR}/gene_families.tsv,${ULTRAMETRIC_TREE}" \
        --outputs="${CAFE_OUTPUT}" \
        --allow-fail \
        ${CAFE} \
        -i "${CAFE_DIR}/gene_families.tsv" \
        -t "${ULTRAMETRIC_TREE}" \
        -o "${CAFE_OUTPUT}"
fi

# --- Step 4: Parse CAFE5 Results ---
log "Step 4: Parsing CAFE5 results"

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import os
import json
from pathlib import Path

cafe_dir = os.environ.get('CAFE_DIR', 'results/cafe')
cafe_output = f"{cafe_dir}/cafe_results"
pvalue_threshold = float(os.environ.get('CAFE_PVALUE_THRESHOLD', 0.05))

# Find CAFE output files
results_found = False
significant_families = []

# Parse family-wise p-values
family_pvalues_file = Path(cafe_output) / "Base_family_results.txt"
if not family_pvalues_file.exists():
    # Try alternative naming
    for f in Path(cafe_output).glob("*family*.txt"):
        family_pvalues_file = f
        break

if family_pvalues_file.exists():
    print(f"Parsing: {family_pvalues_file}")
    results_found = True

    with open(family_pvalues_file, 'r') as f:
        header = None
        for line in f:
            parts = line.strip().split('\t')
            if header is None:
                header = parts
                continue

            if len(parts) >= 2:
                family_id = parts[0]
                try:
                    pvalue = float(parts[1]) if len(parts) > 1 else 1.0
                    if pvalue < pvalue_threshold:
                        significant_families.append({
                            'family_id': family_id,
                            'pvalue': pvalue,
                            'significant': True
                        })
                except ValueError:
                    continue

# Parse branch-specific results
branch_file = Path(cafe_output) / "Base_branch_probabilities.txt"
if not branch_file.exists():
    for f in Path(cafe_output).glob("*branch*.txt"):
        branch_file = f
        break

branch_results = []
if branch_file.exists():
    print(f"Parsing branch results: {branch_file}")
    with open(branch_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                branch_results.append({
                    'family_id': parts[0],
                    'branch': parts[1] if len(parts) > 1 else '',
                    'probability': float(parts[2]) if len(parts) > 2 else 1.0
                })

# Write results
output_file = f"{cafe_dir}/significant_expansions.tsv"
if significant_families:
    pd.DataFrame(significant_families).to_csv(output_file, sep='\t', index=False)
    print(f"\nSignificant families (p < {pvalue_threshold}): {len(significant_families)}")
    print(f"Results written to: {output_file}")
else:
    print("\nNo significant families found or CAFE5 output not parsed.")
    # Create empty file
    with open(output_file, 'w') as f:
        f.write("family_id\tpvalue\tsignificant\n")

# Write summary
summary = {
    'total_families_tested': len(branch_results) if branch_results else 0,
    'significant_families': len(significant_families),
    'pvalue_threshold': pvalue_threshold,
    'results_found': results_found
}

summary_file = f"{cafe_dir}/cafe_summary.json"
with open(summary_file, 'w') as f:
    json.dump(summary, f, indent=2)

print(f"Summary written to: {summary_file}")
PYTHON_SCRIPT

# Mark step as completed
create_checkpoint "cafe_analysis"

log "CAFE5 analysis completed"
log "Output files:"
log "  Ultrametric tree: ${ULTRAMETRIC_TREE}"
log "  Gene families input: ${CAFE_DIR}/gene_families.tsv"
log "  CAFE5 results: ${CAFE_DIR}/cafe_results/"
log "  Significant expansions: ${CAFE_DIR}/significant_expansions.tsv"
