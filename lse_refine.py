#!/usr/bin/env python3
# lse_refine.py
# Purpose: Classify orthogroups into multilevel LSE datasets using BUSCO gene tree reconciliation, ASTRAL species tree, and synteny validation.
# Inputs: Orthogroup FASTA ($1), ID map CSV ($2), ASTRAL species tree ($3), BUSCO gene tree ($4), synteny IDs ($5), output directory ($6)
# Outputs: LSE classification files (${output_dir}/lse_*.txt)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
import os
from ete3 import NCBITaxa, Tree
import pandas as pd

og_file = sys.argv[1]         # Orthogroup FASTA file
id_map_file = sys.argv[2]     # ID mapping CSV
species_tree_file = sys.argv[3]  # ASTRAL species tree file
gene_tree_file = sys.argv[4]  # Representative BUSCO gene tree file
synteny_file = sys.argv[5]    # Synteny IDs file
output_dir = sys.argv[6]      # Output directory

# Load taxonomic and tree data
ncbi = NCBITaxa()
id_map = pd.read_csv(id_map_file)
species_tree = Tree(species_tree_file)
gene_tree = Tree(gene_tree_file)
synteny_ids = pd.read_csv(synteny_file, names=['id'])['id'].tolist() if os.path.exists(synteny_file) else []

# Extract TaxIDs from orthogroup FASTA
taxid_strings = [line.split()[0][1:].split('_')[0] for line in open(og_file) if line.startswith('>')]

# Convert to integers, handling non-numeric taxids gracefully
taxids = []
for t in taxid_strings:
    try:
        taxids.append(int(t))
    except ValueError:
        # If taxid is not numeric (e.g., 'berghia'), try to look it up or skip
        print(f"Warning: Non-numeric taxid '{t}', skipping lineage lookup", file=sys.stderr)
        continue

if not taxids:
    print(f"Warning: No valid taxids found in {og_file}", file=sys.stderr)
    sys.exit(0)

# Get lineages for valid taxids
lineages = []
for taxid in taxids:
    try:
        lineage = ncbi.get_lineage(taxid)
        if lineage:
            lineages.append(set(lineage))
    except Exception as e:
        print(f"Warning: Could not get lineage for taxid {taxid}: {e}", file=sys.stderr)
        continue

if not lineages:
    print(f"Warning: No valid lineages found for {og_file}", file=sys.stderr)
    sys.exit(0)

common_lineage = set.intersection(*lineages)

# Determine taxonomic level for LSE classification
if 54397 in common_lineage and all(t in [1263399, 54397] or ncbi.get_rank([t])[t] == 'species' for t in taxids):
    level = 'aeolids'
elif 13843 in common_lineage and all(t in [13843] or ncbi.get_rank([t])[t] == 'species' for t in taxids):
    level = 'nudibranchs'
elif 644 in common_lineage and all(t in [644] or ncbi.get_rank([t])[t] == 'species' for t in taxids):
    level = 'gastropods'
else:
    level = None

# Process orthogroup if it matches a taxonomic level
if level:
    # Check for duplications using gene tree reconciliation
    taxid_count = pd.Series(taxids).value_counts()
    has_duplication = False
    for node in gene_tree.traverse():
        if not node.is_leaf():
            children_taxids = [leaf.name.split('_')[0] for leaf in node.get_leaves()]
            if len(set(children_taxids)) < len(children_taxids):  # Duplication detected
                has_duplication = True
                break
    if has_duplication and any(count > 1 for count in taxid_count):
        # Validate with synteny data
        og_ids = [line.split()[0][1:] for line in open(og_file) if line.startswith('>')]
        synteny_overlap = any(id in synteny_ids for id in og_ids)
        if synteny_overlap:
            base = og_file.split('/')[-1].replace('.fa', '')
            with open(f"{output_dir}/lse_{level}.txt", 'a') as f:
                f.write(f"{base}\n")
