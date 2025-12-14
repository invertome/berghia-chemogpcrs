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

# Configurable taxonomy IDs from environment (with NCBI defaults)
# These can be set in config.sh: export LSE_AEOLID_TAXID=54397
LSE_AEOLID_TAXID = int(os.getenv('LSE_AEOLID_TAXID', 54397))      # Aeolidida
LSE_NUDIBRANCH_TAXID = int(os.getenv('LSE_NUDIBRANCH_TAXID', 13843))  # Nudibranchia
LSE_GASTROPOD_TAXID = int(os.getenv('LSE_GASTROPOD_TAXID', 644))   # Gastropoda

# Load taxonomic and tree data
ncbi = NCBITaxa()
id_map = pd.read_csv(id_map_file)
species_tree = Tree(species_tree_file, format=1)  # format=1 for trees with internal node names
gene_tree = Tree(gene_tree_file, format=1)
synteny_ids = pd.read_csv(synteny_file, names=['id'])['id'].tolist() if os.path.exists(synteny_file) else []

# Lineage cache for performance
_lineage_cache = {}
_rank_cache = {}

def get_cached_lineage(taxid):
    """Get lineage with caching to avoid repeated NCBI queries."""
    if taxid not in _lineage_cache:
        try:
            _lineage_cache[taxid] = ncbi.get_lineage(taxid)
        except Exception:
            _lineage_cache[taxid] = None
    return _lineage_cache[taxid]

def get_cached_rank(taxid):
    """Get rank with caching."""
    if taxid not in _rank_cache:
        try:
            ranks = ncbi.get_rank([taxid])
            _rank_cache[taxid] = ranks.get(taxid, 'unknown')
        except Exception:
            _rank_cache[taxid] = 'unknown'
    return _rank_cache[taxid]

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

# Get lineages for valid taxids (using cache)
lineages = []
for taxid in taxids:
    lineage = get_cached_lineage(taxid)
    if lineage:
        lineages.append(set(lineage))
    else:
        print(f"Warning: Could not get lineage for taxid {taxid}", file=sys.stderr)

if not lineages:
    print(f"Warning: No valid lineages found for {og_file}", file=sys.stderr)
    sys.exit(0)

common_lineage = set.intersection(*lineages)

# Determine taxonomic level for LSE classification (using configurable taxids)
if LSE_AEOLID_TAXID in common_lineage and all(get_cached_rank(t) == 'species' or t == LSE_AEOLID_TAXID for t in taxids):
    level = 'aeolids'
elif LSE_NUDIBRANCH_TAXID in common_lineage and all(get_cached_rank(t) == 'species' or t == LSE_NUDIBRANCH_TAXID for t in taxids):
    level = 'nudibranchs'
elif LSE_GASTROPOD_TAXID in common_lineage and all(get_cached_rank(t) == 'species' or t == LSE_GASTROPOD_TAXID for t in taxids):
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
