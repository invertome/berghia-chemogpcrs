#!/usr/bin/env python3
# lse_refine.py
# Purpose: Classify orthogroups into multilevel LSE datasets with STRIDE-based duplication refinement and synteny validation.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from ete3 import NCBITaxa, Tree
import pandas as pd

og_file = sys.argv[1]
id_map_file = sys.argv[2]
species_tree_file = sys.argv[3]
gene_tree_file = sys.argv[4]
synteny_file = sys.argv[5]
output_dir = sys.argv[6]

ncbi = NCBITaxa()
id_map = pd.read_csv(id_map_file)
species_tree = Tree(species_tree_file)
gene_tree = Tree(gene_tree_file)
synteny_ids = pd.read_csv(synteny_file, names=['id'])['id'].tolist() if os.path.exists(synteny_file) else []

taxids = [line.split()[0][1:].split('_')[0] for line in open(og_file) if line.startswith('>')]
lineages = [set(ncbi.get_lineage(int(taxid))) for taxid in taxids]
common_lineage = set.intersection(*lineages)

if 54397 in common_lineage and all(t in [1263399, 54397] or ncbi.get_rank([t])[t] == 'species' for t in taxids):
    level = 'aeolids'
elif 13843 in common_lineage and all(t in [13843] or ncbi.get_rank([t])[t] == 'species' for t in taxids):
    level = 'nudibranchs'
elif 644 in common_lineage and all(t in [644] or ncbi.get_rank([t])[t] == 'species' for t in taxids):
    level = 'gastropods'
else:
    level = None

if level:
    # Check for duplications using gene tree reconciliation (STRIDE-like)
    taxid_count = pd.Series(taxids).value_counts()
    has_duplication = False
    for node in gene_tree.traverse():
        if not node.is_leaf():
            children_taxids = [leaf.name.split('_')[0] for leaf in node.get_leaves()]
            if len(set(children_taxids)) < len(children_taxids):  # Duplication detected
                has_duplication = True
                break
    if has_duplication and any(count > 1 for count in taxid_count):
        og_ids = [line.split()[0][1:] for line in open(og_file) if line.startswith('>')]
        synteny_overlap = any(id in synteny_ids for id in og_ids)
        if synteny_overlap:  # Validate with synteny
            base = og_file.split('/')[-1].replace('.fa', '')
            with open(f"{output_dir}/lse_{level}.txt", 'a') as f:
                f.write(f"{base}\n")
