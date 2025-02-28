#!/usr/bin/env python3
# fetch_ligands.py
# Purpose: Fetch high-quality GPCR chemoreceptor references from GPCRdb with invertebrate prioritization and estimate putative ligands for Berghia GPCRs.
# Inputs: Reference directory ($1), output CSV ($2), TM-align scores ($3), UPGMA tree ($4), FoldTree tree ($5), putative ligands CSV ($6)
# Outputs: Reference ligands CSV, putative ligands CSV for Berghia GPCRs
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
import requests
import pandas as pd
from Bio import PDB
from ete3 import Tree

# Command-line arguments
ref_dir = sys.argv[1]  # Directory for reference PDB files
output_csv = sys.argv[2]  # Output CSV for reference ligands
tmalign_scores = sys.argv[3] if len(sys.argv) > 3 else None  # TM-align scores CSV
upgma_tree = sys.argv[4] if len(sys.argv) > 4 else None  # UPGMA tree file
foldtree_tree = sys.argv[5] if len(sys.argv) > 5 else None  # FoldTree tree file
putative_ligands_csv = sys.argv[6] if len(sys.argv) > 6 else None  # Output CSV for putative ligands

# Fetch GPCRs from GPCRdb, prioritizing invertebrates and chemoreceptors
response = requests.get("https://gpcrdb.org/services/receptorlist/?family=Class_A")
if response.status_code == 200:
    gpcrs = response.json()
    ref_ids = [g['uniprot'] for g in gpcrs if ('OR' in g['uniprot'] or 'TAS' in g['uniprot'] or 'VR' in g['uniprot'] or 'peptide' in g['name'].lower()) and ('Aplysia' in g['species'] or 'Lottia' in g['species'] or 'invertebrate' in g['species'])]
    ref_ids += ['P30953', 'Q8NGD0', 'OR1A1', 'OR2J2', 'TAS2R38']  # Additional known chemoreceptors
    ref_ids = list(set(ref_ids))[:50]  # Limit to 50 unique entries
else:
    ref_ids = ['OR1A1', 'OR2J2', 'TAS2R38', 'P30953']  # Fallback list

# Download PDB structures and extract ligands
pdb_list = PDB.PDBList()
parser = PDB.PDBParser()
ref_data = []

for ref_id in ref_ids:
    pdb_response = requests.get(f"https://gpcrdb.org/services/structure/{ref_id}")
    if pdb_response.status_code == 200 and pdb_response.json():
        pdb_id = pdb_response.json()[0]['pdb_code']
        pdb_file = pdb_list.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=ref_dir)
        structure = parser.get_structure(pdb_id, pdb_file)
        ligands = [res.resname for res in structure[0].get_residues() if res.id[0].startswith('H_')]
        if ligands:
            ref_data.append({'ref_id': ref_id, 'pdb_id': pdb_id, 'ligands': ','.join(ligands)})

# Save reference data to CSV
pd.DataFrame(ref_data).to_csv(output_csv, index=False)

# Estimate putative ligands if additional inputs are provided
if tmalign_scores and upgma_tree and foldtree_tree and putative_ligands_csv:
    scores = pd.read_csv(tmalign_scores, names=['file', 'tm_score'], sep=',')
    scores['berghia_id'] = scores['file'].apply(lambda x: x.split('tmalign_')[1].split('_')[0])
    scores['ref_id'] = scores['file'].apply(lambda x: '_'.join(x.split('tmalign_')[1].split('_')[1:]).replace('.txt', ''))

    upgma_t = Tree(upgma_tree)
    foldtree_t = Tree(foldtree_tree)
    ref_ligands = pd.read_csv(output_csv)

    def get_closest(tree, berghia_id, ref_ids):
        distances = {ref_id: tree.get_distance(berghia_id, ref_id) if ref_id in tree else float('inf') for ref_id in ref_ids}
        return min(distances, key=distances.get)

    proximity_data = []
    for berghia_id in scores['berghia_id'].unique():
        tm_top = scores[scores['berghia_id'] == berghia_id].sort_values('tm_score', ascending=False).iloc[0]['ref_id']
        upgma_ref = get_closest(upgma_t, berghia_id, ref_ligands['ref_id'])
        foldtree_ref = get_closest(foldtree_t, berghia_id, ref_ligands['ref_id'])
        refs = [tm_top, upgma_ref, foldtree_ref]
        closest_ref = max(set(refs), key=refs.count) if len(set(refs)) > 1 else tm_top
        ligands = ref_ligands[ref_ligands['ref_id'] == closest_ref]['ligands'].values[0]
        proximity_data.append({'berghia_id': berghia_id, 'closest_ref': closest_ref, 'putative_ligands': ligands})

    pd.DataFrame(proximity_data).to_csv(putative_ligands_csv, index=False)
