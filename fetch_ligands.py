#!/usr/bin/env python3
# fetch_ligands.py
# Purpose: Fetch GPCR references with known ligands or structures from GPCRdb.
# Inputs: Output dir ($1), output CSV ($2), search terms ($3), species ($4)
# Outputs: Reference ligands CSV (ref_id, pdb_id, ligands), PDB files in output dir
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
import requests
import pandas as pd
from Bio import PDB

ref_dir = sys.argv[1]
output_csv = sys.argv[2]
search_terms = sys.argv[3].split(',')
species = sys.argv[4].split(',')

# Fetch GPCR list
response = requests.get("https://gpcrdb.org/services/receptorlist/?family=Class_A")
if response.status_code != 200:
    print(f"Error fetching GPCRdb data: {response.status_code}")
    sys.exit(1)

gpcrs = response.json()
ref_candidates = [g for g in gpcrs if any(term.lower() in g['name'].lower() for term in search_terms) and 
                  any(sp in g['species'] for sp in species)]
ref_ids = [g['uniprot'] for g in ref_candidates]

# Download PDBs and extract metadata
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

# Save metadata
pd.DataFrame(ref_data).to_csv(output_csv, index=False)
