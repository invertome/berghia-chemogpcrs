#!/usr/bin/env python3
# fetch_ligands.py
# Purpose: Fetch GPCR references with known ligands or structures from GPCRdb.
# Inputs: Output dir ($1), output CSV ($2), search terms ($3), species ($4), [family ($5)]
# Outputs: Reference ligands CSV (ref_id, pdb_id, ligands), PDB files in output dir
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
import os
import requests
import pandas as pd
from Bio import PDB

ref_dir = sys.argv[1]
output_csv = sys.argv[2]
search_terms = sys.argv[3].split(',')
species = sys.argv[4].split(',')
# Optional family parameter: can be "Class_A", "Class_B1", "Class_C", "Adhesion", etc.
# or "all" to fetch all GPCR families, or comma-separated list of families
# Default to "all" for broad chemoreceptor searches
gpcr_families = sys.argv[5].split(',') if len(sys.argv) > 5 else os.getenv('GPCRDB_FAMILIES', 'all').split(',')

# Fetch GPCR list from specified families
all_gpcrs = []
if 'all' in gpcr_families:
    # Fetch all receptors (no family filter)
    response = requests.get("https://gpcrdb.org/services/receptorlist/")
    if response.status_code == 200:
        all_gpcrs = response.json()
    else:
        print(f"Warning: Could not fetch all GPCRs: {response.status_code}", file=sys.stderr)
else:
    for family in gpcr_families:
        family = family.strip()
        response = requests.get(f"https://gpcrdb.org/services/receptorlist/?family={family}")
        if response.status_code == 200:
            all_gpcrs.extend(response.json())
        else:
            print(f"Warning: Could not fetch GPCRdb family {family}: {response.status_code}", file=sys.stderr)

if not all_gpcrs:
    print("Error: Could not fetch any GPCR data from GPCRdb")
    sys.exit(1)

# Filter by search terms and species
ref_candidates = [g for g in all_gpcrs if any(term.lower() in g.get('name', '').lower() for term in search_terms) and
                  any(sp in g.get('species', '') for sp in species)]
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
