#!/usr/bin/env python3
# fetch_ligands.py
# Purpose: Fetch GPCR references with known ligands or structures from GPCRdb.
# Inputs: Output dir ($1), output CSV ($2), search terms ($3), species ($4), [family ($5)]
# Outputs: Reference ligands CSV (ref_id, pdb_id, ligands), PDB files in output dir
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
import os
import time
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import pandas as pd
from Bio import PDB

# Configure robust HTTP session with retry logic
def create_session(retries=3, backoff_factor=0.5, timeout=30):
    """Create requests session with retry logic and timeout."""
    session = requests.Session()
    retry_strategy = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session, timeout

def fetch_with_retry(session, url, timeout, description=""):
    """Fetch URL with retry logic and error handling."""
    try:
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        return response
    except requests.exceptions.Timeout:
        print(f"Warning: Timeout fetching {description or url}", file=sys.stderr)
        return None
    except requests.exceptions.RequestException as e:
        print(f"Warning: Error fetching {description or url}: {e}", file=sys.stderr)
        return None

session, timeout = create_session()

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
    response = fetch_with_retry(session, "https://gpcrdb.org/services/receptorlist/", timeout, "all GPCRs")
    if response:
        all_gpcrs = response.json()
else:
    for family in gpcr_families:
        family = family.strip()
        response = fetch_with_retry(session, f"https://gpcrdb.org/services/receptorlist/?family={family}", timeout, f"family {family}")
        if response:
            all_gpcrs.extend(response.json())
        time.sleep(0.2)  # Rate limiting between requests

if not all_gpcrs:
    print("Error: Could not fetch any GPCR data from GPCRdb")
    sys.exit(1)

# Filter by search terms and species
ref_candidates = [g for g in all_gpcrs if any(term.lower() in g.get('name', '').lower() for term in search_terms) and
                  any(sp in g.get('species', '') for sp in species)]
ref_ids = [g['uniprot'] for g in ref_candidates]

# Download PDBs and extract metadata
pdb_list = PDB.PDBList()
parser = PDB.PDBParser(QUIET=True)
ref_data = []

for ref_id in ref_ids:
    pdb_response = fetch_with_retry(session, f"https://gpcrdb.org/services/structure/{ref_id}", timeout, f"structure for {ref_id}")
    if pdb_response:
        structures = pdb_response.json()
        if structures:
            pdb_id = structures[0]['pdb_code']
            try:
                pdb_file = pdb_list.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=ref_dir)
                structure = parser.get_structure(pdb_id, pdb_file)
                ligands = [res.resname for res in structure[0].get_residues() if res.id[0].startswith('H_')]
                if ligands:
                    ref_data.append({'ref_id': ref_id, 'pdb_id': pdb_id, 'ligands': ','.join(ligands)})
            except Exception as e:
                print(f"Warning: Could not process PDB {pdb_id}: {e}", file=sys.stderr)
    time.sleep(0.2)  # Rate limiting

# Save metadata
pd.DataFrame(ref_data).to_csv(output_csv, index=False)
