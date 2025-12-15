#!/usr/bin/env python3
# fetch_ligands.py
# Purpose: Fetch GPCR references with known ligands or structures from GPCRdb.
# Inputs: Output dir ($1), output CSV ($2), search terms ($3), species ($4), [family ($5)]
# Outputs: Reference ligands CSV (ref_id, pdb_id, ligands), PDB files in output dir
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
import os
import json
import time
from pathlib import Path
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


def load_local_gpcrdb_data(local_db_dir: Path, search_terms: list, species: list, gpcr_families: list):
    """
    Load GPCR data from local database (created by setup_databases.py).

    Returns list of filtered GPCR entries or None if local data unavailable.
    """
    gpcrdb_dir = local_db_dir / "gpcrdb"
    receptors_file = gpcrdb_dir / "receptors.json"

    if not receptors_file.exists():
        return None

    try:
        all_gpcrs = json.loads(receptors_file.read_text())
        print(f"Loaded {len(all_gpcrs)} receptors from local database", file=sys.stderr)

        # Filter by family if not 'all'
        if 'all' not in gpcr_families:
            all_gpcrs = [g for g in all_gpcrs
                        if any(fam.lower() in g.get('family', '').lower() for fam in gpcr_families)]

        return all_gpcrs

    except Exception as e:
        print(f"Warning: Could not load local GPCRdb data: {e}", file=sys.stderr)
        return None


def load_local_structures(local_db_dir: Path):
    """Load structure metadata from local database."""
    structures_file = local_db_dir / "gpcrdb" / "structures.json"

    if structures_file.exists():
        try:
            return json.loads(structures_file.read_text())
        except Exception:
            pass
    return None


def get_local_pdb_file(local_db_dir: Path, pdb_id: str):
    """Check if PDB file exists in local database."""
    pdb_file = local_db_dir / "gpcrdb" / "structures_pdb" / f"{pdb_id.lower()}.pdb"
    return pdb_file if pdb_file.exists() else None


session, timeout = create_session()

ref_dir = sys.argv[1]
output_csv = sys.argv[2]
search_terms = sys.argv[3].split(',')
species = sys.argv[4].split(',')
# Optional family parameter: can be "Class_A", "Class_B1", "Class_C", "Adhesion", etc.
# or "all" to fetch all GPCR families, or comma-separated list of families
# Default to "all" for broad chemoreceptor searches
gpcr_families = sys.argv[5].split(',') if len(sys.argv) > 5 else os.getenv('GPCRDB_FAMILIES', 'all').split(',')

# Check for local database first (API-independent mode)
LOCAL_DB_DIR = os.getenv('LOCAL_DB_DIR', '')
local_db_path = Path(LOCAL_DB_DIR) if LOCAL_DB_DIR else None
use_local_db = local_db_path and (local_db_path / "gpcrdb" / "receptors.json").exists()

# Load GPCR data - prefer local, fallback to API
all_gpcrs = None
local_structures = None

if use_local_db:
    print("Using local GPCRdb database (API-independent mode)", file=sys.stderr)
    all_gpcrs = load_local_gpcrdb_data(local_db_path, search_terms, species, gpcr_families)
    local_structures = load_local_structures(local_db_path)

if all_gpcrs is None:
    print("Fetching from GPCRdb API...", file=sys.stderr)
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
    print("Error: Could not fetch any GPCR data (neither local nor API)", file=sys.stderr)
    sys.exit(1)

# Filter by search terms and species
ref_candidates = [g for g in all_gpcrs if any(term.lower() in g.get('name', '').lower() for term in search_terms) and
                  any(sp in g.get('species', '') for sp in species)]
ref_ids = [g['uniprot'] for g in ref_candidates]

# Download PDBs and extract metadata
pdb_list = PDB.PDBList()
parser = PDB.PDBParser(QUIET=True)
ref_data = []

# Build structure lookup from local data if available
structure_by_uniprot = {}
if local_structures:
    for struct in local_structures:
        uniprot = struct.get('protein', '')
        if uniprot:
            if uniprot not in structure_by_uniprot:
                structure_by_uniprot[uniprot] = []
            structure_by_uniprot[uniprot].append(struct)

for ref_id in ref_ids:
    pdb_id = None
    structures = None

    # Try local structure data first
    if ref_id in structure_by_uniprot:
        structures = structure_by_uniprot[ref_id]
        if structures:
            pdb_id = structures[0].get('pdb_code', '')
    else:
        # Fallback to API
        pdb_response = fetch_with_retry(session, f"https://gpcrdb.org/services/structure/{ref_id}", timeout, f"structure for {ref_id}")
        if pdb_response:
            structures = pdb_response.json()
            if structures:
                pdb_id = structures[0].get('pdb_code', '')

    if pdb_id:
        try:
            # Check local PDB file first
            local_pdb = get_local_pdb_file(local_db_path, pdb_id) if use_local_db else None

            if local_pdb:
                pdb_file = str(local_pdb)
                # Copy to output directory
                import shutil
                dest_file = Path(ref_dir) / f"{pdb_id.lower()}.pdb"
                if not dest_file.exists():
                    shutil.copy(local_pdb, dest_file)
            else:
                # Download from PDB
                pdb_file = pdb_list.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=ref_dir)

            structure = parser.get_structure(pdb_id, pdb_file)
            ligands = [res.resname for res in structure[0].get_residues() if res.id[0].startswith('H_')]
            if ligands:
                ref_data.append({'ref_id': ref_id, 'pdb_id': pdb_id, 'ligands': ','.join(ligands)})
        except Exception as e:
            print(f"Warning: Could not process PDB {pdb_id}: {e}", file=sys.stderr)

    time.sleep(0.1)  # Reduced rate limiting when using local data

# Save metadata
pd.DataFrame(ref_data).to_csv(output_csv, index=False)
print(f"Processed {len(ref_data)} structures with ligand data", file=sys.stderr)
