#!/usr/bin/env python3
# setup_databases.py
# Purpose: Download and setup local databases for API-independent pipeline execution.
# Outputs: Local NCBI Taxonomy database, pre-downloaded GPCRdb reference data
# Author: Jorge L. Perez-Moreno, Ph.D.

"""
Database Setup Script

This script pre-downloads all required databases to enable offline/API-independent
pipeline execution:

1. NCBI Taxonomy database (via ete3)
2. GPCRdb reference sequences and structures metadata
3. Reference structure files from PDB/AlphaFold DB

Run this script ONCE before running the pipeline to ensure all data is available locally.
"""

import os
import sys
import json
import time
import argparse
from pathlib import Path
from typing import Dict, List, Optional
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


def create_robust_session(retries: int = 5, backoff_factor: float = 1.0, timeout: int = 60):
    """Create a requests session with retry logic."""
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


def setup_ncbi_taxonomy(db_dir: Path, verbose: bool = True) -> bool:
    """
    Download and setup NCBI Taxonomy database for ete3.

    This enables offline taxonomy lookups without API calls.
    """
    if verbose:
        print("=" * 60)
        print("Setting up NCBI Taxonomy database...")
        print("=" * 60)

    try:
        from ete3 import NCBITaxa

        # Set custom database path
        db_file = db_dir / "taxa.sqlite"

        if verbose:
            print(f"Database location: {db_file}")

        # Initialize NCBITaxa - this will download if not present
        # or update if database is old
        ncbi = NCBITaxa(dbfile=str(db_file))

        # Check if database needs updating (older than 30 days)
        if db_file.exists():
            import datetime
            mod_time = datetime.datetime.fromtimestamp(db_file.stat().st_mtime)
            age_days = (datetime.datetime.now() - mod_time).days

            if age_days > 30:
                if verbose:
                    print(f"Database is {age_days} days old. Updating...")
                ncbi.update_taxonomy_database()
            else:
                if verbose:
                    print(f"Database is {age_days} days old. No update needed.")
        else:
            if verbose:
                print("Downloading NCBI Taxonomy database (this may take a few minutes)...")
            ncbi.update_taxonomy_database()

        # Verify database works
        test_lineage = ncbi.get_lineage(9606)  # Human
        if test_lineage:
            if verbose:
                print(f"Database verification successful (tested with taxid 9606)")
                print(f"Database contains taxonomy data up to NCBI's latest dump")
            return True
        else:
            print("ERROR: Database verification failed", file=sys.stderr)
            return False

    except ImportError:
        print("ERROR: ete3 not installed. Run: pip install ete3", file=sys.stderr)
        return False
    except Exception as e:
        print(f"ERROR: Failed to setup taxonomy database: {e}", file=sys.stderr)
        return False


def fetch_gpcrdb_families(session, timeout: int, verbose: bool = True) -> List[Dict]:
    """Fetch all GPCR family information from GPCRdb."""
    families = []

    try:
        # Get receptor families
        url = "https://gpcrdb.org/services/receptorfamily/"
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        families = response.json()

        if verbose:
            print(f"Retrieved {len(families)} GPCR families from GPCRdb")

    except Exception as e:
        if verbose:
            print(f"Warning: Could not fetch GPCRdb families: {e}")

    return families


def fetch_gpcrdb_receptors(session, timeout: int, species_filter: Optional[List[str]] = None,
                          verbose: bool = True) -> List[Dict]:
    """
    Fetch receptor information from GPCRdb.

    Args:
        session: requests session
        timeout: request timeout
        species_filter: Optional list of species to filter (e.g., ['Homo sapiens', 'Aplysia californica'])
        verbose: Print progress

    Returns:
        List of receptor dictionaries
    """
    receptors = []

    try:
        # Get all proteins
        url = "https://gpcrdb.org/services/proteinfamily/proteins/"
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        all_receptors = response.json()

        if species_filter:
            receptors = [r for r in all_receptors
                        if any(sp.lower() in r.get('species', '').lower()
                              for sp in species_filter)]
        else:
            receptors = all_receptors

        if verbose:
            print(f"Retrieved {len(receptors)} receptors from GPCRdb")
            if species_filter:
                print(f"  (filtered by species: {species_filter})")

    except Exception as e:
        if verbose:
            print(f"Warning: Could not fetch GPCRdb receptors: {e}")

    return receptors


def fetch_gpcrdb_structures(session, timeout: int, verbose: bool = True) -> List[Dict]:
    """Fetch all available structure information from GPCRdb."""
    structures = []

    try:
        url = "https://gpcrdb.org/services/structure/"
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        structures = response.json()

        if verbose:
            print(f"Retrieved {len(structures)} structures from GPCRdb")

    except Exception as e:
        if verbose:
            print(f"Warning: Could not fetch GPCRdb structures: {e}")

    return structures


def fetch_gpcrdb_ligands(session, timeout: int, verbose: bool = True) -> List[Dict]:
    """Fetch ligand information from GPCRdb."""
    ligands = []

    try:
        url = "https://gpcrdb.org/services/ligands/"
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        ligands = response.json()

        if verbose:
            print(f"Retrieved {len(ligands)} ligands from GPCRdb")

    except Exception as e:
        if verbose:
            print(f"Warning: Could not fetch GPCRdb ligands: {e}")

    return ligands


def download_pdb_structure(session, pdb_id: str, output_dir: Path, timeout: int) -> bool:
    """Download a PDB structure file."""
    output_file = output_dir / f"{pdb_id.lower()}.pdb"

    if output_file.exists():
        return True  # Already downloaded

    try:
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        response = session.get(url, timeout=timeout)
        response.raise_for_status()

        output_file.write_text(response.text)
        return True

    except Exception:
        return False


def setup_gpcrdb_data(db_dir: Path, download_structures: bool = False,
                     species_filter: Optional[List[str]] = None,
                     max_structures: int = 100,
                     verbose: bool = True) -> bool:
    """
    Download GPCRdb data for offline use.

    Args:
        db_dir: Directory to store data
        download_structures: Whether to download actual PDB files
        species_filter: Optional species to focus on
        max_structures: Maximum structures to download
        verbose: Print progress
    """
    if verbose:
        print("=" * 60)
        print("Setting up GPCRdb local database...")
        print("=" * 60)

    gpcrdb_dir = db_dir / "gpcrdb"
    gpcrdb_dir.mkdir(parents=True, exist_ok=True)

    session, timeout = create_robust_session()

    # Fetch and save all data types
    data_files = {}

    # 1. Receptor families
    if verbose:
        print("\n[1/4] Fetching receptor families...")
    families = fetch_gpcrdb_families(session, timeout, verbose)
    if families:
        data_files['families'] = gpcrdb_dir / "families.json"
        data_files['families'].write_text(json.dumps(families, indent=2))

    # 2. Receptors
    if verbose:
        print("\n[2/4] Fetching receptors...")
    receptors = fetch_gpcrdb_receptors(session, timeout, species_filter, verbose)
    if receptors:
        data_files['receptors'] = gpcrdb_dir / "receptors.json"
        data_files['receptors'].write_text(json.dumps(receptors, indent=2))

    # 3. Structures metadata
    if verbose:
        print("\n[3/4] Fetching structure metadata...")
    structures = fetch_gpcrdb_structures(session, timeout, verbose)
    if structures:
        data_files['structures'] = gpcrdb_dir / "structures.json"
        data_files['structures'].write_text(json.dumps(structures, indent=2))

    # 4. Ligands
    if verbose:
        print("\n[4/4] Fetching ligand information...")
    ligands = fetch_gpcrdb_ligands(session, timeout, verbose)
    if ligands:
        data_files['ligands'] = gpcrdb_dir / "ligands.json"
        data_files['ligands'].write_text(json.dumps(ligands, indent=2))

    # 5. Optionally download actual structure files
    if download_structures and structures:
        if verbose:
            print(f"\nDownloading up to {max_structures} PDB structure files...")

        pdb_dir = gpcrdb_dir / "structures_pdb"
        pdb_dir.mkdir(exist_ok=True)

        downloaded = 0
        failed = 0

        for struct in structures[:max_structures]:
            pdb_id = struct.get('pdb_code', '')
            if pdb_id:
                if download_pdb_structure(session, pdb_id, pdb_dir, timeout):
                    downloaded += 1
                else:
                    failed += 1

                if verbose and (downloaded + failed) % 10 == 0:
                    print(f"  Progress: {downloaded} downloaded, {failed} failed")

                # Rate limiting
                time.sleep(0.5)

        if verbose:
            print(f"  Final: {downloaded} downloaded, {failed} failed")

    # Create metadata file
    metadata = {
        'download_date': time.strftime('%Y-%m-%d %H:%M:%S'),
        'species_filter': species_filter,
        'n_families': len(families),
        'n_receptors': len(receptors),
        'n_structures': len(structures),
        'n_ligands': len(ligands),
        'files': {k: str(v) for k, v in data_files.items()}
    }

    metadata_file = gpcrdb_dir / "metadata.json"
    metadata_file.write_text(json.dumps(metadata, indent=2))

    if verbose:
        print(f"\nGPCRdb data saved to: {gpcrdb_dir}")
        print(f"Metadata file: {metadata_file}")

    return True


def create_chemoreceptor_reference_list(db_dir: Path, verbose: bool = True) -> bool:
    """
    Create a curated list of known chemoreceptor references from GPCRdb data.

    This file is used by the ranking algorithm to weight references appropriately.
    """
    gpcrdb_dir = db_dir / "gpcrdb"
    receptors_file = gpcrdb_dir / "receptors.json"

    if not receptors_file.exists():
        if verbose:
            print("Warning: GPCRdb receptors file not found. Run setup_gpcrdb_data first.")
        return False

    receptors = json.loads(receptors_file.read_text())

    # Keywords that indicate chemoreceptor function
    chemoreceptor_keywords = [
        'olfactory', 'odorant', 'taste', 'gustatory', 'vomeronasal',
        'pheromone', 'chemosensory', 'odorant receptor', 'trace amine',
        'formyl peptide', 'chemokine', 'chemoattractant'
    ]

    chemoreceptors = []
    other_gpcrs = []

    for receptor in receptors:
        name = receptor.get('name', '').lower()
        family = receptor.get('family', '').lower()

        is_chemoreceptor = any(kw in name or kw in family for kw in chemoreceptor_keywords)

        entry = {
            'entry_name': receptor.get('entry_name', ''),
            'name': receptor.get('name', ''),
            'family': receptor.get('family', ''),
            'species': receptor.get('species', ''),
            'accession': receptor.get('accession', '')
        }

        if is_chemoreceptor:
            chemoreceptors.append(entry)
        else:
            other_gpcrs.append(entry)

    # Save categorized references
    ref_file = gpcrdb_dir / "reference_categories.json"
    ref_data = {
        'chemoreceptors': chemoreceptors,
        'other_gpcrs': other_gpcrs,
        'chemoreceptor_keywords': chemoreceptor_keywords
    }
    ref_file.write_text(json.dumps(ref_data, indent=2))

    if verbose:
        print(f"\nReference categorization:")
        print(f"  Chemoreceptors: {len(chemoreceptors)}")
        print(f"  Other GPCRs: {len(other_gpcrs)}")
        print(f"  Saved to: {ref_file}")

    return True


def verify_setup(db_dir: Path, verbose: bool = True) -> Dict[str, bool]:
    """Verify that all required databases are properly set up."""
    results = {}

    if verbose:
        print("=" * 60)
        print("Verifying database setup...")
        print("=" * 60)

    # Check NCBI Taxonomy
    taxa_db = db_dir / "taxa.sqlite"
    results['ncbi_taxonomy'] = taxa_db.exists()
    if verbose:
        status = "OK" if results['ncbi_taxonomy'] else "MISSING"
        print(f"  NCBI Taxonomy database: {status}")

    # Check GPCRdb data
    gpcrdb_dir = db_dir / "gpcrdb"
    gpcrdb_files = ['families.json', 'receptors.json', 'structures.json', 'ligands.json']

    results['gpcrdb_data'] = all((gpcrdb_dir / f).exists() for f in gpcrdb_files)
    if verbose:
        status = "OK" if results['gpcrdb_data'] else "INCOMPLETE"
        print(f"  GPCRdb data: {status}")

    # Check reference categories
    ref_cat = gpcrdb_dir / "reference_categories.json"
    results['reference_categories'] = ref_cat.exists()
    if verbose:
        status = "OK" if results['reference_categories'] else "MISSING"
        print(f"  Reference categories: {status}")

    # Overall status
    results['all_ok'] = all(results.values())
    if verbose:
        print("-" * 60)
        overall = "READY" if results['all_ok'] else "INCOMPLETE"
        print(f"  Overall status: {overall}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Setup local databases for API-independent pipeline execution",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full setup with default options
  python setup_databases.py

  # Setup with structure downloads for specific species
  python setup_databases.py --download-structures --species "Aplysia" "Lottia"

  # Verify existing setup
  python setup_databases.py --verify-only

  # Specify custom database directory
  python setup_databases.py --db-dir /path/to/databases
        """
    )

    parser.add_argument('--db-dir', type=Path, default=Path('databases'),
                       help='Directory to store databases (default: ./databases)')
    parser.add_argument('--download-structures', action='store_true',
                       help='Download PDB structure files (slower, requires more space)')
    parser.add_argument('--max-structures', type=int, default=100,
                       help='Maximum structures to download (default: 100)')
    parser.add_argument('--species', nargs='+', default=None,
                       help='Species to filter receptors (e.g., "Aplysia" "Lottia")')
    parser.add_argument('--verify-only', action='store_true',
                       help='Only verify existing setup, do not download')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress output')

    args = parser.parse_args()
    verbose = not args.quiet

    # Create database directory
    args.db_dir.mkdir(parents=True, exist_ok=True)

    if args.verify_only:
        results = verify_setup(args.db_dir, verbose)
        sys.exit(0 if results['all_ok'] else 1)

    # Run setup
    success = True

    # 1. NCBI Taxonomy
    if not setup_ncbi_taxonomy(args.db_dir, verbose):
        success = False

    # 2. GPCRdb data
    if not setup_gpcrdb_data(args.db_dir,
                            download_structures=args.download_structures,
                            species_filter=args.species,
                            max_structures=args.max_structures,
                            verbose=verbose):
        success = False

    # 3. Reference categorization
    if not create_chemoreceptor_reference_list(args.db_dir, verbose):
        success = False

    # Final verification
    print()
    results = verify_setup(args.db_dir, verbose)

    if results['all_ok']:
        print("\nSetup complete! You can now run the pipeline offline.")
        print(f"\nAdd this to your config.sh:")
        print(f'  export LOCAL_DB_DIR="{args.db_dir.absolute()}"')
    else:
        print("\nSetup incomplete. Check errors above.", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
