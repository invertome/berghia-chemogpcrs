#!/usr/bin/env python3
# rename_references_with_taxid.py
# Purpose: Rename reference FASTA files to include taxonomy ID prefix
# Input: Directory containing reference files named like "Species_name.faa"
# Output: Files renamed to "taxid_Species_name.faa"
# Author: Jorge L. Perez-Moreno, Ph.D.

"""
Rename Reference Files with Taxonomy IDs

This script renames reference FASTA files to include their NCBI taxonomy ID as a prefix.
It extracts species names from filenames (e.g., "Biomphalaria_straminea.faa") and looks
up the corresponding taxonomy ID from NCBI.

Features:
- Uses local NCBI taxonomy database (via ete3) when available for fast lookups
- Falls back to NCBI Entrez API when local database is unavailable
- Handles species name variations (underscores, spaces)
- Generates a mapping file for reference
- Supports dry-run mode to preview changes
"""

import os
import sys
import re
import time
import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()


def get_ncbi_credentials() -> Tuple[str, str]:
    """Get NCBI API credentials from environment."""
    api_key = os.environ.get('NCBI_API_KEY', '')
    email = os.environ.get('NCBI_EMAIL', 'user@example.com')
    return api_key, email


def setup_entrez():
    """Configure Biopython Entrez with credentials."""
    from Bio import Entrez
    api_key, email = get_ncbi_credentials()
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    return Entrez


def get_local_taxonomy_db(db_dir: Optional[Path]) -> Optional[object]:
    """
    Get local NCBI taxonomy database if available.

    Args:
        db_dir: Directory containing the local database

    Returns:
        NCBITaxa object or None if not available
    """
    if db_dir is None:
        return None

    db_file = db_dir / "taxa.sqlite"
    if not db_file.exists():
        return None

    try:
        from ete3 import NCBITaxa
        ncbi = NCBITaxa(dbfile=str(db_file))
        return ncbi
    except ImportError:
        print("Warning: ete3 not installed. Using API lookups instead.", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Warning: Could not load local taxonomy database: {e}", file=sys.stderr)
        return None


def species_name_from_filename(filename: str) -> str:
    """
    Extract species name from filename.

    Args:
        filename: Like "Biomphalaria_straminea.faa"

    Returns:
        Species name like "Biomphalaria straminea"
    """
    # Remove extension
    base = Path(filename).stem
    # Replace underscores with spaces
    species = base.replace('_', ' ')
    return species


def lookup_taxid_local(ncbi_taxa, species_name: str) -> Optional[int]:
    """
    Look up taxonomy ID using local database.

    Args:
        ncbi_taxa: ete3 NCBITaxa object
        species_name: Species name like "Biomphalaria straminea"

    Returns:
        Taxonomy ID or None if not found
    """
    try:
        name_dict = ncbi_taxa.get_name_translator([species_name])
        if species_name in name_dict:
            taxids = name_dict[species_name]
            if taxids:
                return taxids[0]

        # Try with different name formats
        # Sometimes names have additional qualifiers
        parts = species_name.split()
        if len(parts) >= 2:
            binomial = f"{parts[0]} {parts[1]}"
            name_dict = ncbi_taxa.get_name_translator([binomial])
            if binomial in name_dict:
                taxids = name_dict[binomial]
                if taxids:
                    return taxids[0]

        return None
    except Exception:
        return None


def lookup_taxid_api(species_name: str, entrez) -> Optional[int]:
    """
    Look up taxonomy ID using NCBI Entrez API.

    Args:
        species_name: Species name like "Biomphalaria straminea"
        entrez: Configured Entrez module

    Returns:
        Taxonomy ID or None if not found
    """
    try:
        # Search taxonomy database
        handle = entrez.esearch(db="taxonomy", term=f'"{species_name}"[Scientific Name]')
        record = entrez.read(handle)
        handle.close()

        if record["IdList"]:
            return int(record["IdList"][0])

        # Try without quotes for partial matches
        parts = species_name.split()
        if len(parts) >= 2:
            binomial = f"{parts[0]} {parts[1]}"
            handle = entrez.esearch(db="taxonomy", term=f'{binomial}[Scientific Name]')
            record = entrez.read(handle)
            handle.close()

            if record["IdList"]:
                return int(record["IdList"][0])

        return None
    except Exception as e:
        print(f"  API error for '{species_name}': {e}", file=sys.stderr)
        return None


def lookup_taxid(species_name: str, ncbi_taxa=None, entrez=None,
                 use_api: bool = True) -> Optional[int]:
    """
    Look up taxonomy ID, trying local database first then API.

    Args:
        species_name: Species name to look up
        ncbi_taxa: Local NCBITaxa object (optional)
        entrez: Entrez module (optional)
        use_api: Whether to fall back to API

    Returns:
        Taxonomy ID or None
    """
    # Try local database first
    if ncbi_taxa is not None:
        taxid = lookup_taxid_local(ncbi_taxa, species_name)
        if taxid:
            return taxid

    # Fall back to API
    if use_api and entrez is not None:
        time.sleep(0.15)  # Rate limiting (10/sec with API key)
        return lookup_taxid_api(species_name, entrez)

    return None


def find_reference_files(ref_dir: Path) -> List[Path]:
    """
    Find all FASTA files in reference directory recursively.

    Args:
        ref_dir: Root directory to search

    Returns:
        List of FASTA file paths
    """
    fasta_extensions = ['.faa', '.fa', '.fasta', '.pep']
    files = []

    for ext in fasta_extensions:
        files.extend(ref_dir.rglob(f'*{ext}'))

    # Filter out already renamed files (those starting with a number)
    files = [f for f in files if not re.match(r'^\d+_', f.name)]

    return sorted(files)


def rename_file(old_path: Path, taxid: int, dry_run: bool = False) -> Tuple[Path, bool]:
    """
    Rename a file to include taxid prefix.

    Args:
        old_path: Original file path
        taxid: Taxonomy ID to add as prefix
        dry_run: If True, don't actually rename

    Returns:
        Tuple of (new_path, success)
    """
    new_name = f"{taxid}_{old_path.name}"
    new_path = old_path.parent / new_name

    if new_path.exists():
        print(f"  Warning: {new_path.name} already exists, skipping", file=sys.stderr)
        return old_path, False

    if not dry_run:
        old_path.rename(new_path)

    return new_path, True


def process_references(ref_dir: Path, db_dir: Optional[Path] = None,
                       dry_run: bool = False, use_api: bool = True,
                       verbose: bool = True) -> Dict[str, dict]:
    """
    Process all reference files and rename with taxid prefix.

    Args:
        ref_dir: Directory containing reference files
        db_dir: Directory with local taxonomy database
        dry_run: Preview changes without renaming
        use_api: Fall back to NCBI API when local lookup fails
        verbose: Print progress

    Returns:
        Dictionary with results for each file
    """
    results = {}

    # Setup databases
    ncbi_taxa = get_local_taxonomy_db(db_dir)
    entrez = None

    if use_api:
        try:
            entrez = setup_entrez()
            if verbose:
                api_key, _ = get_ncbi_credentials()
                print(f"NCBI API configured (API key: {'yes' if api_key else 'no'})")
        except ImportError:
            print("Warning: Biopython not installed. API lookups disabled.", file=sys.stderr)
            entrez = None

    if ncbi_taxa:
        print(f"Using local taxonomy database: {db_dir / 'taxa.sqlite'}")
    elif not entrez:
        print("ERROR: No taxonomy lookup method available.", file=sys.stderr)
        print("Either setup local database (python setup_databases.py) or install Biopython", file=sys.stderr)
        return results

    # Find all reference files
    files = find_reference_files(ref_dir)

    if verbose:
        print(f"\nFound {len(files)} reference files to process")
        if dry_run:
            print("DRY RUN - no files will be renamed\n")
        print("-" * 60)

    success_count = 0
    failed_count = 0

    for fpath in files:
        species = species_name_from_filename(fpath.name)

        if verbose:
            print(f"Processing: {fpath.name}")
            print(f"  Species: {species}")

        # Look up taxid
        taxid = lookup_taxid(species, ncbi_taxa, entrez, use_api)

        if taxid:
            if verbose:
                print(f"  TaxID: {taxid}")

            new_path, renamed = rename_file(fpath, taxid, dry_run)

            results[str(fpath)] = {
                'species': species,
                'taxid': taxid,
                'old_path': str(fpath),
                'new_path': str(new_path),
                'renamed': renamed,
                'status': 'success'
            }

            if renamed:
                success_count += 1
                if verbose:
                    action = "Would rename" if dry_run else "Renamed"
                    print(f"  {action} to: {new_path.name}")
        else:
            if verbose:
                print(f"  ERROR: Could not find taxid for '{species}'")

            results[str(fpath)] = {
                'species': species,
                'taxid': None,
                'old_path': str(fpath),
                'new_path': None,
                'renamed': False,
                'status': 'taxid_not_found'
            }
            failed_count += 1

        if verbose:
            print()

    if verbose:
        print("-" * 60)
        print(f"Summary: {success_count} renamed, {failed_count} failed")

    return results


def save_mapping(results: Dict[str, dict], output_file: Path):
    """Save the rename mapping to CSV for reference."""
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'species', 'taxid', 'old_filename', 'new_filename', 'status'
        ])
        writer.writeheader()

        for old_path, info in results.items():
            old_name = Path(old_path).name
            new_name = Path(info['new_path']).name if info['new_path'] else ''

            writer.writerow({
                'species': info['species'],
                'taxid': info['taxid'] or '',
                'old_filename': old_name,
                'new_filename': new_name,
                'status': info['status']
            })


def main():
    parser = argparse.ArgumentParser(
        description="Rename reference FASTA files to include taxonomy ID prefix",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Preview renames (dry run)
  python rename_references_with_taxid.py references/nath_et_al --dry-run

  # Rename all files
  python rename_references_with_taxid.py references/nath_et_al

  # Use local database only (no API calls)
  python rename_references_with_taxid.py references/nath_et_al --local-only --db-dir databases

  # Save mapping to file
  python rename_references_with_taxid.py references/nath_et_al --output-map ref_mapping.csv
        """
    )

    parser.add_argument('ref_dir', type=Path,
                       help='Directory containing reference files')
    parser.add_argument('--db-dir', type=Path, default=None,
                       help='Directory with local taxonomy database (default: from LOCAL_DB_DIR env)')
    parser.add_argument('--dry-run', '-n', action='store_true',
                       help='Preview changes without renaming')
    parser.add_argument('--local-only', action='store_true',
                       help='Only use local database, no API calls')
    parser.add_argument('--output-map', '-o', type=Path, default=None,
                       help='Save rename mapping to this CSV file')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress progress output')

    args = parser.parse_args()

    # Check reference directory
    if not args.ref_dir.exists():
        print(f"ERROR: Reference directory not found: {args.ref_dir}", file=sys.stderr)
        sys.exit(1)

    # Get database directory from env if not specified
    db_dir = args.db_dir
    if db_dir is None:
        db_dir_env = os.environ.get('LOCAL_DB_DIR')
        if db_dir_env:
            db_dir = Path(db_dir_env)

    # Process files
    results = process_references(
        args.ref_dir,
        db_dir=db_dir,
        dry_run=args.dry_run,
        use_api=not args.local_only,
        verbose=not args.quiet
    )

    # Save mapping if requested
    if args.output_map:
        save_mapping(results, args.output_map)
        if not args.quiet:
            print(f"\nMapping saved to: {args.output_map}")

    # Exit with error if any lookups failed
    failed = sum(1 for r in results.values() if r['status'] != 'success')
    if failed > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
