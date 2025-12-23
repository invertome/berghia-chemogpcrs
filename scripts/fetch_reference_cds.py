#!/usr/bin/env python3
# fetch_reference_cds.py
# Purpose: Fetch CDS (coding sequences) for reference protein sequences
# Input: Reference FASTA files with accession IDs in headers
# Output: CDS FASTA files corresponding to each protein file
# Author: Jorge L. Perez-Moreno, Ph.D.

"""
Fetch CDS Sequences for Reference Proteins

This script retrieves the corresponding coding sequences (CDS) for protein
reference files by:

1. Parsing protein FASTA headers to extract accession IDs
2. Determining the source database (NCBI, ENA, UniProt, DDBJ)
3. Fetching the corresponding nucleotide CDS for each protein
4. Generating CDS FASTA files with matching headers

Accession patterns handled:
- XP_/NP_ (NCBI RefSeq protein) -> linked via XM_/NM_
- TSA contigs (GABXXX, JABXXX, etc.) -> NCBI/ENA nuccore
- GenBank protein (ABC12345.1) -> linked via GenBank CDS
- UniProt (P/Q/O) -> cross-reference to EMBL/GenBank

Note: Some proteins may not have available CDS (e.g., predicted from genomic
scaffolds without explicit CDS annotation). These are logged and skipped.
"""

import os
import sys
import re
import time
import argparse
import json
import csv
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict
from dotenv import load_dotenv

# Load environment variables
load_dotenv()


@dataclass
class AccessionInfo:
    """Information about an accession ID."""
    accession: str
    database: str  # 'ncbi_refseq', 'ncbi_genbank', 'ena', 'uniprot', 'ddbj', 'unknown'
    protein_id: str
    original_header: str
    cds_accession: Optional[str] = None
    cds_sequence: Optional[str] = None
    status: str = 'pending'
    error_message: str = ''


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


def classify_accession(accession: str) -> str:
    """
    Classify an accession to determine its source database.

    Args:
        accession: The accession ID

    Returns:
        Database identifier
    """
    # NCBI RefSeq proteins (AP_, NP_, XP_, YP_, WP_)
    # AP_ = Annotated on alternate assembly
    # NP_ = Protein (curated)
    # XP_ = Protein (predicted/model)
    # YP_ = Protein (bacterial/viral)
    # WP_ = Protein (non-redundant prokaryotic)
    if re.match(r'^[ANXYW]P_\d+', accession):
        return 'ncbi_refseq_protein'

    # NCBI RefSeq transcripts/RNA (NM_, XM_, NR_, XR_)
    # NM_ = mRNA (curated)
    # XM_ = mRNA (predicted/model)
    # NR_ = non-coding RNA (curated)
    # XR_ = non-coding RNA (predicted)
    if re.match(r'^[NX][MR]_\d+', accession):
        return 'ncbi_refseq_transcript'

    # NCBI RefSeq genomic (NC_, NG_, NT_, NW_, NZ_, AC_)
    # NC_ = Complete genomic molecule (chromosome, complete genome, plasmid)
    # NG_ = Incomplete genomic region
    # NT_ = Contig (genomic)
    # NW_ = WGS contig (genomic)
    # NZ_ = WGS supercontig/scaffold
    # AC_ = Alternate assembly (genomic)
    # NS_ = Environmental sequence
    if re.match(r'^(N[CGTWZ]|AC|NS)_\d+', accession):
        return 'ncbi_refseq_genomic'

    # TSA (Transcriptome Shotgun Assembly) - typically 4-6 letter prefix + numbers
    # Examples: GABXXX, GECXXX, JABXXX, etc.
    if re.match(r'^[A-Z]{3,6}\d{8,}', accession):
        return 'ncbi_tsa'

    # GenBank protein (3 letters + 5 digits or 3 letters + 7 digits)
    if re.match(r'^[A-Z]{3}\d{5}(\.\d+)?$', accession):
        return 'ncbi_genbank_protein'

    # GenBank nucleotide
    if re.match(r'^[A-Z]{1,2}\d{5,6}(\.\d+)?$', accession):
        return 'ncbi_genbank_nuc'

    # UniProt (standard and extended formats)
    # Standard: [OPQ][0-9][A-Z0-9]{3}[0-9] (e.g., P12345, Q9XYZ1)
    # Extended: [A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9] (e.g., A0A123BC45)
    if re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', accession):
        return 'uniprot'

    # DDBJ/ENA nucleotide (single letter + 5 numbers)
    if re.match(r'^[A-Z]\d{5}(\.\d+)?$', accession):
        return 'ena_ddbj'

    # WGS/TSA contigs with version
    if re.match(r'^[A-Z]{4,}\d+\.\d+$', accession):
        return 'ncbi_wgs'

    return 'unknown'


def extract_accession_from_header(header: str) -> Optional[str]:
    """
    Extract accession ID from FASTA header.

    Expected format: >prefix_ACCESSION_suffix or >prefix_ACCESSION.version_suffix

    Args:
        header: FASTA header line (with or without >)

    Returns:
        Accession ID or None
    """
    # Remove leading >
    header = header.lstrip('>')

    # Split by underscore
    parts = header.split('_')

    if len(parts) < 2:
        return None

    # The accession is typically the second part
    # Handle cases like: aplcal_XP_005089002.1 (3 parts with XP_)
    # or: alvmar_JABMCL010001063.1_3302

    # All known NCBI RefSeq prefixes that contain underscores
    # Protein prefixes:
    #   AP_ - Annotated on alternate assembly
    #   NP_ - Protein (curated RefSeq)
    #   XP_ - Protein (predicted/model RefSeq)
    #   YP_ - Protein (bacterial/viral)
    #   WP_ - Protein (non-redundant prokaryotic)
    # Transcript/RNA prefixes:
    #   NM_ - mRNA (curated RefSeq)
    #   XM_ - mRNA (predicted/model RefSeq)
    #   NR_ - non-coding RNA (curated)
    #   XR_ - non-coding RNA (predicted)
    # Genomic prefixes:
    #   NC_ - Complete genomic molecule (chromosome, complete genome, plasmid)
    #   NG_ - Incomplete genomic region
    #   NT_ - Contig (genomic)
    #   NW_ - WGS contig (genomic)
    #   NZ_ - WGS supercontig/scaffold
    #   AC_ - Alternate assembly (genomic)
    # Other:
    #   NS_ - Environmental sequence
    REFSEQ_PREFIXES = {
        # Proteins
        'AP', 'NP', 'XP', 'YP', 'WP',
        # Transcripts/RNA
        'NM', 'XM', 'NR', 'XR',
        # Genomic
        'NC', 'NG', 'NT', 'NW', 'NZ', 'AC',
        # Other
        'NS'
    }

    # Check if second part looks like a RefSeq prefix
    if len(parts) >= 3 and parts[1] in REFSEQ_PREFIXES:
        # Reconstruct RefSeq accession: XP_005089002.1
        accession = f"{parts[1]}_{parts[2]}"
        # Remove any trailing suffix after the version
        accession = re.sub(r'(\.\d+).*$', r'\1', accession)
        return accession

    # Otherwise, the accession is the second part
    accession = parts[1]

    # Clean up - remove any trailing non-accession characters
    # But preserve version numbers like .1
    accession = re.match(r'^([A-Za-z0-9_]+\.?\d*)', accession)
    if accession:
        return accession.group(1)

    return None


def parse_fasta_headers(fasta_file: Path) -> List[AccessionInfo]:
    """
    Parse FASTA file and extract accession information.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        List of AccessionInfo objects
    """
    accessions = []

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()
                accession = extract_accession_from_header(header)

                if accession:
                    db = classify_accession(accession)
                    accessions.append(AccessionInfo(
                        accession=accession,
                        database=db,
                        protein_id=header.lstrip('>').split()[0],
                        original_header=header
                    ))

    return accessions


def fetch_cds_ncbi_refseq(accession: str, entrez) -> Optional[Tuple[str, str]]:
    """
    Fetch CDS for NCBI RefSeq protein accession.

    Handles all RefSeq protein prefixes:
    - XP_ (predicted) -> typically linked to XM_ transcript
    - NP_ (curated) -> typically linked to NM_ transcript
    - AP_ (alternate assembly) -> linked nucleotide varies
    - YP_ (bacterial/viral) -> linked to NC_ or other genomic
    - WP_ (non-redundant prokaryotic) -> may link to multiple sources

    Args:
        accession: RefSeq protein accession
        entrez: Entrez module

    Returns:
        Tuple of (CDS accession, sequence) or None
    """
    try:
        # Use elink to find linked nucleotide - works for all protein prefixes
        # Try mRNA link first (best for XP_/NP_)
        handle = entrez.elink(
            dbfrom='protein',
            db='nuccore',
            id=accession,
            linkname='protein_nuccore_mrna'
        )
        result = entrez.read(handle)
        handle.close()

        if not result or not result[0].get('LinkSetDb'):
            # Try general nucleotide link (works for YP_, WP_, AP_)
            handle = entrez.elink(
                dbfrom='protein',
                db='nuccore',
                id=accession,
                linkname='protein_nuccore'
            )
            result = entrez.read(handle)
            handle.close()

        if result and result[0].get('LinkSetDb'):
            nuc_id = result[0]['LinkSetDb'][0]['Link'][0]['Id']

            # Fetch the CDS sequence
            handle = entrez.efetch(
                db='nuccore',
                id=nuc_id,
                rettype='fasta_cds_na',
                retmode='text'
            )
            cds_data = handle.read()
            handle.close()

            # Parse the FASTA
            if cds_data and '>' in cds_data:
                lines = cds_data.strip().split('\n')
                cds_header = lines[0]
                cds_seq = ''.join(lines[1:])

                # Extract accession from header
                cds_acc_match = re.search(r'lcl\|([^\s]+)', cds_header)
                if cds_acc_match:
                    cds_acc = cds_acc_match.group(1)
                else:
                    cds_acc = nuc_id

                return (cds_acc, cds_seq)

        return None

    except Exception as e:
        return None


def fetch_cds_ncbi_tsa(accession: str, entrez) -> Optional[Tuple[str, str]]:
    """
    Fetch sequence for TSA/WGS contig accession.

    These are already nucleotide sequences, so we fetch directly.

    Args:
        accession: TSA accession
        entrez: Entrez module

    Returns:
        Tuple of (accession, sequence) or None
    """
    try:
        handle = entrez.efetch(
            db='nuccore',
            id=accession,
            rettype='fasta',
            retmode='text'
        )
        fasta_data = handle.read()
        handle.close()

        if fasta_data and '>' in fasta_data:
            lines = fasta_data.strip().split('\n')
            seq = ''.join(lines[1:])
            return (accession, seq)

        return None

    except Exception:
        # Try nucleotide est database
        try:
            handle = entrez.efetch(
                db='nucest',
                id=accession,
                rettype='fasta',
                retmode='text'
            )
            fasta_data = handle.read()
            handle.close()

            if fasta_data and '>' in fasta_data:
                lines = fasta_data.strip().split('\n')
                seq = ''.join(lines[1:])
                return (accession, seq)

        except Exception:
            pass

        return None


def fetch_cds_ena(accession: str) -> Optional[Tuple[str, str]]:
    """
    Fetch sequence from ENA.

    Args:
        accession: ENA/DDBJ accession

    Returns:
        Tuple of (accession, sequence) or None
    """
    import requests

    try:
        url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}"
        response = requests.get(url, timeout=30)

        if response.status_code == 200 and '>' in response.text:
            lines = response.text.strip().split('\n')
            seq = ''.join(lines[1:])
            return (accession, seq)

        return None

    except Exception:
        return None


def fetch_cds_uniprot(accession: str) -> Optional[Tuple[str, str]]:
    """
    Fetch CDS for UniProt protein via cross-references.

    Args:
        accession: UniProt accession

    Returns:
        Tuple of (CDS accession, sequence) or None
    """
    import requests

    try:
        # Get UniProt entry with cross-references
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        response = requests.get(url, timeout=30)

        if response.status_code != 200:
            return None

        data = response.json()

        # Look for EMBL cross-references
        xrefs = data.get('uniProtKBCrossReferences', [])

        for xref in xrefs:
            if xref.get('database') == 'EMBL':
                properties = {p['key']: p['value'] for p in xref.get('properties', [])}

                # Get the nucleotide accession
                nuc_acc = xref.get('id')
                if nuc_acc:
                    # Try to fetch from ENA
                    result = fetch_cds_ena(nuc_acc)
                    if result:
                        return result

        return None

    except Exception:
        return None


def fetch_cds_for_accession(info: AccessionInfo, entrez) -> AccessionInfo:
    """
    Fetch CDS for a single accession.

    Args:
        info: AccessionInfo object
        entrez: Entrez module

    Returns:
        Updated AccessionInfo with CDS data or error
    """
    result = None

    try:
        if info.database == 'ncbi_refseq_protein':
            result = fetch_cds_ncbi_refseq(info.accession, entrez)

        elif info.database in ['ncbi_tsa', 'ncbi_wgs', 'ncbi_genbank_nuc']:
            result = fetch_cds_ncbi_tsa(info.accession, entrez)

        elif info.database in ['ena_ddbj']:
            result = fetch_cds_ena(info.accession)

        elif info.database == 'uniprot':
            result = fetch_cds_uniprot(info.accession)

        elif info.database == 'ncbi_refseq_genomic':
            # Genomic scaffolds don't have direct CDS - need gene coordinates
            info.status = 'skipped'
            info.error_message = 'Genomic scaffold - no direct CDS available'
            return info

        elif info.database == 'unknown':
            # Try NCBI first, then ENA
            result = fetch_cds_ncbi_tsa(info.accession, entrez)
            if not result:
                result = fetch_cds_ena(info.accession)

        if result:
            info.cds_accession, info.cds_sequence = result
            info.status = 'success'
        else:
            info.status = 'failed'
            info.error_message = f'Could not retrieve CDS from {info.database}'

    except Exception as e:
        info.status = 'error'
        info.error_message = str(e)

    return info


def process_fasta_file(fasta_file: Path, output_dir: Path, entrez,
                       rate_limit: float = 0.15,
                       verbose: bool = True) -> Dict[str, Any]:
    """
    Process a single FASTA file and generate corresponding CDS file.

    Args:
        fasta_file: Input protein FASTA file
        output_dir: Directory to write CDS FASTA
        entrez: Entrez module
        rate_limit: Seconds between API calls
        verbose: Print progress

    Returns:
        Statistics dictionary
    """
    stats = {
        'total': 0,
        'success': 0,
        'failed': 0,
        'skipped': 0,
        'by_database': defaultdict(int)
    }

    # Parse headers
    accessions = parse_fasta_headers(fasta_file)
    stats['total'] = len(accessions)

    if verbose:
        print(f"  Found {len(accessions)} sequences")

    # Fetch CDS for each
    cds_sequences = []

    for i, info in enumerate(accessions):
        stats['by_database'][info.database] += 1

        # Rate limiting
        time.sleep(rate_limit)

        info = fetch_cds_for_accession(info, entrez)

        if info.status == 'success':
            stats['success'] += 1
            cds_sequences.append((info.protein_id, info.cds_sequence))
        elif info.status == 'skipped':
            stats['skipped'] += 1
        else:
            stats['failed'] += 1

        if verbose and (i + 1) % 50 == 0:
            print(f"    Processed {i + 1}/{len(accessions)}")

    # Write CDS file
    if cds_sequences:
        output_file = output_dir / fasta_file.name.replace('.faa', '_cds.fna')
        output_file = output_file.with_suffix('.fna')

        with open(output_file, 'w') as f:
            for protein_id, cds_seq in cds_sequences:
                f.write(f">{protein_id}\n")
                # Write sequence in 60-character lines
                for i in range(0, len(cds_seq), 60):
                    f.write(cds_seq[i:i+60] + '\n')

        if verbose:
            print(f"  Wrote {len(cds_sequences)} CDS to {output_file.name}")

    return stats


def find_protein_files(ref_dir: Path) -> List[Path]:
    """Find all protein FASTA files in reference directory."""
    files = list(ref_dir.rglob('*.faa'))
    # Also check for .fa and .fasta
    files.extend(ref_dir.rglob('*.fa'))
    files.extend(ref_dir.rglob('*.fasta'))

    return sorted(set(files))


def main():
    parser = argparse.ArgumentParser(
        description="Fetch CDS sequences for reference protein files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all reference files
  python fetch_reference_cds.py references/nath_et_al

  # Process specific file
  python fetch_reference_cds.py references/nath_et_al/lse/gastropoda/Aplysia_californica.faa

  # Specify output directory
  python fetch_reference_cds.py references/nath_et_al -o references/nath_et_al_cds

  # Save failed accessions log
  python fetch_reference_cds.py references/nath_et_al --log-failed failed_accessions.csv
        """
    )

    parser.add_argument('input', type=Path,
                       help='Reference directory or single FASTA file')
    parser.add_argument('-o', '--output-dir', type=Path, default=None,
                       help='Output directory (default: same as input with _cds suffix on files)')
    parser.add_argument('--log-failed', type=Path, default=None,
                       help='Log failed accessions to CSV file')
    parser.add_argument('--rate-limit', type=float, default=0.12,
                       help='Seconds between API calls (default: 0.12 for ~8/sec)')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress progress output')

    args = parser.parse_args()

    # Validate input
    if not args.input.exists():
        print(f"ERROR: Input not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Setup Entrez
    try:
        entrez = setup_entrez()
        api_key, _ = get_ncbi_credentials()
        if not args.quiet:
            print(f"NCBI API configured (API key: {'yes' if api_key else 'no'})")
            if not api_key:
                print("  Note: Add NCBI_API_KEY to .env for faster requests (10/sec vs 3/sec)")
    except ImportError:
        print("ERROR: Biopython not installed. Run: pip install biopython", file=sys.stderr)
        sys.exit(1)

    # Find files to process
    if args.input.is_file():
        files = [args.input]
        output_dir = args.output_dir or args.input.parent
    else:
        files = find_protein_files(args.input)
        output_dir = args.output_dir or args.input

    if not files:
        print("ERROR: No protein FASTA files found", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print(f"\nProcessing {len(files)} protein files")
        print("-" * 60)

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process files
    all_stats = {
        'files': 0,
        'total_seqs': 0,
        'success': 0,
        'failed': 0,
        'skipped': 0
    }

    failed_accessions = []

    for fpath in files:
        if not args.quiet:
            print(f"\nProcessing: {fpath.name}")

        # Determine output subdirectory (preserve structure)
        if args.input.is_dir():
            rel_path = fpath.parent.relative_to(args.input)
            file_output_dir = output_dir / rel_path
        else:
            file_output_dir = output_dir

        file_output_dir.mkdir(parents=True, exist_ok=True)

        stats = process_fasta_file(
            fpath,
            file_output_dir,
            entrez,
            rate_limit=args.rate_limit,
            verbose=not args.quiet
        )

        all_stats['files'] += 1
        all_stats['total_seqs'] += stats['total']
        all_stats['success'] += stats['success']
        all_stats['failed'] += stats['failed']
        all_stats['skipped'] += stats['skipped']

    # Summary
    if not args.quiet:
        print("\n" + "=" * 60)
        print("Summary")
        print("=" * 60)
        print(f"Files processed: {all_stats['files']}")
        print(f"Total sequences: {all_stats['total_seqs']}")
        print(f"CDS retrieved:   {all_stats['success']}")
        print(f"Failed:          {all_stats['failed']}")
        print(f"Skipped:         {all_stats['skipped']}")

    # Log failed accessions if requested
    if args.log_failed and failed_accessions:
        with open(args.log_failed, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'accession', 'database', 'protein_id', 'error'
            ])
            writer.writeheader()
            for acc_info in failed_accessions:
                writer.writerow({
                    'accession': acc_info.accession,
                    'database': acc_info.database,
                    'protein_id': acc_info.protein_id,
                    'error': acc_info.error_message
                })

        if not args.quiet:
            print(f"\nFailed accessions logged to: {args.log_failed}")

    # Exit with error if many failures
    if all_stats['failed'] > all_stats['success']:
        sys.exit(1)


if __name__ == '__main__':
    main()
