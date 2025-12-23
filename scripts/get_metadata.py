#!/usr/bin/env python3
"""
get_metadata.py
Purpose: Centralized metadata lookup utility for sequence information.
         Provides a single source of truth for taxid, source type, and other
         sequence metadata, eliminating fragile FASTA header parsing.

Usage:
    # As a module
    from get_metadata import MetadataLookup
    lookup = MetadataLookup('path/to/id_map.csv')
    taxid = lookup.get_taxid('ref_9606_1')

    # As a CLI tool
    python3 get_metadata.py --metadata id_map.csv --seq-id ref_9606_1
    python3 get_metadata.py --metadata id_map.csv --fasta orthogroup.fa --field taxid
    python3 get_metadata.py --metadata id_map.csv --fasta orthogroup.fa --unique-taxids

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
"""

import sys
import os
import csv
import argparse
from typing import Dict, List, Set, Optional, Any
from functools import lru_cache


class MetadataLookup:
    """
    Centralized metadata lookup for sequence IDs.

    Provides efficient lookup of sequence metadata from the ID map CSV,
    with caching for repeated lookups and support for both exact and
    fuzzy matching.
    """

    def __init__(self, metadata_file: str):
        """
        Initialize the metadata lookup.

        Args:
            metadata_file: Path to the metadata CSV file
        """
        self.metadata_file = metadata_file
        self._metadata: Dict[str, Dict[str, Any]] = {}
        self._taxid_to_seqs: Dict[str, List[str]] = {}
        self._loaded = False

    def _load(self):
        """Load metadata from CSV file."""
        if self._loaded:
            return

        if not os.path.exists(self.metadata_file):
            raise FileNotFoundError(f"Metadata file not found: {self.metadata_file}")

        with open(self.metadata_file, 'r') as f:
            reader = csv.DictReader(f)

            # Validate required headers exist
            required_headers = {'short_id'}  # Minimum required
            optional_headers = {'original_id', 'taxid'}
            if reader.fieldnames is None:
                raise ValueError(f"Metadata file is empty or has no header: {self.metadata_file}")

            headers = set(reader.fieldnames)
            missing_required = required_headers - headers
            if missing_required:
                raise ValueError(f"Missing required headers in {self.metadata_file}: {missing_required}")

            missing_optional = optional_headers - headers
            if missing_optional:
                import sys
                print(f"Warning: Optional headers missing in metadata file: {missing_optional}", file=sys.stderr)

            for row in reader:
                short_id = row.get('short_id', '')
                if short_id:
                    self._metadata[short_id] = row

                    # Index by original_id too for reverse lookups
                    original_id = row.get('original_id', '')
                    if original_id and original_id != short_id:
                        self._metadata[original_id] = row

                    # Build taxid -> sequences index
                    taxid = row.get('taxid', '')
                    if taxid:
                        if taxid not in self._taxid_to_seqs:
                            self._taxid_to_seqs[taxid] = []
                        self._taxid_to_seqs[taxid].append(short_id)

        self._loaded = True

    def get(self, seq_id: str, field: str, default: str = '') -> str:
        """
        Get a metadata field for a sequence ID.

        Args:
            seq_id: Sequence ID (short_id or original_id)
            field: Field name to retrieve
            default: Default value if not found

        Returns:
            Field value or default
        """
        self._load()
        record = self._metadata.get(seq_id, {})
        return record.get(field, default)

    def get_taxid(self, seq_id: str) -> str:
        """Get taxid for a sequence ID."""
        return self.get(seq_id, 'taxid', '')

    def get_source_type(self, seq_id: str) -> str:
        """Get source type (reference/target/outgroup) for a sequence ID."""
        return self.get(seq_id, 'source_type', 'unknown')

    def get_species_name(self, seq_id: str) -> str:
        """Get species name for a sequence ID."""
        return self.get(seq_id, 'species_name', '')

    def get_gene_name(self, seq_id: str) -> str:
        """Get gene name for a sequence ID."""
        return self.get(seq_id, 'gene_name', '')

    def get_original_id(self, seq_id: str) -> str:
        """Get original ID for a short ID."""
        return self.get(seq_id, 'original_id', seq_id)

    def is_reference(self, seq_id: str) -> bool:
        """Check if sequence is from reference database."""
        return self.get_source_type(seq_id) == 'reference'

    def is_target(self, seq_id: str) -> bool:
        """Check if sequence is a target (query) sequence."""
        return self.get_source_type(seq_id) == 'target'

    def get_sequences_for_taxid(self, taxid: str) -> List[str]:
        """Get all sequence IDs for a given taxid."""
        self._load()
        return self._taxid_to_seqs.get(taxid, [])

    def get_all_taxids(self) -> Set[str]:
        """Get all unique taxids in the metadata."""
        self._load()
        return set(self._taxid_to_seqs.keys())

    def get_taxids_from_fasta(self, fasta_file: str, exclude_references: bool = False) -> Set[str]:
        """
        Get unique taxids from sequences in a FASTA file.

        Args:
            fasta_file: Path to FASTA file
            exclude_references: If True, exclude reference sequences

        Returns:
            Set of taxids
        """
        self._load()
        taxids = set()

        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    seq_id = line[1:].split()[0].strip()
                    if exclude_references and self.is_reference(seq_id):
                        continue
                    taxid = self.get_taxid(seq_id)
                    if taxid:
                        taxids.add(taxid)

        return taxids

    def get_non_reference_taxids_from_fasta(self, fasta_file: str) -> Set[str]:
        """Get taxids from non-reference sequences in a FASTA file."""
        return self.get_taxids_from_fasta(fasta_file, exclude_references=True)

    def fallback_taxid_from_header(self, header: str) -> str:
        """
        Fallback method to extract taxid from header if not in metadata.

        Handles the ref_TAXID_N format and plain TAXID_N format.

        Args:
            header: FASTA header (without '>')

        Returns:
            Extracted taxid or empty string
        """
        parts = header.split('_')
        if len(parts) >= 2:
            if parts[0] == 'ref':
                # ref_TAXID_N format
                return parts[1]
            else:
                # TAXID_N format
                return parts[0]
        return ''

    def get_taxid_robust(self, seq_id: str) -> str:
        """
        Get taxid with fallback to header parsing.

        First tries metadata lookup, then falls back to header parsing.
        This provides backwards compatibility during migration.

        Args:
            seq_id: Sequence ID

        Returns:
            Taxid
        """
        taxid = self.get_taxid(seq_id)
        if taxid:
            return taxid
        return self.fallback_taxid_from_header(seq_id)


# Global singleton for convenience
_default_lookup: Optional[MetadataLookup] = None


def init_metadata(metadata_file: str) -> MetadataLookup:
    """Initialize the global metadata lookup."""
    global _default_lookup
    _default_lookup = MetadataLookup(metadata_file)
    return _default_lookup


def get_lookup() -> Optional[MetadataLookup]:
    """Get the global metadata lookup instance."""
    return _default_lookup


def get_taxid(seq_id: str) -> str:
    """Convenience function to get taxid using global lookup."""
    if _default_lookup is None:
        raise RuntimeError("Metadata not initialized. Call init_metadata() first.")
    return _default_lookup.get_taxid_robust(seq_id)


def main():
    parser = argparse.ArgumentParser(
        description='Look up sequence metadata from ID map CSV'
    )
    parser.add_argument('--metadata', '-m', required=True,
                        help='Path to metadata CSV file')
    parser.add_argument('--seq-id', '-s',
                        help='Sequence ID to look up')
    parser.add_argument('--fasta', '-f',
                        help='FASTA file to extract IDs from')
    parser.add_argument('--field', '-F', default='taxid',
                        choices=['taxid', 'source_type', 'species_name', 'gene_name', 'original_id', 'seq_length'],
                        help='Field to retrieve (default: taxid)')
    parser.add_argument('--unique-taxids', action='store_true',
                        help='Output unique taxids from FASTA (space-separated)')
    parser.add_argument('--exclude-refs', action='store_true',
                        help='Exclude reference sequences when processing FASTA')
    parser.add_argument('--json', action='store_true',
                        help='Output as JSON')

    args = parser.parse_args()

    try:
        lookup = MetadataLookup(args.metadata)

        if args.seq_id:
            # Single sequence lookup
            if args.json:
                import json
                lookup._load()
                record = lookup._metadata.get(args.seq_id, {})
                print(json.dumps(record))
            else:
                value = lookup.get(args.seq_id, args.field)
                print(value)

        elif args.fasta:
            if args.unique_taxids:
                # Get unique taxids from FASTA
                taxids = lookup.get_taxids_from_fasta(args.fasta, exclude_references=args.exclude_refs)
                print(' '.join(sorted(taxids)))
            else:
                # Get field for each sequence in FASTA
                with open(args.fasta, 'r') as f:
                    for line in f:
                        if line.startswith('>'):
                            seq_id = line[1:].split()[0].strip()
                            if args.exclude_refs and lookup.is_reference(seq_id):
                                continue
                            value = lookup.get(seq_id, args.field)
                            print(f"{seq_id}\t{value}")

        else:
            parser.error("Either --seq-id or --fasta is required")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
