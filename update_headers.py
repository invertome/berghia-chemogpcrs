#!/usr/bin/env python3
# update_headers.py
# Purpose: Update FASTA headers with short IDs while preserving taxonomic information.
#          Generates a comprehensive metadata CSV for use by downstream scripts.
# Inputs: FASTA file ($1), output ID map CSV ($2), optional: --source-type, --append
# Outputs: Updated FASTA file (${fasta_file}_updated.fa), ID map CSV
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
import os
import re
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Track counters per taxid for unique IDs
taxid_counters = defaultdict(int)


def extract_taxid(header, description=""):
    """
    Extract taxonomic identifier from FASTA header.
    Tries multiple common formats:
    - taxid_proteinid
    - species_proteinid
    - >sp|ACCESSION|NAME OS=Species OX=taxid
    - Plain accession
    """
    full_header = f"{header} {description}"

    # Try to extract NCBI taxid from OX= field (UniProt format)
    ox_match = re.search(r'OX=(\d+)', full_header)
    if ox_match:
        return ox_match.group(1)

    # Try to extract taxid from underscore-separated format
    parts = header.split('_')
    if len(parts) >= 2:
        # Check if first part looks like a taxid (numeric or short species code)
        first_part = parts[0]
        if first_part.isdigit() or (len(first_part) <= 10 and first_part.isalnum()):
            return first_part

    # Try to extract from UniProt format OS=Species
    uniprot_match = re.search(r'OS=([A-Za-z]+)', full_header)
    if uniprot_match:
        return uniprot_match.group(1).lower()

    # Try to extract from NCBI format (gi|xxx|ref|xxx| or similar)
    if '|' in header:
        parts = header.split('|')
        for part in parts:
            if part and not part.isdigit():
                return part[:10]  # First 10 chars of non-numeric part

    # Fallback: use first 10 chars
    return header[:10] if len(header) > 10 else header


def extract_species_name(header, description=""):
    """Extract species name from header if available."""
    full_header = f"{header} {description}"

    # Try UniProt format: OS=Homo sapiens
    os_match = re.search(r'OS=([^=]+?)(?:\s+[A-Z]{2}=|$)', full_header)
    if os_match:
        return os_match.group(1).strip()

    # Try NCBI format: [Species name]
    bracket_match = re.search(r'\[([^\]]+)\]', full_header)
    if bracket_match:
        return bracket_match.group(1).strip()

    return ""


def extract_gene_name(header, description=""):
    """Extract gene name from header if available."""
    full_header = f"{header} {description}"

    # Try UniProt format: GN=GENE
    gn_match = re.search(r'GN=(\S+)', full_header)
    if gn_match:
        return gn_match.group(1)

    # Try common patterns: gene=XXX or gene:XXX
    gene_match = re.search(r'gene[=:](\S+)', full_header, re.IGNORECASE)
    if gene_match:
        return gene_match.group(1)

    return ""


def infer_source_type(fasta_file, header):
    """Infer source type from filename or header patterns."""
    basename = os.path.basename(fasta_file).lower()
    header_lower = header.lower()

    # Check filename patterns
    if 'reference' in basename or '_ref' in basename or 'conserved' in basename:
        return 'reference'
    if 'outgroup' in basename:
        return 'outgroup'
    if 'target' in basename or 'query' in basename:
        return 'target'

    # Check header patterns
    if header.startswith('ref_'):
        return 'reference'
    if 'outgroup' in header_lower:
        return 'outgroup'

    return 'unknown'


def process_sequences(fasta_file, id_map_file, source_type=None, append=False,
                      id_prefix=None, existing_counters=None):
    """
    Process sequences using Biopython for robust FASTA handling.

    Args:
        fasta_file: Path to input FASTA file
        id_map_file: Path to output metadata CSV
        source_type: Override source type (reference/target/outgroup)
        append: If True, append to existing CSV instead of overwriting
        id_prefix: Custom prefix for short IDs (default: 'ref' for reference, taxid for others)
        existing_counters: Dictionary of existing taxid counters (for append mode)

    Returns:
        Tuple of (count, taxid_counters)
    """
    global taxid_counters

    if existing_counters:
        taxid_counters = existing_counters

    output_file = f"{fasta_file}_updated.fa"
    updated_records = []

    mode = 'a' if append else 'w'
    write_header = not append or not os.path.exists(id_map_file)

    with open(id_map_file, mode) as map_out:
        if write_header:
            map_out.write("original_id,short_id,taxid,source_type,species_name,gene_name,seq_length\n")

        # Use SeqIO for robust FASTA parsing (handles multi-line sequences, validates format)
        for record in SeqIO.parse(fasta_file, "fasta"):
            original_id = record.id
            description = record.description.replace(record.id, "").strip()

            taxid = extract_taxid(original_id, description)
            species = extract_species_name(original_id, description)
            gene = extract_gene_name(original_id, description)
            src_type = source_type or infer_source_type(fasta_file, original_id)

            # Increment counter for this taxid
            taxid_counters[taxid] += 1
            counter = taxid_counters[taxid]

            # Create short ID preserving taxonomy
            if id_prefix:
                short_id = f"{id_prefix}_{taxid}_{counter}"
            elif src_type == 'reference':
                short_id = f"ref_{taxid}_{counter}"
            else:
                short_id = f"{taxid}_{counter}"

            # Escape commas in fields for CSV
            species_safe = species.replace(',', ';')
            gene_safe = gene.replace(',', ';')
            original_safe = original_id.replace(',', ';')

            # Write mapping with extended metadata
            map_out.write(f"{original_safe},{short_id},{taxid},{src_type},{species_safe},{gene_safe},{len(record.seq)}\n")

            # Create updated record with new ID
            updated_record = SeqRecord(
                record.seq,
                id=short_id,
                description=""  # Clear description to keep header clean
            )
            updated_records.append(updated_record)

    # Write updated sequences
    SeqIO.write(updated_records, output_file, "fasta")

    return len(updated_records), dict(taxid_counters)


def main():
    parser = argparse.ArgumentParser(
        description='Update FASTA headers with short IDs and generate metadata CSV'
    )
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('id_map_file', help='Output metadata CSV file')
    parser.add_argument('--source-type', choices=['reference', 'target', 'outgroup'],
                        help='Override source type detection')
    parser.add_argument('--append', action='store_true',
                        help='Append to existing CSV instead of overwriting')
    parser.add_argument('--id-prefix', help='Custom prefix for short IDs')

    # Support legacy positional-only usage
    if len(sys.argv) == 3 and not sys.argv[1].startswith('-'):
        fasta_file = sys.argv[1]
        id_map_file = sys.argv[2]
        count, counters = process_sequences(fasta_file, id_map_file)
    else:
        args = parser.parse_args()
        count, counters = process_sequences(
            args.fasta_file,
            args.id_map_file,
            source_type=args.source_type,
            append=args.append,
            id_prefix=args.id_prefix
        )

    print(f"Processed {count} sequences from {len(counters)} taxa")


if __name__ == '__main__':
    main()
