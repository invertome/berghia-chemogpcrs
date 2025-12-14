#!/usr/bin/env python3
# update_headers.py
# Purpose: Update FASTA headers with short IDs while preserving taxonomic information.
# Inputs: FASTA file ($1), output ID map CSV ($2)
# Outputs: Updated FASTA file (${fasta_file}_updated.fa), ID map CSV
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

fasta_file = sys.argv[1]
id_map_file = sys.argv[2]
output_file = f"{fasta_file}_updated.fa"

# Track counters per taxid for unique IDs
taxid_counters = defaultdict(int)


def extract_taxid(header):
    """
    Extract taxonomic identifier from FASTA header.
    Tries multiple common formats:
    - taxid_proteinid
    - species_proteinid
    - >sp|ACCESSION|NAME OS=Species
    - Plain accession
    """
    # Try to extract taxid from underscore-separated format
    parts = header.split('_')
    if len(parts) >= 2:
        # First part is likely taxid or species abbreviation
        return parts[0]

    # Try to extract from UniProt format
    uniprot_match = re.search(r'OS=([A-Za-z]+)', header)
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


def process_sequences():
    """Process sequences using Biopython for robust FASTA handling."""
    updated_records = []

    with open(id_map_file, 'w') as map_out:
        map_out.write("original_id,short_id,taxid\n")

        # Use SeqIO for robust FASTA parsing (handles multi-line sequences, validates format)
        for record in SeqIO.parse(fasta_file, "fasta"):
            original_id = record.id
            taxid = extract_taxid(original_id)

            # Increment counter for this taxid
            taxid_counters[taxid] += 1
            counter = taxid_counters[taxid]

            # Create short ID preserving taxonomy: ref_TAXID_N
            short_id = f"ref_{taxid}_{counter}"

            # Write mapping
            map_out.write(f"{original_id},{short_id},{taxid}\n")

            # Create updated record with new ID
            updated_record = SeqRecord(
                record.seq,
                id=short_id,
                description=""  # Clear description to keep header clean
            )
            updated_records.append(updated_record)

    # Write updated sequences
    SeqIO.write(updated_records, output_file, "fasta")

    return len(updated_records)


if __name__ == '__main__':
    count = process_sequences()
    print(f"Processed {count} sequences from {len(taxid_counters)} taxa")
