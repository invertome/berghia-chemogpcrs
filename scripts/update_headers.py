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

# Opt-in gate for assigning IDs that do not already exist in the map.
# IDs are WRITE-ONCE: once a record has a short ID, that ID denotes that record
# forever. Re-issuing one is silent identity corruption -- both the old and new
# files contain the same ID, so every downstream join SUCCEEDS and returns the
# wrong record. The default is therefore refusal: when a map already exists,
# any record it does not cover is a hard error unless the operator explicitly
# opts in.
ALLOW_NEW_IDS_ENV = 'ALLOW_NEW_IDS'

# A short ID this script previously assigned: ref_<taxid>_<N>. Used to detect
# input that has ALREADY been through the minter (stage 01 moves its output
# over its own input, so a re-run feeds short IDs back in). Re-deriving a taxid
# from such a header yields 'ref' for every taxon and collapses the whole
# namespace, so it is never allowed to fall through to assignment.
MINTED_ID_RE = re.compile(r'^ref_[A-Za-z0-9]+_\d+$')


class IdentityError(Exception):
    """Raised when honoring a request would re-issue or invent an identifier."""


def _sample(items, limit=5):
    """Render a short, actionable excerpt of an offending-ID list."""
    shown = ", ".join(items[:limit])
    if len(items) > limit:
        shown += f", ... (+{len(items) - limit} more)"
    return shown


def new_ids_allowed():
    """True when the operator has explicitly opted in to new-ID assignment."""
    return os.environ.get(ALLOW_NEW_IDS_ENV, '').strip().lower() in {
        '1', 'true', 'yes',
    }


def load_id_map(id_map_file):
    """Load an existing ID map.

    Returns (by_original, by_short, high_water) where by_original maps source
    accession -> short ID, by_short is the inverse, and high_water records the
    largest suffix ever handed out per taxid. The high-water mark is what keeps
    a retired ID retired: pruning a record leaves a GAP, and the gap is never
    reused by a later record.
    """
    by_original = {}
    by_short = {}
    high_water = defaultdict(int)

    if not os.path.exists(id_map_file):
        return by_original, by_short, high_water

    with open(id_map_file) as handle:
        handle.readline()  # discard the CSV header row
        for line in handle:
            fields = line.rstrip('\n').split(',')
            if len(fields) < 3 or not fields[0]:
                continue
            original, short, taxid = fields[0], fields[1], fields[2]
            by_original[original] = short
            by_short[short] = original
            suffix = re.search(r'_(\d+)$', short)
            if suffix:
                high_water[taxid] = max(high_water[taxid], int(suffix.group(1)))

    return by_original, by_short, high_water


def extract_taxid(header, description=""):
    """
    Extract taxonomic identifier from FASTA header.
    Tries multiple common formats:
    - taxid_proteinid (numeric NCBI taxid)
    - species_code_accession (e.g., aplcal_XP_005089826.1, hadi_HDSC00048_47552)
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
        first_part = parts[0]
        # Numeric taxid (e.g., 1287507_Berghia_stephanieae)
        if first_part.isdigit():
            return first_part
        # Short alphanumeric species code (e.g., aplcal, hadi, gima)
        # Accept codes 2-10 chars, all lowercase alpha or with digits
        if 2 <= len(first_part) <= 10 and first_part.isalnum():
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

    output_file = f"{fasta_file}_updated.fa"

    # --- Reuse before assign -------------------------------------------------
    # Look for an existing map FIRST. If one exists it is authoritative: every
    # record must resolve through it, and the file is only ever appended to,
    # never rewritten. That is what makes a re-run a genuine no-op instead of a
    # renumbering.
    existing_by_original, existing_by_short, high_water = load_id_map(id_map_file)
    map_exists = bool(existing_by_original)

    taxid_counters = defaultdict(int, high_water)
    if existing_counters:
        for key, value in existing_counters.items():
            taxid_counters[key] = max(taxid_counters[key], value)

    # With no prior map there is no identity to preserve, so the first run
    # bootstraps freely. Once a map exists, assignment needs the opt-in.
    may_assign = (not map_exists) or new_ids_allowed()

    updated_records = []
    new_rows = []
    unresolved = []
    premapped_without_entry = []

    # Use SeqIO for robust FASTA parsing (handles multi-line sequences, validates format)
    for record in SeqIO.parse(fasta_file, "fasta"):
        original_id = record.id
        description = record.description.replace(record.id, "").strip()

        # 1. The record's source accession is already mapped -> reuse its ID.
        #    This is what keeps the partial CDS subset aligned with the protein
        #    pass: resolution is by IDENTITY, not by position, so a dropped CDS
        #    record cannot shift the ones after it.
        short_id = existing_by_original.get(original_id)

        # 2. The input is itself already-mapped output (stage 01 moves the
        #    updated FASTA over its own input). Carry the ID through unchanged
        #    and add no row -- the record already owns this ID.
        if short_id is None and original_id in existing_by_short:
            short_id = original_id

        if short_id is not None:
            updated_records.append(
                SeqRecord(record.seq, id=short_id, description="")
            )
            continue

        # 3. Unresolved. An already-minted-looking header that the map does not
        #    know about means the map is missing or stale relative to this
        #    FASTA. Deriving a taxid from it would yield 'ref' for every taxon,
        #    so this is a hard stop that the opt-in deliberately cannot waive.
        if MINTED_ID_RE.match(original_id):
            premapped_without_entry.append(original_id)
            continue

        if not may_assign:
            unresolved.append(original_id)
            continue

        taxid = extract_taxid(original_id, description)
        species = extract_species_name(original_id, description)
        gene = extract_gene_name(original_id, description)
        src_type = source_type or infer_source_type(fasta_file, original_id)

        # Continue from the high-water mark so a gap left by a pruned record is
        # never handed to a different record.
        taxid_counters[taxid] += 1
        position = taxid_counters[taxid]

        # Create short ID preserving taxonomy
        # Avoid double-prefixing (e.g., ref_ref_aplcal_1)
        if id_prefix:
            short_id = f"{id_prefix}_{taxid}_{position}"
        elif src_type == 'reference':
            if taxid == 'ref' or taxid.startswith('ref_'):
                short_id = f"{taxid}_{position}"
            else:
                short_id = f"ref_{taxid}_{position}"
        else:
            short_id = f"{taxid}_{position}"

        # Escape commas in fields for CSV
        species_safe = species.replace(',', ';')
        gene_safe = gene.replace(',', ';')
        original_safe = original_id.replace(',', ';')

        new_rows.append(
            f"{original_safe},{short_id},{taxid},{src_type},"
            f"{species_safe},{gene_safe},{len(record.seq)}\n"
        )
        existing_by_short[short_id] = original_id
        updated_records.append(
            SeqRecord(record.seq, id=short_id, description="")
        )

    # --- Abort loudly, before writing anything -------------------------------
    if premapped_without_entry:
        raise IdentityError(
            "input contains short IDs this map does not cover: "
            + _sample(premapped_without_entry)
            + f"\nThe FASTA has already been through {os.path.basename(__file__)}"
            f" but '{id_map_file}' is missing or stale."
            "\nRestore the ID map that assigned those IDs; re-deriving them"
            " would give every taxon the same 'ref' namespace and silently"
            " re-point existing IDs at different records."
        )

    if unresolved:
        raise IdentityError(
            f"{len(unresolved)} record(s) are absent from '{id_map_file}': "
            + _sample(unresolved)
            + "\nIDs are write-once, so this script will not assign new ones by"
            f" default. Set {ALLOW_NEW_IDS_ENV}=1 to extend the map (existing"
            " IDs are preserved; new records continue above the high-water"
            " mark)."
        )

    # --- Commit --------------------------------------------------------------
    # The map is append-only: a pruned record keeps its row so its ID stays
    # retired rather than being recycled. `append` is retained for callers but
    # no longer decides truncation -- nothing truncates the map.
    if new_rows or not os.path.exists(id_map_file):
        need_header = not os.path.exists(id_map_file)
        with open(id_map_file, 'a') as map_out:
            if need_header:
                map_out.write(
                    "original_id,short_id,taxid,source_type,"
                    "species_name,gene_name,seq_length\n"
                )
            map_out.writelines(new_rows)

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
    try:
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
    except IdentityError as exc:
        # Fail hard and visibly. Falling back to assigning IDs here is exactly
        # the silent corruption this guard exists to prevent.
        sys.stderr.write(f"ERROR: refusing to re-issue identifiers.\n{exc}\n")
        return 2

    print(f"Processed {count} sequences from {len(counters)} taxa")
    return 0


if __name__ == '__main__':
    sys.exit(main())
