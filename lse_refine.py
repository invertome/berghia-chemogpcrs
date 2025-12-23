#!/usr/bin/env python3
# lse_refine.py
# Purpose: Classify orthogroups into multilevel LSE datasets using NCBI taxonomy,
#          optional gene tree reconciliation, and synteny validation.
# Inputs: Orthogroup FASTA ($1), optional: ID map CSV ($2), species tree ($3),
#         gene tree ($4), synteny IDs ($5), output directory ($6)
# Outputs: LSE classification files (${output_dir}/lse_*.txt)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
import os
import argparse
import re
from pathlib import Path
from typing import List, Set, Optional, Dict
from ete3 import NCBITaxa, Tree
import pandas as pd


def parse_args():
    """Parse command line arguments with proper defaults."""
    parser = argparse.ArgumentParser(
        description="Classify orthogroups into multilevel LSE datasets"
    )
    parser.add_argument('og_file', help='Orthogroup FASTA file')
    parser.add_argument('--id-map', dest='id_map_file', default=None,
                       help='ID mapping CSV file')
    parser.add_argument('--species-tree', dest='species_tree_file', default=None,
                       help='Species tree file (Newick)')
    parser.add_argument('--gene-tree', dest='gene_tree_file', default=None,
                       help='Gene tree file (Newick)')
    parser.add_argument('--synteny-ids', dest='synteny_file', default=None,
                       help='Synteny IDs file')
    parser.add_argument('--output-dir', dest='output_dir', required=True,
                       help='Output directory')
    parser.add_argument('--exclude-refs', action='store_true',
                       help='Exclude reference sequences from LSE classification')
    parser.add_argument('--require-duplication', action='store_true',
                       help='Require gene tree duplication evidence for LSE')
    parser.add_argument('--require-synteny', action='store_true',
                       help='Require synteny evidence for LSE')
    return parser.parse_args()


# Configurable taxonomy IDs from environment (with NCBI defaults)
# These can be set in config.sh: export LSE_AEOLID_TAXID=54397
LSE_AEOLID_TAXID = int(os.getenv('LSE_AEOLID_TAXID', 54397))      # Aeolidida
LSE_NUDIBRANCH_TAXID = int(os.getenv('LSE_NUDIBRANCH_TAXID', 13843))  # Nudibranchia
LSE_GASTROPOD_TAXID = int(os.getenv('LSE_GASTROPOD_TAXID', 644))   # Gastropoda

# Use local taxonomy database if available (set via setup_databases.py)
LOCAL_DB_DIR = os.getenv('LOCAL_DB_DIR', '')
local_taxa_db = Path(LOCAL_DB_DIR) / "taxa.sqlite" if LOCAL_DB_DIR else None


def get_ncbi_taxa() -> NCBITaxa:
    """Initialize NCBITaxa with local database if available."""
    if local_taxa_db and local_taxa_db.exists():
        ncbi = NCBITaxa(dbfile=str(local_taxa_db))
        print(f"Using local taxonomy database: {local_taxa_db}", file=sys.stderr)
    else:
        ncbi = NCBITaxa()  # Will use default location or download if needed
        if LOCAL_DB_DIR:
            print(f"Warning: Local taxonomy database not found at {local_taxa_db}, using default", file=sys.stderr)
    return ncbi


# Lineage cache for performance
_lineage_cache: Dict[int, Optional[List[int]]] = {}
_rank_cache: Dict[int, str] = {}
_ncbi: Optional[NCBITaxa] = None


def get_ncbi() -> NCBITaxa:
    """Lazy initialization of NCBITaxa."""
    global _ncbi
    if _ncbi is None:
        _ncbi = get_ncbi_taxa()
    return _ncbi


def get_cached_lineage(taxid: int) -> Optional[List[int]]:
    """Get lineage with caching to avoid repeated NCBI queries."""
    if taxid not in _lineage_cache:
        try:
            _lineage_cache[taxid] = get_ncbi().get_lineage(taxid)
        except Exception:
            _lineage_cache[taxid] = None
    return _lineage_cache[taxid]


def get_cached_rank(taxid: int) -> str:
    """Get rank with caching."""
    if taxid not in _rank_cache:
        try:
            ranks = get_ncbi().get_rank([taxid])
            _rank_cache[taxid] = ranks.get(taxid, 'unknown')
        except Exception:
            _rank_cache[taxid] = 'unknown'
    return _rank_cache[taxid]


def extract_taxid_from_header(header: str, exclude_refs: bool = False) -> Optional[int]:
    """
    Extract taxonomy ID from FASTA header.

    Handles multiple formats:
    - ref_TAXID_N (reference sequences) -> TAXID (or None if exclude_refs)
    - TAXID_N (regular sequences) -> TAXID
    - TAXID (just the ID) -> TAXID

    Args:
        header: FASTA header (with or without >)
        exclude_refs: If True, return None for reference sequences

    Returns:
        Taxonomy ID as integer, or None if not found/excluded
    """
    # Remove > and get first field
    seq_id = header.lstrip('>').split()[0]

    # Check if this is a reference sequence
    if seq_id.startswith('ref_'):
        if exclude_refs:
            return None
        # Extract taxid from ref_TAXID_N format
        parts = seq_id.split('_')
        if len(parts) >= 2:
            taxid_str = parts[1]
        else:
            return None
    else:
        # Regular format: TAXID_N or just TAXID
        parts = seq_id.split('_')
        if not parts or not parts[0]:
            return None
        taxid_str = parts[0]

    # Convert to integer
    try:
        return int(taxid_str)
    except ValueError:
        return None


def extract_taxids_from_fasta(og_file: str, id_map: Optional[pd.DataFrame] = None,
                              exclude_refs: bool = False) -> List[int]:
    """
    Extract taxonomy IDs from orthogroup FASTA file.

    Args:
        og_file: Path to FASTA file
        id_map: Optional ID mapping DataFrame for metadata lookup
        exclude_refs: If True, exclude reference sequences

    Returns:
        List of taxonomy IDs
    """
    taxids = []

    with open(og_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                taxid = extract_taxid_from_header(line, exclude_refs)
                if taxid is not None:
                    taxids.append(taxid)

    return taxids


def check_duplication_in_gene_tree(gene_tree: Tree, taxids: List[int]) -> bool:
    """
    Check if gene tree shows evidence of duplication events.

    A duplication is inferred when a node has descendant leaves from the
    same taxon on both sides (paralogs).

    Args:
        gene_tree: ETE3 Tree object
        taxids: List of taxids in the orthogroup

    Returns:
        True if duplication is detected
    """
    taxid_count = pd.Series(taxids).value_counts()

    # First check: are there multiple copies from any taxon?
    if not any(count > 1 for count in taxid_count):
        return False

    # Check gene tree for duplication nodes
    for node in gene_tree.traverse():
        if not node.is_leaf() and len(node.children) >= 2:
            # Get taxids from each child subtree
            child_taxid_sets = []
            for child in node.children:
                child_taxids = set()
                for leaf in child.get_leaves():
                    taxid = extract_taxid_from_header(leaf.name)
                    if taxid:
                        child_taxids.add(taxid)
                child_taxid_sets.append(child_taxids)

            # Duplication: same taxid appears in multiple child subtrees
            if len(child_taxid_sets) >= 2:
                intersection = child_taxid_sets[0]
                for s in child_taxid_sets[1:]:
                    intersection = intersection & s
                if intersection:
                    return True

    return False


def classify_lse_level(taxids: List[int]) -> Optional[str]:
    """
    Determine the LSE taxonomic level based on common lineage.

    Args:
        taxids: List of taxonomy IDs

    Returns:
        Level name ('aeolids', 'nudibranchs', 'gastropods') or None
    """
    if not taxids:
        return None

    # Get lineages for all taxids
    lineages = []
    for taxid in taxids:
        lineage = get_cached_lineage(taxid)
        if lineage:
            lineages.append(set(lineage))

    if not lineages:
        return None

    # Find common ancestors
    common_lineage = set.intersection(*lineages)

    # Check taxonomic levels (most specific first)
    if LSE_AEOLID_TAXID in common_lineage:
        return 'aeolids'
    elif LSE_NUDIBRANCH_TAXID in common_lineage:
        return 'nudibranchs'
    elif LSE_GASTROPOD_TAXID in common_lineage:
        return 'gastropods'

    return None


def main():
    args = parse_args()

    og_file = args.og_file
    output_dir = args.output_dir

    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load optional inputs
    id_map = None
    if args.id_map_file and os.path.exists(args.id_map_file):
        id_map = pd.read_csv(args.id_map_file)

    gene_tree = None
    if args.gene_tree_file and os.path.exists(args.gene_tree_file):
        try:
            gene_tree = Tree(args.gene_tree_file, format=1)
        except Exception as e:
            print(f"Warning: Could not load gene tree: {e}", file=sys.stderr)

    synteny_ids = []
    if args.synteny_file and os.path.exists(args.synteny_file):
        synteny_ids = pd.read_csv(args.synteny_file, names=['id'])['id'].tolist()

    # Extract taxids from orthogroup
    taxids = extract_taxids_from_fasta(og_file, id_map, args.exclude_refs)

    if not taxids:
        print(f"Warning: No valid taxids found in {og_file}", file=sys.stderr)
        sys.exit(0)

    # Classify LSE level
    level = classify_lse_level(taxids)

    if not level:
        # Not an LSE at any defined level
        sys.exit(0)

    # Optional: Check for duplication evidence in gene tree
    if args.require_duplication and gene_tree:
        has_duplication = check_duplication_in_gene_tree(gene_tree, taxids)
        if not has_duplication:
            print(f"No duplication evidence in gene tree for {og_file}", file=sys.stderr)
            sys.exit(0)

    # Optional: Check for synteny evidence
    if args.require_synteny and synteny_ids:
        og_ids = []
        with open(og_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    og_ids.append(line.strip().lstrip('>').split()[0])

        synteny_overlap = any(oid in synteny_ids for oid in og_ids)
        if not synteny_overlap:
            print(f"No synteny evidence for {og_file}", file=sys.stderr)
            sys.exit(0)

    # Write classification result
    base = Path(og_file).stem
    output_file = Path(output_dir) / f"lse_{level}.txt"

    with open(output_file, 'a') as f:
        f.write(f"{base}\n")

    print(f"Classified {base} as LSE level: {level}")


if __name__ == '__main__':
    # Support both old positional args and new argparse style
    if len(sys.argv) >= 7 and not sys.argv[1].startswith('-'):
        # Legacy mode: positional arguments
        # Args: og_file id_map_file species_tree_file gene_tree_file synteny_file output_dir
        og_file = sys.argv[1]
        id_map_file = sys.argv[2]
        species_tree_file = sys.argv[3]
        gene_tree_file = sys.argv[4]
        synteny_file = sys.argv[5]
        output_dir = sys.argv[6]

        # Convert to argparse-style call, only adding args for non-empty file paths
        new_argv = [sys.argv[0], og_file]

        if id_map_file and os.path.exists(id_map_file):
            new_argv.extend(['--id-map', id_map_file])
        if species_tree_file and os.path.exists(species_tree_file):
            new_argv.extend(['--species-tree', species_tree_file])
        if gene_tree_file and os.path.exists(gene_tree_file):
            new_argv.extend(['--gene-tree', gene_tree_file])
        if synteny_file and os.path.exists(synteny_file):
            new_argv.extend(['--synteny-ids', synteny_file])

        new_argv.extend(['--output-dir', output_dir])
        # Note: --exclude-refs is typically desired for LSE classification
        new_argv.append('--exclude-refs')

        sys.argv = new_argv

    main()
