#!/usr/bin/env python3
# convergent_evolution.py
# Purpose: Detect convergent amino acid substitutions across independent GPCR lineages.
# Inputs: Multiple sequence alignment, phylogenetic tree
# Outputs: Convergent substitution sites with functional annotations
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   python convergent_evolution.py <alignment> <tree> <output_prefix> [focal_species]
#
# Arguments:
#   alignment     - Multiple sequence alignment (FASTA format)
#   tree          - Phylogenetic tree (Newick format)
#   output_prefix - Prefix for output files
#   focal_species - Optional: comma-separated list of focal species to compare

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from ete3 import Tree

# Amino acid properties for functional interpretation
AA_PROPERTIES = {
    # Hydrophobic
    'A': {'hydrophobic': True, 'polar': False, 'charged': False, 'aromatic': False, 'small': True},
    'V': {'hydrophobic': True, 'polar': False, 'charged': False, 'aromatic': False, 'small': False},
    'L': {'hydrophobic': True, 'polar': False, 'charged': False, 'aromatic': False, 'small': False},
    'I': {'hydrophobic': True, 'polar': False, 'charged': False, 'aromatic': False, 'small': False},
    'M': {'hydrophobic': True, 'polar': False, 'charged': False, 'aromatic': False, 'small': False},
    'F': {'hydrophobic': True, 'polar': False, 'charged': False, 'aromatic': True, 'small': False},
    'W': {'hydrophobic': True, 'polar': False, 'charged': False, 'aromatic': True, 'small': False},
    'P': {'hydrophobic': True, 'polar': False, 'charged': False, 'aromatic': False, 'small': False},
    # Polar uncharged
    'S': {'hydrophobic': False, 'polar': True, 'charged': False, 'aromatic': False, 'small': True},
    'T': {'hydrophobic': False, 'polar': True, 'charged': False, 'aromatic': False, 'small': False},
    'N': {'hydrophobic': False, 'polar': True, 'charged': False, 'aromatic': False, 'small': False},
    'Q': {'hydrophobic': False, 'polar': True, 'charged': False, 'aromatic': False, 'small': False},
    'Y': {'hydrophobic': False, 'polar': True, 'charged': False, 'aromatic': True, 'small': False},
    'C': {'hydrophobic': False, 'polar': True, 'charged': False, 'aromatic': False, 'small': True},
    'G': {'hydrophobic': False, 'polar': False, 'charged': False, 'aromatic': False, 'small': True},
    # Positively charged
    'K': {'hydrophobic': False, 'polar': True, 'charged': True, 'aromatic': False, 'small': False},
    'R': {'hydrophobic': False, 'polar': True, 'charged': True, 'aromatic': False, 'small': False},
    'H': {'hydrophobic': False, 'polar': True, 'charged': True, 'aromatic': True, 'small': False},
    # Negatively charged
    'D': {'hydrophobic': False, 'polar': True, 'charged': True, 'aromatic': False, 'small': False},
    'E': {'hydrophobic': False, 'polar': True, 'charged': True, 'aromatic': False, 'small': False},
}

# GPCR functional motifs (positions are approximate, relative to alignment)
GPCR_MOTIFS = {
    'DRY': {'pattern': '[DE]R[YF]', 'function': 'G-protein coupling', 'region': 'TM3-ICL2'},
    'NPxxY': {'pattern': 'NP..Y', 'function': 'Activation switch', 'region': 'TM7'},
    'CWxP': {'pattern': 'CW.P', 'function': 'Rotamer toggle', 'region': 'TM6'},
    'ERC': {'pattern': 'ERC', 'function': 'Mollusc DRY variant', 'region': 'TM3-ICL2'},
}


def load_alignment(alignment_file: str) -> MultipleSeqAlignment:
    """Load multiple sequence alignment from FASTA file."""
    alignment = AlignIO.read(alignment_file, "fasta")
    print(f"Loaded alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} positions",
          file=sys.stderr)
    return alignment


def load_tree(tree_file: str) -> Tree:
    """Load phylogenetic tree from Newick file."""
    tree = Tree(tree_file, format=1)
    print(f"Loaded tree: {len(tree)} leaves", file=sys.stderr)
    return tree


def get_ancestral_sequences(tree: Tree, alignment: MultipleSeqAlignment) -> Dict[str, str]:
    """
    Estimate ancestral sequences using parsimony.

    Simple implementation - for production use FastML or PAML.

    Args:
        tree: Phylogenetic tree
        alignment: Multiple sequence alignment

    Returns:
        Dictionary mapping node names to ancestral sequences
    """
    # Map leaf names to sequences
    seq_dict = {rec.id: str(rec.seq) for rec in alignment}
    aln_length = alignment.get_alignment_length()

    ancestral = {}

    # Traverse tree from leaves to root (post-order)
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            if node.name in seq_dict:
                ancestral[node.name] = seq_dict[node.name]
            continue

        # Get child sequences
        child_seqs = []
        for child in node.children:
            child_name = child.name if child.name else f"node_{id(child)}"
            if child_name in ancestral:
                child_seqs.append(ancestral[child_name])

        if not child_seqs:
            continue

        # Simple majority-rule consensus for ancestral state
        anc_seq = []
        for pos in range(aln_length):
            aas = [seq[pos] for seq in child_seqs if pos < len(seq)]
            if aas:
                # Most common amino acid
                aa_counts = defaultdict(int)
                for aa in aas:
                    if aa != '-':
                        aa_counts[aa] += 1
                if aa_counts:
                    anc_seq.append(max(aa_counts, key=aa_counts.get))
                else:
                    anc_seq.append('-')
            else:
                anc_seq.append('-')

        node_name = node.name if node.name else f"node_{id(node)}"
        ancestral[node_name] = ''.join(anc_seq)

    return ancestral


def identify_convergent_sites(alignment: MultipleSeqAlignment,
                              tree: Tree,
                              focal_species: Optional[List[str]] = None) -> List[Dict]:
    """
    Identify sites with convergent amino acid substitutions.

    Convergent sites are where the same amino acid appears in
    phylogenetically independent lineages.

    Args:
        alignment: Multiple sequence alignment
        tree: Phylogenetic tree
        focal_species: Optional list of species to focus on

    Returns:
        List of convergent site information
    """
    seq_dict = {rec.id: str(rec.seq) for rec in alignment}
    aln_length = alignment.get_alignment_length()

    # Get ancestral sequences
    print("Estimating ancestral sequences...", file=sys.stderr)
    ancestral = get_ancestral_sequences(tree, alignment)

    # Identify independent lineages
    leaf_names = [leaf.name for leaf in tree]

    if focal_species:
        # Compare focal species against others
        focal_leaves = [l for l in leaf_names if any(f in l for f in focal_species)]
        other_leaves = [l for l in leaf_names if l not in focal_leaves]
    else:
        # Compare first half vs second half (simple split)
        mid = len(leaf_names) // 2
        focal_leaves = leaf_names[:mid]
        other_leaves = leaf_names[mid:]

    convergent_sites = []

    for pos in range(aln_length):
        # Get amino acids at this position
        focal_aas = set()
        other_aas = set()

        for leaf in focal_leaves:
            if leaf in seq_dict and pos < len(seq_dict[leaf]):
                aa = seq_dict[leaf][pos]
                if aa != '-' and aa != 'X':
                    focal_aas.add(aa)

        for leaf in other_leaves:
            if leaf in seq_dict and pos < len(seq_dict[leaf]):
                aa = seq_dict[leaf][pos]
                if aa != '-' and aa != 'X':
                    other_aas.add(aa)

        # Check for convergent amino acids
        convergent_aas = focal_aas & other_aas

        if convergent_aas and len(focal_aas) > 0 and len(other_aas) > 0:
            # Check if this is a derived state (different from ancestral)
            root_aa = None
            if tree.name and tree.name in ancestral:
                root_seq = ancestral[tree.name]
                if pos < len(root_seq):
                    root_aa = root_seq[pos]

            for conv_aa in convergent_aas:
                # Calculate frequencies
                focal_freq = sum(1 for l in focal_leaves
                                if l in seq_dict and pos < len(seq_dict[l])
                                and seq_dict[l][pos] == conv_aa) / len(focal_leaves)
                other_freq = sum(1 for l in other_leaves
                                if l in seq_dict and pos < len(seq_dict[l])
                                and seq_dict[l][pos] == conv_aa) / len(other_leaves)

                # Only report if not ancestral and reasonably frequent
                is_derived = root_aa is None or conv_aa != root_aa
                is_significant = focal_freq > 0.1 and other_freq > 0.1

                if is_derived and is_significant:
                    convergent_sites.append({
                        'position': pos + 1,  # 1-indexed
                        'convergent_aa': conv_aa,
                        'ancestral_aa': root_aa,
                        'focal_frequency': focal_freq,
                        'other_frequency': other_freq,
                        'aa_properties': AA_PROPERTIES.get(conv_aa, {}),
                        'property_change': get_property_change(root_aa, conv_aa) if root_aa else None
                    })

    return convergent_sites


def get_property_change(aa1: str, aa2: str) -> Dict:
    """Determine what physicochemical properties changed between amino acids."""
    props1 = AA_PROPERTIES.get(aa1, {})
    props2 = AA_PROPERTIES.get(aa2, {})

    changes = {}
    for prop in ['hydrophobic', 'polar', 'charged', 'aromatic', 'small']:
        if props1.get(prop) != props2.get(prop):
            changes[prop] = f"{props1.get(prop, '?')} -> {props2.get(prop, '?')}"

    return changes


def annotate_functional_regions(convergent_sites: List[Dict],
                               alignment: MultipleSeqAlignment) -> List[Dict]:
    """
    Annotate convergent sites with potential functional significance.

    Args:
        convergent_sites: List of convergent site dictionaries
        alignment: Multiple sequence alignment

    Returns:
        Annotated convergent sites
    """
    import re

    # Get consensus sequence for motif searching
    consensus = []
    for pos in range(alignment.get_alignment_length()):
        aas = [str(rec.seq)[pos] for rec in alignment if str(rec.seq)[pos] != '-']
        if aas:
            aa_counts = defaultdict(int)
            for aa in aas:
                aa_counts[aa] += 1
            consensus.append(max(aa_counts, key=aa_counts.get))
        else:
            consensus.append('-')
    consensus_seq = ''.join(consensus)

    # Find motif positions
    motif_positions = {}
    for motif_name, motif_info in GPCR_MOTIFS.items():
        pattern = motif_info['pattern']
        for match in re.finditer(pattern, consensus_seq):
            for pos in range(match.start(), match.end()):
                motif_positions[pos] = {
                    'motif': motif_name,
                    'function': motif_info['function'],
                    'region': motif_info['region']
                }

    # Annotate convergent sites
    for site in convergent_sites:
        pos = site['position'] - 1  # Convert to 0-indexed

        if pos in motif_positions:
            site['in_functional_motif'] = True
            site['motif_info'] = motif_positions[pos]
        else:
            site['in_functional_motif'] = False
            site['motif_info'] = None

        # Estimate if in TM region (rough approximation based on hydrophobicity)
        window_size = 20
        start = max(0, pos - window_size // 2)
        end = min(len(consensus_seq), pos + window_size // 2)
        window = consensus_seq[start:end]

        hydrophobic_count = sum(1 for aa in window if AA_PROPERTIES.get(aa, {}).get('hydrophobic', False))
        site['likely_tm_region'] = hydrophobic_count / len(window) > 0.5 if window else False

    return convergent_sites


def calculate_convergence_score(sites: List[Dict]) -> float:
    """
    Calculate overall convergence score for a set of sites.

    Higher scores indicate more significant convergent evolution.
    """
    if not sites:
        return 0.0

    score = 0.0
    for site in sites:
        # Base score for convergent site
        site_score = 1.0

        # Bonus for functional motif
        if site.get('in_functional_motif'):
            site_score *= 2.0

        # Bonus for TM region
        if site.get('likely_tm_region'):
            site_score *= 1.5

        # Bonus for property-changing substitution
        if site.get('property_change'):
            site_score *= 1.2

        # Weight by frequency
        site_score *= (site['focal_frequency'] + site['other_frequency']) / 2

        score += site_score

    return score


def main():
    """Main execution function."""
    if len(sys.argv) < 4:
        print("Usage: python convergent_evolution.py <alignment> <tree> <output_prefix> [focal_species]")
        sys.exit(1)

    alignment_file = sys.argv[1]
    tree_file = sys.argv[2]
    output_prefix = sys.argv[3]
    focal_species = sys.argv[4].split(',') if len(sys.argv) > 4 else None

    # Validate inputs
    if not os.path.exists(alignment_file):
        print(f"Error: Alignment file not found: {alignment_file}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(tree_file):
        print(f"Error: Tree file not found: {tree_file}", file=sys.stderr)
        sys.exit(1)

    # Load data
    print("Loading alignment and tree...", file=sys.stderr)
    alignment = load_alignment(alignment_file)
    tree = load_tree(tree_file)

    # Identify convergent sites
    print("\nIdentifying convergent amino acid substitutions...", file=sys.stderr)
    convergent_sites = identify_convergent_sites(alignment, tree, focal_species)
    print(f"Found {len(convergent_sites)} convergent sites", file=sys.stderr)

    # Annotate with functional information
    print("Annotating functional regions...", file=sys.stderr)
    convergent_sites = annotate_functional_regions(convergent_sites, alignment)

    # Calculate convergence score
    convergence_score = calculate_convergence_score(convergent_sites)

    # Create output directory
    output_dir = Path(output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Write Outputs ---

    # 1. Convergent sites table
    if convergent_sites:
        # Flatten for CSV output
        flat_sites = []
        for site in convergent_sites:
            flat_site = {
                'position': site['position'],
                'convergent_aa': site['convergent_aa'],
                'ancestral_aa': site['ancestral_aa'],
                'focal_frequency': site['focal_frequency'],
                'other_frequency': site['other_frequency'],
                'in_functional_motif': site['in_functional_motif'],
                'motif_name': site['motif_info']['motif'] if site['motif_info'] else '',
                'likely_tm_region': site['likely_tm_region']
            }
            # Add property changes
            if site.get('property_change'):
                for prop, change in site['property_change'].items():
                    flat_site[f'change_{prop}'] = change
            flat_sites.append(flat_site)

        sites_df = pd.DataFrame(flat_sites)
        sites_file = f"{output_prefix}_convergent_sites.tsv"
        sites_df.to_csv(sites_file, sep='\t', index=False)
        print(f"\nConvergent sites written to: {sites_file}", file=sys.stderr)

    # 2. Summary
    functional_sites = [s for s in convergent_sites if s.get('in_functional_motif')]
    tm_sites = [s for s in convergent_sites if s.get('likely_tm_region')]

    summary = {
        'alignment_file': alignment_file,
        'tree_file': tree_file,
        'n_sequences': len(alignment),
        'alignment_length': alignment.get_alignment_length(),
        'focal_species': focal_species,
        'n_convergent_sites': len(convergent_sites),
        'n_in_functional_motifs': len(functional_sites),
        'n_in_tm_regions': len(tm_sites),
        'convergence_score': convergence_score,
        'top_convergent_positions': [s['position'] for s in sorted(
            convergent_sites, key=lambda x: x['focal_frequency'] * x['other_frequency'],
            reverse=True
        )[:10]]
    }

    summary_file = f"{output_prefix}_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Summary written to: {summary_file}", file=sys.stderr)

    # Print summary
    print("\n=== Convergent Evolution Summary ===", file=sys.stderr)
    print(f"Total convergent sites: {len(convergent_sites)}", file=sys.stderr)
    print(f"Sites in functional motifs: {len(functional_sites)}", file=sys.stderr)
    print(f"Sites in TM regions: {len(tm_sites)}", file=sys.stderr)
    print(f"Convergence score: {convergence_score:.2f}", file=sys.stderr)

    if functional_sites:
        print("\nConvergent sites in functional motifs:", file=sys.stderr)
        for site in functional_sites[:5]:
            print(f"  Position {site['position']}: {site['ancestral_aa'] or '?'} -> {site['convergent_aa']} "
                  f"({site['motif_info']['motif']}, {site['motif_info']['function']})", file=sys.stderr)


if __name__ == "__main__":
    main()
