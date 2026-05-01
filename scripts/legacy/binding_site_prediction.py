#!/usr/bin/env python3
# binding_site_prediction.py
# Purpose: Predict ligand binding sites in GPCR candidates using structural information.
# Maps AlphaFold structures to GPCRdb templates and identifies putative binding residues.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   python binding_site_prediction.py <structure_dir> <sequence_file> <output_prefix>
#
# Arguments:
#   structure_dir  - Directory containing AlphaFold PDB files
#   sequence_file  - FASTA file with candidate sequences
#   output_prefix  - Prefix for output files

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.PDB import PDBParser, NeighborSearch, Selection

# GPCRdb generic residue numbering for Class A GPCRs
# Ballesteros-Weinstein numbering: X.YY where X=TM helix, YY=position relative to conserved residue
GPCR_BINDING_POCKET = {
    # Core orthosteric binding pocket residues (Ballesteros-Weinstein)
    '3.28': {'helix': 3, 'function': 'Orthosteric pocket floor'},
    '3.32': {'helix': 3, 'function': 'Key binding residue'},
    '3.33': {'helix': 3, 'function': 'Ligand recognition'},
    '3.36': {'helix': 3, 'function': 'Toggle switch'},
    '3.37': {'helix': 3, 'function': 'Binding pocket'},
    '5.42': {'helix': 5, 'function': 'Orthosteric pocket'},
    '5.43': {'helix': 5, 'function': 'Binding pocket'},
    '5.46': {'helix': 5, 'function': 'Orthosteric interaction'},
    '6.48': {'helix': 6, 'function': 'Rotamer toggle'},
    '6.51': {'helix': 6, 'function': 'Orthosteric pocket'},
    '6.52': {'helix': 6, 'function': 'Binding pocket'},
    '6.55': {'helix': 6, 'function': 'Ligand recognition'},
    '7.35': {'helix': 7, 'function': 'Binding pocket'},
    '7.39': {'helix': 7, 'function': 'Orthosteric pocket'},
    '7.43': {'helix': 7, 'function': 'Binding pocket'},
}

# Conserved positions for TM helix identification
CONSERVED_MOTIFS = {
    'TM1': {'N': 1.50},  # Most conserved in TM1
    'TM2': {'D': 2.50},  # Conserved Asp in TM2
    'TM3': {'R': 3.50},  # DRY motif - Arg
    'TM4': {'W': 4.50},  # Conserved Trp
    'TM5': {'P': 5.50},  # Conserved Pro
    'TM6': {'P': 6.50},  # CWxP motif
    'TM7': {'P': 7.50},  # NPxxY motif - Pro
}

# Reference binding pocket positions from Class A GPCRs (absolute positions in canonical sequence)
# These are approximate and should be refined with structural alignment
REFERENCE_POCKET_RESIDUES = {
    'rhodopsin': [113, 117, 118, 121, 122, 207, 208, 211, 265, 268, 269, 272, 293, 297, 301],
    'beta2_adrenergic': [113, 117, 118, 121, 122, 203, 204, 207, 286, 289, 290, 293, 312, 316, 320],
}

# Kyte-Doolittle hydrophobicity scale (negative = hydrophilic, positive = hydrophobic)
HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Three-letter to one-letter amino acid code mapping
AA_3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# Get pocket analysis mode from environment
POCKET_ANALYSIS_MODE = os.getenv('POCKET_ANALYSIS_MODE', 'aquatic')
MIN_POCKET_VOLUME = float(os.getenv('MIN_POCKET_VOLUME', 100))
MAX_POCKET_VOLUME = float(os.getenv('MAX_POCKET_VOLUME', 1000))


def parse_pdb_structure(pdb_file: str) -> Optional[Dict]:
    """
    Parse PDB structure and extract relevant information.

    Args:
        pdb_file: Path to PDB file

    Returns:
        Dictionary with structure information
    """
    parser = PDBParser(QUIET=True)

    try:
        structure = parser.get_structure('protein', pdb_file)
        model = structure[0]

        # Get chain (assume single chain for AlphaFold)
        chain = list(model.get_chains())[0]

        # Extract CA coordinates and residues
        residues = []
        for res in chain.get_residues():
            if res.get_resname() in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE',
                                     'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER',
                                     'THR', 'VAL', 'TRP', 'TYR']:
                ca = res['CA'] if 'CA' in res else None
                if ca:
                    residues.append({
                        'position': res.get_id()[1],
                        'resname': res.get_resname(),
                        'ca_coord': ca.get_coord().tolist(),
                        'chain': chain.id
                    })

        # Get pLDDT scores from B-factor (AlphaFold stores pLDDT there)
        plddt_scores = []
        for res in chain.get_residues():
            try:
                ca = res['CA']
                plddt_scores.append(ca.get_bfactor())
            except KeyError:
                # Residue doesn't have a CA atom (e.g., glycine or incomplete residue)
                continue

        return {
            'file': pdb_file,
            'n_residues': len(residues),
            'residues': residues,
            'mean_plddt': np.mean(plddt_scores) if plddt_scores else 0,
            'structure': structure
        }

    except Exception as e:
        print(f"Warning: Could not parse {pdb_file}: {e}", file=sys.stderr)
        return None


def identify_tm_helices(residues: List[Dict], sequence: str) -> Dict[int, int]:
    """
    Identify TM helices based on sequence patterns.

    Simple hydrophobicity-based prediction.

    Args:
        residues: List of residue information
        sequence: Amino acid sequence

    Returns:
        Dictionary mapping residue positions to TM helix numbers
    """
    # Hydrophobicity scale (Kyte-Doolittle)
    hydrophobicity = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
        'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
        'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }

    # Calculate hydrophobicity profile
    window_size = 19  # Typical TM helix length
    profile = []

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        score = np.mean([hydrophobicity.get(aa, 0) for aa in window])
        profile.append(score)

    # Find TM helices (regions with high hydrophobicity)
    tm_assignments = {}
    threshold = 1.0
    current_tm = 0
    in_tm = False

    for i, score in enumerate(profile):
        if score > threshold and not in_tm:
            in_tm = True
            current_tm += 1
        elif score <= threshold and in_tm:
            in_tm = False

        if in_tm and current_tm <= 7:  # GPCRs have 7 TM helices
            for pos in range(i, min(i + window_size, len(sequence))):
                if pos + 1 not in tm_assignments:  # 1-indexed
                    tm_assignments[pos + 1] = current_tm

    return tm_assignments


def predict_binding_pocket(structure_info: Dict,
                          sequence: str,
                          tm_assignments: Dict[int, int]) -> List[Dict]:
    """
    Predict ligand binding pocket residues.

    Uses TM helix assignments and known binding pocket positions.

    Args:
        structure_info: Parsed PDB structure
        sequence: Amino acid sequence
        tm_assignments: TM helix assignments

    Returns:
        List of predicted binding pocket residues
    """
    binding_residues = []
    residues = structure_info['residues']

    # Get centroid of structure
    ca_coords = np.array([r['ca_coord'] for r in residues])
    centroid = np.mean(ca_coords, axis=0)

    # Find residues facing the central cavity
    for res in residues:
        pos = res['position']
        tm_helix = tm_assignments.get(pos, 0)

        if tm_helix == 0:
            continue  # Skip non-TM residues

        # Check if residue faces central cavity
        coord = np.array(res['ca_coord'])
        dist_to_center = np.linalg.norm(coord - centroid)

        # Estimate side chain orientation (simplified)
        # Residues closer to center are more likely to face the binding pocket
        tm_distances = [np.linalg.norm(np.array(r['ca_coord']) - centroid)
                       for r in residues if tm_assignments.get(r['position'], 0) > 0]
        # Need at least 2 residues for meaningful percentile calculation
        if len(tm_distances) >= 2:
            is_pocket_facing = dist_to_center < np.percentile(tm_distances, 50)
        else:
            is_pocket_facing = True  # Default to pocket-facing if insufficient data

        # Check if in binding pocket region based on TM position
        pocket_probability = 0.0

        # TM3 and TM6 are most important for binding
        if tm_helix == 3:
            pocket_probability += 0.3
        elif tm_helix == 6:
            pocket_probability += 0.3
        elif tm_helix in [5, 7]:
            pocket_probability += 0.2

        if is_pocket_facing:
            pocket_probability += 0.3

        # Aromatic residues often important for binding
        if res['resname'] in ['PHE', 'TYR', 'TRP', 'HIS']:
            pocket_probability += 0.2

        if pocket_probability > 0.4:
            binding_residues.append({
                'position': pos,
                'resname': res['resname'],
                'tm_helix': tm_helix,
                'distance_to_center': float(dist_to_center),
                'pocket_facing': is_pocket_facing,
                'pocket_probability': pocket_probability,
                'ca_coord': res['ca_coord']
            })

    # Sort by probability
    binding_residues.sort(key=lambda x: x['pocket_probability'], reverse=True)

    return binding_residues


def calculate_pocket_conservation(binding_residues: List[Dict],
                                 alignment_file: Optional[str] = None) -> List[Dict]:
    """
    Calculate conservation scores for predicted binding residues.

    Args:
        binding_residues: List of binding residue predictions
        alignment_file: Optional MSA file for conservation calculation

    Returns:
        Updated binding residues with conservation scores
    """
    # If no alignment, return with placeholder conservation
    if not alignment_file or not os.path.exists(alignment_file):
        for res in binding_residues:
            res['conservation_score'] = None
        return binding_residues

    # TODO: Implement proper conservation calculation from MSA
    # For now, return placeholder
    return binding_residues


def find_cavity_volume(structure_info: Dict, binding_residues: List[Dict]) -> float:
    """
    Estimate binding cavity volume.

    Args:
        structure_info: Parsed structure
        binding_residues: List of binding pocket residues

    Returns:
        Estimated cavity volume (Å³)
    """
    if len(binding_residues) < 4:
        return 0.0

    coords = np.array([r['ca_coord'] for r in binding_residues])

    # Simple convex hull volume approximation
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(coords)
        return hull.volume
    except Exception:
        # Fallback: bounding box (ConvexHull can fail with coplanar points)
        min_coords = coords.min(axis=0)
        max_coords = coords.max(axis=0)
        return np.prod(max_coords - min_coords)


def calculate_pocket_hydrophilicity(pocket_residues: List[Dict]) -> Dict:
    """
    Calculate hydrophilicity metrics for binding pocket residues.

    For aquatic chemoreceptors: expect negative mean (hydrophilic)
    For terrestrial olfactory: expect positive mean (hydrophobic)

    Args:
        pocket_residues: List of pocket residue dicts with 'resname'

    Returns:
        Dict with hydrophilicity metrics
    """
    if not pocket_residues:
        return {
            'mean_hydrophobicity': 0.0,
            'hydrophilicity': 0.0,
            'n_charged': 0,
            'n_aromatic': 0,
            'charged_fraction': 0.0,
            'aromatic_fraction': 0.0
        }

    hydrophobicity_scores = []
    n_charged = 0
    n_aromatic = 0

    charged_residues = {'ARG', 'LYS', 'ASP', 'GLU', 'HIS'}
    aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}

    for res in pocket_residues:
        resname = res.get('resname', '')

        # Convert to 1-letter code if needed
        if len(resname) == 3:
            aa = AA_3TO1.get(resname.upper(), 'X')
        else:
            aa = resname[0] if resname else 'X'

        score = HYDROPHOBICITY.get(aa, 0.0)
        hydrophobicity_scores.append(score)

        if resname.upper() in charged_residues:
            n_charged += 1
        if resname.upper() in aromatic_residues:
            n_aromatic += 1

    mean_hydrophobicity = np.mean(hydrophobicity_scores) if hydrophobicity_scores else 0.0
    n_residues = len(pocket_residues)

    return {
        'mean_hydrophobicity': mean_hydrophobicity,
        'hydrophilicity': -mean_hydrophobicity,  # Invert for hydrophilicity
        'n_charged': n_charged,
        'n_aromatic': n_aromatic,
        'charged_fraction': n_charged / n_residues if n_residues > 0 else 0.0,
        'aromatic_fraction': n_aromatic / n_residues if n_residues > 0 else 0.0
    }


def score_pocket_for_mode(hydrophilicity: float, mode: str) -> float:
    """
    Score pocket based on expected properties for ligand type.

    Aquatic: water-soluble ligands (amino acids, nucleotides) -> hydrophilic pocket preferred
    Terrestrial: volatile ligands -> hydrophobic pocket preferred
    Both: return raw hydrophilicity value

    Args:
        hydrophilicity: Mean hydrophilicity of pocket
        mode: 'aquatic', 'terrestrial', or 'both'

    Returns:
        Score (higher = better match for expected ligand type)
    """
    if mode == 'aquatic':
        # Higher hydrophilicity = better for aquatic
        return max(0, hydrophilicity)
    elif mode == 'terrestrial':
        # Lower hydrophilicity (higher hydrophobicity) = better
        return max(0, -hydrophilicity)
    else:  # 'both'
        # Return raw value for interpretation
        return hydrophilicity


def check_pocket_volume_suitable(volume: float, min_vol: float, max_vol: float) -> bool:
    """Check if pocket volume is suitable for small molecule binding."""
    return min_vol <= volume <= max_vol


def main():
    """Main execution function."""
    if len(sys.argv) < 4:
        print("Usage: python binding_site_prediction.py <structure_dir> <sequence_file> <output_prefix>")
        sys.exit(1)

    structure_dir = Path(sys.argv[1])
    sequence_file = sys.argv[2]
    output_prefix = sys.argv[3]

    # Validate inputs
    if not structure_dir.exists():
        print(f"Error: Structure directory not found: {structure_dir}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(sequence_file):
        print(f"Error: Sequence file not found: {sequence_file}", file=sys.stderr)
        sys.exit(1)

    # Load sequences
    print("Loading sequences...", file=sys.stderr)
    sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(sequence_file, "fasta")}
    print(f"Loaded {len(sequences)} sequences", file=sys.stderr)

    # Find PDB files
    pdb_files = list(structure_dir.glob("*.pdb")) + list(structure_dir.glob("*/*.pdb"))
    print(f"Found {len(pdb_files)} PDB files", file=sys.stderr)

    # Process each structure
    results = []
    all_binding_residues = []

    for pdb_file in pdb_files:
        pdb_name = pdb_file.stem
        print(f"\nProcessing: {pdb_name}", file=sys.stderr)

        # Parse structure
        structure_info = parse_pdb_structure(str(pdb_file))
        if not structure_info:
            continue

        # Find matching sequence
        seq_id = None
        sequence = None
        for sid, seq in sequences.items():
            if pdb_name in sid or sid in pdb_name:
                seq_id = sid
                sequence = seq
                break

        if not sequence:
            # Use sequence from structure
            sequence = ''.join([r['resname'][0] if len(r['resname']) == 3 else r['resname']
                               for r in structure_info['residues']])

        # Identify TM helices
        tm_assignments = identify_tm_helices(structure_info['residues'], sequence)

        # Predict binding pocket
        binding_residues = predict_binding_pocket(structure_info, sequence, tm_assignments)

        # Estimate cavity volume
        cavity_volume = find_cavity_volume(structure_info, binding_residues[:20])

        # Calculate pocket hydrophilicity (Phase 4 enhancement)
        hydrophilicity_metrics = calculate_pocket_hydrophilicity(binding_residues[:20])
        pocket_score = score_pocket_for_mode(
            hydrophilicity_metrics['hydrophilicity'],
            POCKET_ANALYSIS_MODE
        )
        volume_suitable = check_pocket_volume_suitable(
            cavity_volume, MIN_POCKET_VOLUME, MAX_POCKET_VOLUME
        )

        # Store results
        result = {
            'structure': pdb_name,
            'sequence_id': seq_id or pdb_name,
            'n_residues': structure_info['n_residues'],
            'mean_plddt': structure_info['mean_plddt'],
            'n_tm_helices': len(set(tm_assignments.values())),
            'n_binding_residues': len(binding_residues),
            'top_binding_residues': [r['position'] for r in binding_residues[:15]],
            'cavity_volume': cavity_volume,
            # Phase 4: Hydrophilicity metrics
            'pocket_hydrophilicity': hydrophilicity_metrics['hydrophilicity'],
            'pocket_hydrophobicity': hydrophilicity_metrics['mean_hydrophobicity'],
            'n_charged_residues': hydrophilicity_metrics['n_charged'],
            'n_aromatic_residues': hydrophilicity_metrics['n_aromatic'],
            'charged_fraction': hydrophilicity_metrics['charged_fraction'],
            'aromatic_fraction': hydrophilicity_metrics['aromatic_fraction'],
            'pocket_score': pocket_score,
            'volume_suitable': volume_suitable,
            'analysis_mode': POCKET_ANALYSIS_MODE
        }
        results.append(result)

        # Store detailed binding residues
        for res in binding_residues[:20]:  # Top 20 residues
            all_binding_residues.append({
                'structure': pdb_name,
                **res
            })

    # Create output directory
    output_dir = Path(output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Write Outputs ---

    # 1. Structure summary
    summary_df = pd.DataFrame(results)
    summary_file = f"{output_prefix}_structure_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"\nStructure summary written to: {summary_file}", file=sys.stderr)

    # 2. Binding residues
    if all_binding_residues:
        binding_df = pd.DataFrame(all_binding_residues)
        # Remove coordinate column for TSV output
        binding_df_out = binding_df.drop(columns=['ca_coord'], errors='ignore')
        binding_file = f"{output_prefix}_binding_residues.tsv"
        binding_df_out.to_csv(binding_file, sep='\t', index=False)
        print(f"Binding residues written to: {binding_file}", file=sys.stderr)

    # 3. Summary JSON
    summary = {
        'n_structures_analyzed': len(results),
        'mean_cavity_volume': np.mean([r['cavity_volume'] for r in results]) if results else 0,
        'mean_binding_residues': np.mean([r['n_binding_residues'] for r in results]) if results else 0,
        'high_confidence_structures': len([r for r in results if r['mean_plddt'] > 70]),
        'complete_tm_structures': len([r for r in results if r['n_tm_helices'] >= 6]),
        # Phase 4: Hydrophilicity summary
        'analysis_mode': POCKET_ANALYSIS_MODE,
        'mean_pocket_hydrophilicity': np.mean([r['pocket_hydrophilicity'] for r in results]) if results else 0,
        'mean_charged_fraction': np.mean([r['charged_fraction'] for r in results]) if results else 0,
        'mean_aromatic_fraction': np.mean([r['aromatic_fraction'] for r in results]) if results else 0,
        'suitable_volume_count': len([r for r in results if r.get('volume_suitable', False)])
    }

    summary_json_file = f"{output_prefix}_summary.json"
    with open(summary_json_file, 'w') as f:
        json.dump(summary, f, indent=2)

    # Print summary
    print("\n=== Binding Site Prediction Summary ===", file=sys.stderr)
    print(f"Structures analyzed: {summary['n_structures_analyzed']}", file=sys.stderr)
    print(f"High confidence (pLDDT > 70): {summary['high_confidence_structures']}", file=sys.stderr)
    print(f"Complete TM (≥6 helices): {summary['complete_tm_structures']}", file=sys.stderr)
    print(f"Mean binding pocket residues: {summary['mean_binding_residues']:.1f}", file=sys.stderr)
    print(f"Mean cavity volume: {summary['mean_cavity_volume']:.1f} Å³", file=sys.stderr)
    print(f"\n=== Pocket Hydrophilicity Analysis (Mode: {POCKET_ANALYSIS_MODE}) ===", file=sys.stderr)
    print(f"Mean pocket hydrophilicity: {summary['mean_pocket_hydrophilicity']:.2f}", file=sys.stderr)
    print(f"Mean charged residue fraction: {summary['mean_charged_fraction']:.2f}", file=sys.stderr)
    print(f"Mean aromatic residue fraction: {summary['mean_aromatic_fraction']:.2f}", file=sys.stderr)
    print(f"Structures with suitable pocket volume: {summary['suitable_volume_count']}", file=sys.stderr)

    if results:
        print("\nTop structures by binding pocket definition:", file=sys.stderr)
        for r in sorted(results, key=lambda x: x['n_binding_residues'], reverse=True)[:5]:
            print(f"  {r['structure']}: {r['n_binding_residues']} residues, "
                  f"pLDDT={r['mean_plddt']:.1f}, volume={r['cavity_volume']:.1f} Å³",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
