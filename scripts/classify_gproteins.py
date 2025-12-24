#!/usr/bin/env python3
# classify_gproteins.py
# Purpose: Identify and classify G-proteins in transcriptomes for co-expression analysis.
# Uses BLAST/DIAMOND against reference G-proteins to classify Gα subunits.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   python classify_gproteins.py --transcriptomes <files> --reference <fasta>
#                                --classes <tsv> --output <csv>

import os
import sys
import argparse
import subprocess
import tempfile
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
from Bio import SeqIO


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Classify G-proteins in transcriptomes'
    )
    parser.add_argument(
        '--transcriptomes', required=True,
        help='Glob pattern for transcriptome protein files (e.g., "*.aa")'
    )
    parser.add_argument(
        '--reference', required=True,
        help='Reference G-protein FASTA file'
    )
    parser.add_argument(
        '--classes', required=True,
        help='TSV file mapping reference IDs to G-protein classes'
    )
    parser.add_argument(
        '--mollusc-reference', default='',
        help='Optional mollusc-specific G-protein FASTA'
    )
    parser.add_argument(
        '--mollusc-classes', default='',
        help='Optional mollusc-specific class mapping TSV'
    )
    parser.add_argument(
        '--evalue', type=float, default=1e-10,
        help='E-value threshold for BLAST hits (default: 1e-10)'
    )
    parser.add_argument(
        '--output', required=True,
        help='Output CSV file'
    )
    return parser.parse_args()


def load_class_mapping(class_file: str) -> Dict[str, Dict]:
    """
    Load G-protein class mapping from TSV file.

    Expected format:
    sequence_id    class    subtype    species    notes

    Returns dict: sequence_id -> {class, subtype, species, notes}
    """
    mapping = {}

    if not os.path.exists(class_file):
        return mapping

    try:
        df = pd.read_csv(class_file, sep='\t')
        for _, row in df.iterrows():
            seq_id = str(row['sequence_id'])
            mapping[seq_id] = {
                'class': row.get('class', ''),
                'subtype': row.get('subtype', ''),
                'species': row.get('species', ''),
                'notes': row.get('notes', '')
            }
    except Exception as e:
        print(f"Warning: Could not load class mapping from {class_file}: {e}",
              file=sys.stderr)

    return mapping


def find_transcriptome_files(pattern: str) -> List[str]:
    """
    Find transcriptome files matching pattern.

    Args:
        pattern: Glob pattern (may include wildcards)

    Returns:
        List of matching file paths
    """
    import glob
    files = glob.glob(pattern)

    # Also try with directory expansion
    if not files:
        base_dir = os.path.dirname(pattern) or '.'
        file_pattern = os.path.basename(pattern)
        if os.path.isdir(base_dir):
            files = glob.glob(os.path.join(base_dir, file_pattern))

    return [f for f in files if os.path.isfile(f)]


def run_blast_search(query_file: str, subject_file: str, evalue: float) -> List[Dict]:
    """
    Run BLAST search (blastp or diamond).

    Returns list of hits with: query_id, subject_id, evalue, pident, bitscore
    """
    hits = []

    # Check for diamond first (faster), then blastp
    diamond_path = subprocess.run(['which', 'diamond'], capture_output=True, text=True)
    use_diamond = diamond_path.returncode == 0

    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp_out:
        tmp_outfile = tmp_out.name

    try:
        if use_diamond:
            # Create diamond database
            with tempfile.NamedTemporaryFile(suffix='.dmnd', delete=False) as tmp_db:
                tmp_dbfile = tmp_db.name

            subprocess.run([
                'diamond', 'makedb',
                '--in', subject_file,
                '--db', tmp_dbfile
            ], capture_output=True)

            # Run diamond blastp
            result = subprocess.run([
                'diamond', 'blastp',
                '--query', query_file,
                '--db', tmp_dbfile,
                '--evalue', str(evalue),
                '--outfmt', '6', 'qseqid', 'sseqid', 'evalue', 'pident', 'bitscore',
                '--out', tmp_outfile,
                '--max-target-seqs', '1'
            ], capture_output=True, text=True)

            os.unlink(tmp_dbfile)
        else:
            # Run blastp
            result = subprocess.run([
                'blastp',
                '-query', query_file,
                '-subject', subject_file,
                '-evalue', str(evalue),
                '-outfmt', '6 qseqid sseqid evalue pident bitscore',
                '-out', tmp_outfile,
                '-max_target_seqs', '1'
            ], capture_output=True, text=True)

        # Parse results
        if os.path.exists(tmp_outfile):
            with open(tmp_outfile) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        hits.append({
                            'query_id': parts[0],
                            'subject_id': parts[1],
                            'evalue': float(parts[2]),
                            'pident': float(parts[3]),
                            'bitscore': float(parts[4])
                        })

    except Exception as e:
        print(f"Warning: BLAST search failed: {e}", file=sys.stderr)
    finally:
        if os.path.exists(tmp_outfile):
            os.unlink(tmp_outfile)

    return hits


def classify_gproteins(transcriptome_files: List[str],
                       reference_file: str,
                       class_mapping: Dict[str, Dict],
                       mollusc_reference: Optional[str],
                       mollusc_mapping: Optional[Dict[str, Dict]],
                       evalue: float) -> pd.DataFrame:
    """
    Classify G-proteins in transcriptomes.

    Args:
        transcriptome_files: List of protein FASTA files
        reference_file: Reference G-protein FASTA
        class_mapping: Reference ID to class mapping
        mollusc_reference: Optional mollusc-specific reference
        mollusc_mapping: Optional mollusc-specific class mapping
        evalue: E-value threshold

    Returns:
        DataFrame with classified G-proteins
    """
    all_results = []

    for trans_file in transcriptome_files:
        species = Path(trans_file).stem
        print(f"Processing {species}...", file=sys.stderr)

        # Search against mollusc reference first (if available)
        if mollusc_reference and os.path.exists(mollusc_reference) and mollusc_mapping:
            hits = run_blast_search(trans_file, mollusc_reference, evalue)
            for hit in hits:
                subject = hit['subject_id']
                if subject in mollusc_mapping:
                    info = mollusc_mapping[subject]
                    all_results.append({
                        'gene_id': hit['query_id'],
                        'species': species,
                        'gprotein_class': info['class'],
                        'subtype': info['subtype'],
                        'evalue': hit['evalue'],
                        'pident': hit['pident'],
                        'reference_hit': subject,
                        'classification_source': 'mollusc_reference'
                    })

        # Also search against general reference
        hits = run_blast_search(trans_file, reference_file, evalue)
        for hit in hits:
            subject = hit['subject_id']
            # Skip if already classified from mollusc reference
            if any(r['gene_id'] == hit['query_id'] for r in all_results):
                continue

            if subject in class_mapping:
                info = class_mapping[subject]
                all_results.append({
                    'gene_id': hit['query_id'],
                    'species': species,
                    'gprotein_class': info['class'],
                    'subtype': info['subtype'],
                    'evalue': hit['evalue'],
                    'pident': hit['pident'],
                    'reference_hit': subject,
                    'classification_source': 'general_reference'
                })

    return pd.DataFrame(all_results)


def main():
    """Main execution function."""
    args = parse_arguments()

    # Check reference exists
    if not os.path.exists(args.reference):
        print(f"Error: Reference file not found: {args.reference}", file=sys.stderr)
        sys.exit(1)

    # Load class mappings
    class_mapping = load_class_mapping(args.classes)
    if not class_mapping:
        print(f"Error: Could not load class mapping from {args.classes}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(class_mapping)} reference G-protein classes", file=sys.stderr)

    # Load mollusc-specific mapping if available
    mollusc_mapping = None
    if args.mollusc_classes and os.path.exists(args.mollusc_classes):
        mollusc_mapping = load_class_mapping(args.mollusc_classes)
        print(f"Loaded {len(mollusc_mapping)} mollusc-specific G-protein classes",
              file=sys.stderr)

    # Find transcriptome files
    trans_files = find_transcriptome_files(args.transcriptomes)
    if not trans_files:
        print(f"Warning: No transcriptome files found matching {args.transcriptomes}",
              file=sys.stderr)
        # Create empty output
        pd.DataFrame(columns=[
            'gene_id', 'species', 'gprotein_class', 'subtype', 'evalue',
            'pident', 'reference_hit', 'classification_source'
        ]).to_csv(args.output, index=False)
        return

    print(f"Found {len(trans_files)} transcriptome files", file=sys.stderr)

    # Classify G-proteins
    result_df = classify_gproteins(
        trans_files,
        args.reference,
        class_mapping,
        args.mollusc_reference if args.mollusc_reference else None,
        mollusc_mapping,
        args.evalue
    )

    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write output
    result_df.to_csv(args.output, index=False)

    # Print summary
    print(f"\n=== G-protein Classification Summary ===", file=sys.stderr)
    print(f"Total G-proteins identified: {len(result_df)}", file=sys.stderr)
    if not result_df.empty:
        print(f"By class:", file=sys.stderr)
        for cls, count in result_df['gprotein_class'].value_counts().items():
            print(f"  {cls}: {count}", file=sys.stderr)
    print(f"Output written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
