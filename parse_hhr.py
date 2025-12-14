#!/usr/bin/env python3
# parse_hhr.py
# Purpose: Parse HHblits/HHsearch HHR output files to extract hit information.
# Inputs: HHR file ($1), E-value threshold ($2), output file ($3)
# Outputs: List of query IDs with hits below E-value threshold
# Author: Jorge L. Perez-Moreno, Ph.D.

"""
HHR (HHblits/HHsearch Result) File Parser

HHR format structure:
- Header lines with query info
- Summary table with columns: No, Hit, Prob, E-value, P-value, Score, SS, Cols, Query HMM, Template HMM
- Detailed alignments

This module provides robust parsing that handles format variations across HHblits versions.
"""

import sys
import re
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass


@dataclass
class HHRHit:
    """Represents a single hit from HHR output."""
    hit_num: int
    hit_id: str
    probability: float
    evalue: float
    pvalue: float
    score: float
    ss_score: float
    cols: int
    query_range: Tuple[int, int]
    template_range: Tuple[int, int]


@dataclass
class HHRResult:
    """Represents parsed results for a single query."""
    query_id: str
    query_length: int
    hits: List[HHRHit]


def parse_hhr_file(filepath: str) -> List[HHRResult]:
    """
    Parse an HHR file and return all query results with their hits.

    Args:
        filepath: Path to the HHR file

    Returns:
        List of HHRResult objects, one per query in the file
    """
    results = []
    current_query = None
    current_length = 0
    current_hits = []
    in_hit_table = False

    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip()

            # Parse query line
            if line.startswith('Query'):
                # Save previous query if exists
                if current_query:
                    results.append(HHRResult(
                        query_id=current_query,
                        query_length=current_length,
                        hits=current_hits
                    ))

                # Extract query ID (second field after "Query")
                parts = line.split()
                current_query = parts[1] if len(parts) > 1 else "unknown"
                current_hits = []
                in_hit_table = False

            # Parse match_columns for query length
            elif line.startswith('Match_columns'):
                try:
                    current_length = int(line.split()[1])
                except (IndexError, ValueError):
                    current_length = 0

            # Detect start of hit table
            elif line.strip().startswith('No Hit'):
                in_hit_table = True
                continue

            # Parse hit lines in the summary table
            elif in_hit_table and current_query:
                # Hit lines start with a number
                match = re.match(r'^\s*(\d+)\s+(\S+)\s+', line)
                if match:
                    hit = parse_hit_line(line)
                    if hit:
                        current_hits.append(hit)
                elif line.strip() == '' or line.startswith('>'):
                    # End of hit table
                    in_hit_table = False

    # Don't forget the last query
    if current_query:
        results.append(HHRResult(
            query_id=current_query,
            query_length=current_length,
            hits=current_hits
        ))

    return results


def parse_hit_line(line: str) -> Optional[HHRHit]:
    """
    Parse a single hit line from the HHR summary table.

    Expected format (may vary slightly between versions):
    No Hit                             Prob E-value P-value  Score    SS  Cols Query HMM  Template HMM
     1 sp|P12345|PROT_HUMAN Desc...   99.9 1.2E-50 3.4E-55  234.5  12.3   150    1-150     10-160 (200)

    Args:
        line: A single line from the hit table

    Returns:
        HHRHit object or None if parsing fails
    """
    # Split by whitespace but preserve structure
    parts = line.split()

    if len(parts) < 10:
        return None

    try:
        hit_num = int(parts[0])
        hit_id = parts[1]

        # Find numeric columns (working backwards from fixed positions)
        # Format: Prob E-value P-value Score SS Cols Query_HMM Template_HMM
        # The last columns are most reliable

        # Try to find E-value by looking for scientific notation
        evalue = None
        prob = None

        for i, part in enumerate(parts[2:], start=2):
            # E-value is typically in scientific notation
            if 'E' in part.upper() or 'e' in part:
                try:
                    val = float(part)
                    if evalue is None:
                        evalue = val
                    break
                except ValueError:
                    continue
            # Probability is typically 0-100
            elif evalue is None:
                try:
                    val = float(part)
                    if 0 <= val <= 100:
                        prob = val
                except ValueError:
                    continue

        if evalue is None:
            # Fallback: try to parse from expected position
            try:
                evalue = float(parts[3])
            except (IndexError, ValueError):
                return None

        if prob is None:
            try:
                prob = float(parts[2])
            except (IndexError, ValueError):
                prob = 0.0

        # Parse remaining fields with defaults
        try:
            pvalue = float(parts[4])
        except (IndexError, ValueError):
            pvalue = evalue  # Use E-value as fallback

        try:
            score = float(parts[5])
        except (IndexError, ValueError):
            score = 0.0

        try:
            ss_score = float(parts[6])
        except (IndexError, ValueError):
            ss_score = 0.0

        try:
            cols = int(parts[7])
        except (IndexError, ValueError):
            cols = 0

        # Parse range fields (e.g., "1-150")
        query_range = (0, 0)
        template_range = (0, 0)

        for part in parts[-4:]:
            if '-' in part and not part.startswith('-'):
                try:
                    start, end = part.split('-')[:2]
                    start = int(start)
                    end = int(end.split('(')[0])  # Remove trailing (length)
                    if query_range == (0, 0):
                        query_range = (start, end)
                    else:
                        template_range = (start, end)
                except ValueError:
                    continue

        return HHRHit(
            hit_num=hit_num,
            hit_id=hit_id,
            probability=prob,
            evalue=evalue,
            pvalue=pvalue,
            score=score,
            ss_score=ss_score,
            cols=cols,
            query_range=query_range,
            template_range=template_range
        )

    except (ValueError, IndexError) as e:
        return None


def filter_hits_by_evalue(results: List[HHRResult], evalue_threshold: float) -> Dict[str, List[HHRHit]]:
    """
    Filter hits by E-value threshold.

    Args:
        results: List of HHRResult objects
        evalue_threshold: Maximum E-value to include

    Returns:
        Dictionary mapping query IDs to lists of passing hits
    """
    filtered = {}
    for result in results:
        passing_hits = [h for h in result.hits if h.evalue < evalue_threshold]
        if passing_hits:
            filtered[result.query_id] = passing_hits
    return filtered


def get_query_ids_with_hits(filepath: str, evalue_threshold: float) -> List[str]:
    """
    Get list of query IDs that have at least one hit below E-value threshold.

    Args:
        filepath: Path to HHR file
        evalue_threshold: Maximum E-value

    Returns:
        List of query IDs with qualifying hits
    """
    results = parse_hhr_file(filepath)
    filtered = filter_hits_by_evalue(results, evalue_threshold)
    return list(filtered.keys())


def main():
    """Command-line interface for HHR parsing."""
    if len(sys.argv) < 3:
        print("Usage: parse_hhr.py <input.hhr> <evalue_threshold> [output.txt]", file=sys.stderr)
        print("  Extracts query IDs with hits below E-value threshold", file=sys.stderr)
        sys.exit(1)

    hhr_file = sys.argv[1]
    evalue_threshold = float(sys.argv[2])
    output_file = sys.argv[3] if len(sys.argv) > 3 else None

    try:
        query_ids = get_query_ids_with_hits(hhr_file, evalue_threshold)

        if output_file:
            with open(output_file, 'w') as f:
                for qid in sorted(set(query_ids)):
                    f.write(f"{qid}\n")
            print(f"Wrote {len(set(query_ids))} query IDs to {output_file}")
        else:
            for qid in sorted(set(query_ids)):
                print(qid)

    except FileNotFoundError:
        print(f"Error: File not found: {hhr_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing {hhr_file}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
