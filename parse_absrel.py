#!/usr/bin/env python3
# parse_absrel.py
# Purpose: Parse HyPhy aBSREL JSON output to CSV format for downstream analysis.
# Inputs: aBSREL JSON file ($1), output CSV ($2)
# Outputs: CSV with columns: branch_id, omega, p_value, corrected_p_value, is_selected
# Author: Jorge L. Perez-Moreno, Ph.D.

import json
import csv
import sys
import os

def parse_absrel_json(json_file, output_csv):
    """
    Parse HyPhy aBSREL JSON output.

    aBSREL tests for episodic diversifying selection on individual branches.
    Output includes per-branch omega (dN/dS) estimates and p-values.
    """
    with open(json_file, 'r') as f:
        data = json.load(f)

    results = []

    # aBSREL stores branch-specific results in "branch attributes"
    branch_attrs = data.get('branch attributes', {})

    # Get the tested branches (usually indexed by node number or name)
    for branch_key, branch_data in branch_attrs.items():
        if isinstance(branch_data, dict):
            # Each branch may have multiple omega rate classes
            # We extract the overall or weighted omega

            branch_id = branch_data.get('original name', branch_key)

            # Get omega (dN/dS) - aBSREL may report multiple rate classes
            # Use the "Rate Distributions" if available
            rate_dist = branch_data.get('Rate Distributions', [])

            if rate_dist:
                # Calculate weighted average omega across rate classes
                total_weight = 0
                weighted_omega = 0
                for rate_class in rate_dist:
                    if isinstance(rate_class, (list, tuple)) and len(rate_class) >= 2:
                        omega, weight = rate_class[0], rate_class[1]
                        weighted_omega += omega * weight
                        total_weight += weight
                omega = weighted_omega / total_weight if total_weight > 0 else 1.0
            else:
                # Fallback to direct omega value if present
                omega = branch_data.get('omega', branch_data.get('Omega', 1.0))

            # Get p-value for positive selection
            # Track whether HyPhy already provided corrected p-values
            hyphy_corrected_p = branch_data.get('Corrected P-value')
            raw_p_value = branch_data.get('Uncorrected P-value',
                          branch_data.get('p-value', 1.0))

            # Use corrected if available, otherwise raw
            if hyphy_corrected_p is not None:
                p_value = hyphy_corrected_p
                is_already_corrected = 1
            else:
                p_value = raw_p_value
                is_already_corrected = 0

            corrected_p = hyphy_corrected_p if hyphy_corrected_p is not None else p_value

            # Determine if branch shows evidence of positive selection
            is_selected = 1 if corrected_p < 0.05 and omega > 1 else 0

            results.append({
                'branch_id': branch_id,
                'omega': omega,
                'p_value': raw_p_value,  # Always store raw p-value
                'corrected_p_value': corrected_p,
                'is_already_corrected': is_already_corrected,  # Flag for downstream
                'is_selected': is_selected
            })

    # Also check for "tested" section in newer HyPhy versions
    tested = data.get('tested', {})
    for branch_id, test_result in tested.items():
        if isinstance(test_result, dict) and branch_id not in [r['branch_id'] for r in results]:
            omega = test_result.get('omega', test_result.get('Omega', 1.0))
            raw_p_value = test_result.get('p', test_result.get('p-value', 1.0))
            hyphy_corrected_p = test_result.get('corrected p')

            if hyphy_corrected_p is not None:
                corrected_p = hyphy_corrected_p
                is_already_corrected = 1
            else:
                corrected_p = raw_p_value
                is_already_corrected = 0

            is_selected = 1 if corrected_p < 0.05 and omega > 1 else 0

            results.append({
                'branch_id': branch_id,
                'omega': omega,
                'p_value': raw_p_value,
                'corrected_p_value': corrected_p,
                'is_already_corrected': is_already_corrected,
                'is_selected': is_selected
            })

    # Write or append to CSV
    file_exists = os.path.exists(output_csv)
    mode = 'a' if file_exists else 'w'

    fieldnames = ['branch_id', 'omega', 'p_value', 'corrected_p_value', 'is_already_corrected', 'is_selected']
    with open(output_csv, mode, newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerows(results)

    print(f"Parsed {len(results)} branches from {json_file}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: parse_absrel.py <input.json> <output.csv>", file=sys.stderr)
        sys.exit(1)

    json_file = sys.argv[1]
    output_csv = sys.argv[2]

    if not os.path.exists(json_file):
        print(f"Error: Input file {json_file} not found", file=sys.stderr)
        sys.exit(1)

    parse_absrel_json(json_file, output_csv)
