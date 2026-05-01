#!/usr/bin/env python3
# parse_absrel.py
# Purpose: Parse HyPhy aBSREL JSON output to CSV format for downstream analysis.
# Inputs: aBSREL JSON file ($1), output CSV ($2)
# Outputs: CSV with columns:
#   branch_id, omega, omega_max, omega_mean, weight_at_max, n_rate_classes,
#   p_value, corrected_p_value, is_already_corrected, is_selected
# Author: Jorge L. Perez-Moreno, Ph.D.

"""Parse HyPhy aBSREL output JSON into a tabular CSV.

Reports per-branch:
  - omega: legacy weighted-mean omega (kept for backward compatibility)
  - omega_max: maximum omega across rate classes — the chemoreceptor-relevant
    signal that aBSREL is designed to detect (episodic positive selection
    on a fraction of sites). Bead -ea9.
  - omega_mean: weighted mean across rate classes (same as legacy 'omega')
  - weight_at_max: site-fraction at the maximum-omega rate class
  - n_rate_classes: number of mixture components

Supports atomic per-OG CSV writes plus a final concatenation step (avoids the
column-misalignment risk from raw appends; bead -mqt).
"""

import argparse
import json
import csv
import sys
import os
from pathlib import Path

# Allow imports of the shared helpers from sibling module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _rank_candidates_lib import extract_branch_omega  # noqa: E402

FIELDNAMES = [
    'branch_id',
    'omega',                # legacy weighted-mean (= omega_mean)
    'omega_max',            # bead -ea9: episodic-selection signal
    'omega_mean',
    'weight_at_max',
    'n_rate_classes',
    'p_value',
    'corrected_p_value',
    'is_already_corrected',
    'is_selected',
]


def parse_absrel_json(json_file, output_csv, append_mode='atomic'):
    """Parse one aBSREL JSON file and write a per-OG CSV.

    Args:
        json_file: input HyPhy aBSREL JSON path
        output_csv: output CSV path (per-OG, then concatenate at end)
        append_mode: 'atomic' (default; write new file or fail if exists)
                     or 'append' (legacy behavior; risk column drift)
    """
    with open(json_file, 'r') as f:
        data = json.load(f)

    results = []
    branch_attrs = data.get('branch attributes', {})

    # In aBSREL output, 'branch attributes' has nested "PARTITION_ID" keys (e.g. "0")
    # which contain dicts of branch_name -> branch_data. Detect partition layout
    # by checking if the top-level keys are numeric (partition IDs) AND values
    # are dicts of dicts. Branch names are typically alphanumeric with underscores.
    branch_data_iter = []
    if branch_attrs and all(isinstance(v, dict) for v in branch_attrs.values()):
        keys_look_like_partitions = all(
            k.isdigit() for k in branch_attrs.keys()
        )
        values_are_dicts_of_dicts = all(
            v and all(isinstance(vv, dict) for vv in v.values())
            for v in branch_attrs.values()
        )
        if keys_look_like_partitions and values_are_dicts_of_dicts:
            for partition_id, partition in branch_attrs.items():
                for branch_key, branch_data in partition.items():
                    branch_data_iter.append((branch_key, branch_data))
        else:
            for branch_key, branch_data in branch_attrs.items():
                branch_data_iter.append((branch_key, branch_data))

    seen_ids = set()
    for branch_key, branch_data in branch_data_iter:
        if not isinstance(branch_data, dict):
            continue
        branch_id = branch_data.get('original name', branch_key)
        if branch_id in seen_ids:
            continue
        seen_ids.add(branch_id)

        omega_stats = extract_branch_omega(branch_data)
        omega_max = omega_stats['omega_max']
        omega_mean = omega_stats['omega_mean']
        weight_at_max = omega_stats['weight_at_max']
        n_rate_classes = omega_stats['n_rate_classes']

        # Legacy 'omega' column = mean for backward compatibility
        if n_rate_classes > 0:
            omega = omega_mean
        else:
            # Fallback: scalar omega field if no rate distributions
            omega = branch_data.get('omega', branch_data.get('Omega', 1.0))
            try:
                omega = float(omega)
            except (TypeError, ValueError):
                omega = 1.0
            omega_max = omega
            omega_mean = omega

        # P-values
        hyphy_corrected_p = branch_data.get('Corrected P-value')
        raw_p_value = branch_data.get('Uncorrected P-value',
                                      branch_data.get('p-value', 1.0))
        try:
            raw_p_value = float(raw_p_value)
        except (TypeError, ValueError):
            raw_p_value = 1.0

        if hyphy_corrected_p is not None:
            try:
                corrected_p = float(hyphy_corrected_p)
                is_already_corrected = 1
            except (TypeError, ValueError):
                corrected_p = raw_p_value
                is_already_corrected = 0
        else:
            corrected_p = raw_p_value
            is_already_corrected = 0

        # Use omega_max for "is_selected" (chemoreceptor-relevant signal)
        is_selected = 1 if (corrected_p < 0.05 and omega_max > 1.0) else 0

        results.append({
            'branch_id': branch_id,
            'omega': omega,
            'omega_max': omega_max,
            'omega_mean': omega_mean,
            'weight_at_max': weight_at_max,
            'n_rate_classes': n_rate_classes,
            'p_value': raw_p_value,
            'corrected_p_value': corrected_p,
            'is_already_corrected': is_already_corrected,
            'is_selected': is_selected,
        })

    # Also check legacy "tested" section
    tested = data.get('tested', {})
    for branch_id, test_result in tested.items():
        if not isinstance(test_result, dict) or branch_id in seen_ids:
            continue
        omega = test_result.get('omega', test_result.get('Omega', 1.0))
        try:
            omega = float(omega)
        except (TypeError, ValueError):
            omega = 1.0
        raw_p = test_result.get('p', test_result.get('p-value', 1.0))
        try:
            raw_p = float(raw_p)
        except (TypeError, ValueError):
            raw_p = 1.0
        hyphy_p = test_result.get('corrected p')
        if hyphy_p is not None:
            try:
                corrected_p = float(hyphy_p)
                is_already_corrected = 1
            except (TypeError, ValueError):
                corrected_p = raw_p
                is_already_corrected = 0
        else:
            corrected_p = raw_p
            is_already_corrected = 0

        results.append({
            'branch_id': branch_id,
            'omega': omega,
            'omega_max': omega,
            'omega_mean': omega,
            'weight_at_max': float('nan'),
            'n_rate_classes': 0,
            'p_value': raw_p,
            'corrected_p_value': corrected_p,
            'is_already_corrected': is_already_corrected,
            'is_selected': 1 if (corrected_p < 0.05 and omega > 1) else 0,
        })

    # Write CSV. Atomic mode: per-OG file, no append (avoids schema-drift bug).
    if append_mode == 'append':
        file_exists = os.path.exists(output_csv)
        if file_exists:
            # Validate header matches before appending (bead -mqt)
            with open(output_csv) as f:
                existing_header = next(csv.reader(f), None)
            if existing_header != FIELDNAMES:
                raise ValueError(
                    f"Cannot append to {output_csv}: column schema differs.\n"
                    f"Existing: {existing_header}\nExpected: {FIELDNAMES}\n"
                    "Hint: regenerate the per-OG CSV with the new schema "
                    "or use --append-mode atomic."
                )
        mode = 'a' if file_exists else 'w'
    else:
        mode = 'w'

    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, mode, newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        if mode == 'w':
            writer.writeheader()
        writer.writerows(results)

    print(f"Parsed {len(results)} branches from {json_file} -> {output_csv}",
          file=sys.stderr)
    return len(results)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n', 1)[0])
    ap.add_argument('json_file', help='HyPhy aBSREL JSON output')
    ap.add_argument('output_csv', help='Output CSV path')
    ap.add_argument('--append-mode', choices=['atomic', 'append'], default='atomic',
                    help="atomic = per-OG CSV (recommended; concatenate later); "
                         "append = legacy in-place append (validates schema first)")
    args = ap.parse_args()

    if not os.path.exists(args.json_file):
        print(f"Error: Input file {args.json_file} not found", file=sys.stderr)
        sys.exit(1)

    parse_absrel_json(args.json_file, args.output_csv, append_mode=args.append_mode)
