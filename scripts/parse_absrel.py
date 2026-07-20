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

Writes one per-OG CSV, later joined by a concatenation step (avoids the
column-misalignment risk from raw appends; bead -mqt).

The publish is ATOMIC in the write sense as well as the per-OG-isolation
sense: every CSV goes through write_csv_atomic() -- a per-task-unique temp
in the destination directory, fsync, then os.replace. This docstring
previously said "Supports atomic per-OG CSV writes", which meant only the
per-orthogroup isolation; the write itself was a bare open(path, mode)
truncate-and-refill, and stage 05 globs this file from up to 50 concurrent
array tasks. Keep the two meanings distinct when editing this line.
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


def write_csv_atomic(path, fieldnames, rows) -> None:
    """Publish a CSV atomically: unique temp in the same dir, then os.replace.

    Stage 05 runs as ``#SBATCH --array=0-999%50``, so up to 50 tasks on
    DIFFERENT nodes write per-OG CSVs into one shared directory while other
    tasks glob that same directory to rebuild the cumulative results file.
    A plain ``open(path, "w")`` truncates the published file at open() and
    refills it incrementally, so a concurrent globber can read an empty or
    half-written CSV and concatenate it into the cumulative -- silently,
    with no error anywhere. ``os.replace`` is atomic within a filesystem,
    so readers see either the whole old file or the whole new one.

    The staging name folds in the SLURM array identity, not just the PID:
    PIDs are unique per node only, and this directory is shared storage, so
    two tasks on different nodes can hold the same PID (same reasoning as
    TASK_TAG in scripts/hpc/run_selection_stack.sh). The ".tmp.<tag>" suffix
    also keeps staging files out of the consumers' "*_absrel.csv" globs.

    Byte-identical to the helper in parse_busted.py / parse_meme.py: the
    three stage-05 per-OG writers deliberately share ONE mechanism.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tag = "{}.{}.{}".format(
        os.environ.get("SLURM_JOB_ID", "nojob"),
        os.environ.get("SLURM_ARRAY_TASK_ID", "0"),
        os.getpid(),
    )
    tmp = path.with_name(f"{path.name}.tmp.{tag}")
    try:
        with open(tmp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        # Never leave debris in a directory that concurrent tasks glob.
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise

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
    'cds_source',           # bead -325: native | miniprot | unknown
                            #   (joined from cds_provenance.csv if available)
]


def load_cds_provenance(path):
    """Return dict seq_id -> source (native|miniprot|unknown).

    Reads the CSV produced by recover_cds_from_assemblies.py --manifest.
    Sequences not in the manifest default to 'native' (i.e. they came from
    the original efetch / fetch_reference_cds.py path or the Berghia genome
    annotation, both of which are considered native).
    """
    if not path or not os.path.exists(path):
        return {}
    out = {}
    try:
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                sid = row.get('seq_id', '').strip()
                src = row.get('source', '').strip() or 'unknown'
                if sid:
                    out[sid] = src
    except Exception as e:
        print(f"Warning: could not read CDS provenance manifest {path}: {e}",
              file=sys.stderr)
    return out


def parse_absrel_json(json_file, output_csv, append_mode='atomic',
                       cds_provenance=None):
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

        cds_source = ((cds_provenance or {}).get(branch_id, 'native')
                      if cds_provenance is not None else 'unknown')
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
            'cds_source': cds_source,
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

        cds_source = ((cds_provenance or {}).get(branch_id, 'native')
                      if cds_provenance is not None else 'unknown')
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
            'cds_source': cds_source,
        })

    # Write CSV. Atomic mode: per-OG file, no append (avoids schema-drift bug).
    #
    # Both modes publish through write_csv_atomic(). Legacy append mode reads
    # the already-published rows and republishes old+new in ONE os.replace
    # rather than opening the live file in "a": an in-place append has the
    # same visibility problem as a truncate-and-refill, because a task
    # globbing this directory can read the file between the two writes.
    # Per-OG CSVs are a few hundred rows, so holding them in memory is free.
    rows_to_write = results
    if append_mode == 'append':
        if os.path.exists(output_csv):
            # Validate header matches before appending (bead -mqt)
            with open(output_csv, newline='') as f:
                reader = csv.reader(f)
                existing_header = next(reader, None)
                if existing_header != FIELDNAMES:
                    raise ValueError(
                        f"Cannot append to {output_csv}: column schema differs.\n"
                        f"Existing: {existing_header}\nExpected: {FIELDNAMES}\n"
                        "Hint: regenerate the per-OG CSV with the new schema "
                        "or use --append-mode atomic."
                    )
                prior = [dict(zip(FIELDNAMES, row)) for row in reader if row]
            rows_to_write = prior + results

    write_csv_atomic(output_csv, FIELDNAMES, rows_to_write)

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
    ap.add_argument('--cds-provenance',
                    default=os.environ.get('CDS_PROVENANCE_CSV', ''),
                    help="Optional path to cds_provenance.csv (from "
                         "recover_cds_from_assemblies.py --manifest). When given, "
                         "each output row carries a cds_source column "
                         "(native | miniprot) so dN/dS results can be stratified.")
    args = ap.parse_args()

    if not os.path.exists(args.json_file):
        print(f"Error: Input file {args.json_file} not found", file=sys.stderr)
        sys.exit(1)

    cds_prov = load_cds_provenance(args.cds_provenance) if args.cds_provenance else None
    parse_absrel_json(args.json_file, args.output_csv,
                       append_mode=args.append_mode, cds_provenance=cds_prov)
