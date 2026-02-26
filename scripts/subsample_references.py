#!/usr/bin/env python3
"""
subsample_references.py
Purpose: Subsample reference GPCR sequences for phylogenetic tree construction.
         Two-stage approach: CD-HIT clustering + taxonomy-weighted proportional allocation.
         Operates directly on raw .faa files from nath_et_al/ directory structure.
Inputs:  --ref-dir (nath_et_al/ path), --target-size, --cluster-identity, etc.
Outputs: Subsampled FASTA + JSON report
Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def build_sequence_index(ref_dir):
    """Walk nath_et_al/{lse,one_to_one_ortholog}/{group}/*.faa to build sequence index.

    Returns:
        dict: header -> {taxid, species, category, group, seq_length, file_path}
        list: ordered list of (header, file_path) for deterministic iteration
    """
    index = {}
    ordered = []
    ref_path = Path(ref_dir)

    for category in ["lse", "one_to_one_ortholog"]:
        cat_dir = ref_path / category
        if not cat_dir.is_dir():
            continue
        for group_dir in sorted(cat_dir.iterdir()):
            if not group_dir.is_dir():
                continue
            group = group_dir.name
            for faa_file in sorted(group_dir.glob("*.faa")):
                # Parse taxid and species from filename: {taxid}_{Genus}_{species}.faa
                stem = faa_file.stem
                parts = stem.split("_", 2)
                taxid = parts[0] if parts else stem
                species = "_".join(parts[1:]) if len(parts) > 1 else stem

                for record in SeqIO.parse(str(faa_file), "fasta"):
                    header = record.id
                    if header in index:
                        # Same sequence in both lse and one_to_one: prefer lse
                        if category == "lse":
                            index[header]["category"] = "lse"
                        continue
                    index[header] = {
                        "taxid": taxid,
                        "species": species,
                        "category": category,
                        "group": group,
                        "seq_length": len(record.seq),
                        "file_path": str(faa_file),
                    }
                    ordered.append((header, str(faa_file)))

    return index, ordered


def concatenate_fasta(ref_dir, output_path):
    """Concatenate all .faa files into a single FASTA, deduplicating by header."""
    seen = set()
    ref_path = Path(ref_dir)
    with open(output_path, "w") as out:
        for category in ["lse", "one_to_one_ortholog"]:
            cat_dir = ref_path / category
            if not cat_dir.is_dir():
                continue
            for group_dir in sorted(cat_dir.iterdir()):
                if not group_dir.is_dir():
                    continue
                for faa_file in sorted(group_dir.glob("*.faa")):
                    for record in SeqIO.parse(str(faa_file), "fasta"):
                        if record.id not in seen:
                            seen.add(record.id)
                            SeqIO.write([record], out, "fasta")
    return len(seen)


def run_cdhit(input_fasta, output_fasta, identity, memory, threads):
    """Run CD-HIT and return set of representative sequence IDs."""
    # Auto-select word size based on identity threshold
    if identity >= 0.7:
        word_size = 5
    elif identity >= 0.6:
        word_size = 4
    elif identity >= 0.5:
        word_size = 3
    else:
        word_size = 2

    cmd = [
        "cd-hit",
        "-i", input_fasta,
        "-o", output_fasta,
        "-c", str(identity),
        "-n", str(word_size),
        "-M", str(memory),
        "-T", str(threads),
        "-d", "0",  # Use full header in .clstr file
    ]

    print(f"Running CD-HIT: identity={identity}, word_size={word_size}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"CD-HIT stderr: {result.stderr}", file=sys.stderr)
        raise RuntimeError(f"CD-HIT failed with return code {result.returncode}")

    # Parse .clstr file to get representative IDs
    representatives = set()
    clstr_file = output_fasta + ".clstr"
    with open(clstr_file) as f:
        for line in f:
            if line.strip().endswith("*"):
                # Representative sequence line: "0\t500aa, >header... *"
                start = line.index(">") + 1
                end = line.index("...", start)
                rep_id = line[start:end]
                representatives.add(rep_id)

    return representatives


def allocate_quotas(bin_counts, bin_weights, target):
    """Allocate target slots across bins using weighted proportional allocation.

    Args:
        bin_counts: dict of bin_key -> count of available sequences
        bin_weights: dict of bin_key -> weight
        target: total number of sequences to select

    Returns:
        dict of bin_key -> quota (number to select)
    """
    if not bin_counts:
        return {}

    total_available = sum(bin_counts.values())
    if total_available <= target:
        return dict(bin_counts)

    # Compute weighted totals
    total_weighted = sum(
        bin_counts[k] * bin_weights.get(k, 1.0) for k in bin_counts
    )

    # Initial allocation: proportional to weighted count
    quotas = {}
    for k, count in bin_counts.items():
        weight = bin_weights.get(k, 1.0)
        raw_quota = (count * weight / total_weighted) * target
        quotas[k] = max(1, int(raw_quota))  # Minimum 1 per non-empty bin

    # Adjust to hit target exactly
    current_total = sum(quotas.values())

    if current_total > target:
        # Trim from lowest-weight bins first
        bins_by_weight = sorted(quotas.keys(), key=lambda k: bin_weights.get(k, 1.0))
        for k in bins_by_weight:
            if current_total <= target:
                break
            reducible = quotas[k] - 1  # Keep at least 1
            reduce_by = min(reducible, current_total - target)
            if reduce_by > 0:
                quotas[k] -= reduce_by
                current_total -= reduce_by

    elif current_total < target:
        # Fill from highest-weight bins first
        bins_by_weight = sorted(
            quotas.keys(), key=lambda k: bin_weights.get(k, 1.0), reverse=True
        )
        for k in bins_by_weight:
            if current_total >= target:
                break
            fillable = bin_counts[k] - quotas[k]
            fill_by = min(fillable, target - current_total)
            if fill_by > 0:
                quotas[k] += fill_by
                current_total += fill_by

    return quotas


def subsample(index, representatives, target_size, taxonomy_weights, lse_weight):
    """Select sequences using taxonomy-weighted proportional allocation.

    Args:
        index: sequence index from build_sequence_index
        representatives: set of CD-HIT representative IDs
        target_size: desired number of output sequences
        taxonomy_weights: dict of group -> weight
        lse_weight: multiplier for LSE vs one_to_one_ortholog

    Returns:
        set of selected sequence headers
    """
    # Filter index to representatives only
    rep_index = {h: info for h, info in index.items() if h in representatives}

    if len(rep_index) <= target_size:
        return set(rep_index.keys())

    # Group representatives into bins by (group, category)
    bins = defaultdict(list)
    for header, info in rep_index.items():
        bin_key = (info["group"], info["category"])
        bins[bin_key].append((header, info))

    # Compute bin counts and weights
    bin_counts = {k: len(v) for k, v in bins.items()}
    bin_weights = {}
    for (group, category), _ in bins.items():
        tax_w = taxonomy_weights.get(group, 1.0)
        cat_w = lse_weight if category == "lse" else 1.0
        bin_weights[(group, category)] = tax_w * cat_w

    # Allocate quotas
    quotas = allocate_quotas(bin_counts, bin_weights, target_size)

    # Within each bin, sort by sequence length descending (deterministic), take top N
    selected = set()
    for bin_key, entries in bins.items():
        quota = quotas.get(bin_key, 0)
        # Sort by length descending, then by header for determinism
        entries.sort(key=lambda x: (-x[1]["seq_length"], x[0]))
        for header, _ in entries[:quota]:
            selected.add(header)

    return selected


def extract_selected(ref_dir, selected_headers, output_path):
    """Extract selected sequences from .faa files and write to output FASTA."""
    remaining = set(selected_headers)
    ref_path = Path(ref_dir)
    written = 0

    with open(output_path, "w") as out:
        for category in ["lse", "one_to_one_ortholog"]:
            cat_dir = ref_path / category
            if not cat_dir.is_dir():
                continue
            for group_dir in sorted(cat_dir.iterdir()):
                if not group_dir.is_dir():
                    continue
                for faa_file in sorted(group_dir.glob("*.faa")):
                    for record in SeqIO.parse(str(faa_file), "fasta"):
                        if record.id in remaining:
                            SeqIO.write([record], out, "fasta")
                            remaining.discard(record.id)
                            written += 1
                    if not remaining:
                        return written

    return written


def write_report(report_path, input_count, post_cluster_count, final_count,
                 group_breakdown, taxonomy_weights, lse_weight, cluster_identity):
    """Write JSON report with subsampling statistics."""
    report = {
        "input_sequences": input_count,
        "post_clustering_representatives": post_cluster_count,
        "final_selected": final_count,
        "cluster_identity": cluster_identity,
        "lse_weight": lse_weight,
        "taxonomy_weights": taxonomy_weights,
        "per_group_breakdown": group_breakdown,
    }
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)


def parse_taxonomy_weights(weights_str):
    """Parse 'gastropoda:3.0,bivalvia:1.5,...' into dict."""
    weights = {}
    if not weights_str:
        return weights
    for pair in weights_str.split(","):
        pair = pair.strip()
        if ":" in pair:
            group, weight = pair.rsplit(":", 1)
            weights[group.strip()] = float(weight.strip())
    return weights


def main():
    parser = argparse.ArgumentParser(
        description="Subsample reference GPCR sequences for phylogenetic trees"
    )
    parser.add_argument("--ref-dir", required=True,
                        help="Path to nath_et_al/ directory")
    parser.add_argument("--target-size", type=int, default=2000,
                        help="Desired number of reference sequences (default: 2000)")
    parser.add_argument("--cluster-identity", type=float, default=0.7,
                        help="CD-HIT clustering identity threshold (default: 0.7)")
    parser.add_argument("--taxonomy-weights", type=str, default="",
                        help="Group weights string, e.g. 'gastropoda:3.0,bivalvia:1.5'")
    parser.add_argument("--lse-weight", type=float, default=1.5,
                        help="Multiplier for LSE vs one-to-one orthologs (default: 1.5)")
    parser.add_argument("--cdhit-path", default="cd-hit",
                        help="Path to CD-HIT binary")
    parser.add_argument("--cdhit-memory", type=int, default=8000,
                        help="CD-HIT memory limit in MB (default: 8000)")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads for CD-HIT (default: 4)")
    parser.add_argument("--output", required=True,
                        help="Output subsampled FASTA path")
    parser.add_argument("--report", default=None,
                        help="JSON report output path")

    args = parser.parse_args()

    ref_dir = args.ref_dir
    if not os.path.isdir(ref_dir):
        print(f"Error: reference directory not found: {ref_dir}", file=sys.stderr)
        sys.exit(1)

    taxonomy_weights = parse_taxonomy_weights(args.taxonomy_weights)

    # Step 1: Build sequence index
    print("Building sequence index...", file=sys.stderr)
    index, ordered = build_sequence_index(ref_dir)
    input_count = len(index)
    print(f"  Indexed {input_count} unique sequences", file=sys.stderr)

    if input_count == 0:
        print("Error: no sequences found in reference directory", file=sys.stderr)
        sys.exit(1)

    if input_count <= args.target_size:
        print(f"  Input ({input_count}) <= target ({args.target_size}), keeping all",
              file=sys.stderr)
        extract_selected(ref_dir, set(index.keys()), args.output)
        if args.report:
            group_breakdown = defaultdict(lambda: {"lse": 0, "one_to_one_ortholog": 0})
            for info in index.values():
                group_breakdown[info["group"]][info["category"]] += 1
            write_report(args.report, input_count, input_count, input_count,
                         dict(group_breakdown), taxonomy_weights, args.lse_weight,
                         args.cluster_identity)
        print(f"Wrote {input_count} sequences to {args.output}", file=sys.stderr)
        return

    # Step 2: Concatenate and run CD-HIT
    with tempfile.TemporaryDirectory() as tmpdir:
        concat_path = os.path.join(tmpdir, "all_refs.fa")
        cdhit_out = os.path.join(tmpdir, "cdhit_out")

        print("Concatenating reference sequences...", file=sys.stderr)
        dedup_count = concatenate_fasta(ref_dir, concat_path)
        print(f"  {dedup_count} unique sequences after deduplication", file=sys.stderr)

        print(f"Running CD-HIT at {args.cluster_identity} identity...", file=sys.stderr)
        representatives = run_cdhit(
            concat_path, cdhit_out,
            args.cluster_identity, args.cdhit_memory, args.threads
        )
        post_cluster_count = len(representatives)
        print(f"  {post_cluster_count} cluster representatives", file=sys.stderr)

        # Step 3: Taxonomy-weighted selection
        if post_cluster_count <= args.target_size:
            print(f"  Representatives ({post_cluster_count}) <= target ({args.target_size}), keeping all",
                  file=sys.stderr)
            selected = representatives
        else:
            print(f"Selecting {args.target_size} from {post_cluster_count} representatives...",
                  file=sys.stderr)
            selected = subsample(
                index, representatives, args.target_size,
                taxonomy_weights, args.lse_weight
            )

    # Step 4: Extract selected sequences from original files
    print("Extracting selected sequences...", file=sys.stderr)
    written = extract_selected(ref_dir, selected, args.output)
    print(f"Wrote {written} sequences to {args.output}", file=sys.stderr)

    # Step 5: Write report
    if args.report:
        group_breakdown = {}
        for header in selected:
            if header in index:
                info = index[header]
                group = info["group"]
                cat = info["category"]
                if group not in group_breakdown:
                    group_breakdown[group] = {"lse": 0, "one_to_one_ortholog": 0}
                group_breakdown[group][cat] += 1
        write_report(args.report, input_count, post_cluster_count, written,
                     group_breakdown, taxonomy_weights, args.lse_weight,
                     args.cluster_identity)
        print(f"Report written to {args.report}", file=sys.stderr)


if __name__ == "__main__":
    main()
