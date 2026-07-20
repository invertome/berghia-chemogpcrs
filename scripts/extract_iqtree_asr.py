#!/usr/bin/env python3
"""extract_iqtree_asr.py — Extract a single ancestral sequence from IQ-TREE
.state output.

IQ-TREE 3 ``--ancestral`` writes a tab-separated .state file with columns:
    Node    Site    State   p_A  p_C  p_G  p_T  ...   (or AA states for AA models)

For each node + site, the most-probable state is the column with the
maximum posterior probability. We pick the argmax per site for the
requested node and write a single FASTA record.

Bead -j44: this is the ASR extractor used by 05_selective_pressure_and_asr.sh
when ASR_BACKEND=iqtree (the new default).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def extract_sites_for_node(state_file: str, node_name: str):
    """Yield (site, state) pairs for node ``node_name`` in ``state_file``."""
    header_seen = False
    with open(state_file) as f:
        for line in f:
            if line.startswith("# Node") or line.startswith("#"):
                continue
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if not header_seen:
                # First non-comment line is the header (Node, Site, State, ...)
                header_seen = True
                continue
            if len(parts) < 3:
                continue
            if parts[0] != node_name:
                continue
            try:
                site_num = int(parts[1])
            except ValueError:
                continue
            best_state = parts[2]  # IQ-TREE writes the marginal ML state here
            yield site_num, best_state


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("state_file", help="IQ-TREE .state file (from --ancestral)")
    ap.add_argument("node_name", help="Internal node name to extract "
                                      "(the .state lookup key, e.g. 'Node5')")
    ap.add_argument("output_fasta", help="Output FASTA path")
    # The .state lookup key and the emitted record id are two different jobs.
    # IQ-TREE numbers internal nodes per tree, so 'Node5' is only unique WITHIN
    # one orthogroup's tree; two orthogroups would both emit '>Node5' for
    # completely different ancestral sequences, and stage 08 concatenates every
    # *_asr.fa into one extraction source. Callers therefore pass a globally
    # unique id (the '{orthogroup}_{node}' form already used for the filename),
    # while the lookup below keeps using the bare label IQ-TREE actually wrote.
    # Optional, defaulting to the bare node name, so the three-positional CLI
    # stays backward compatible.
    ap.add_argument("--record-id", default=None,
                    help="Record id for the emitted FASTA header "
                         "(default: node_name). Use an orthogroup-namespaced "
                         "id such as 'OG0000001_Node5' when the output is "
                         "pooled across orthogroups.")
    args = ap.parse_args()

    record_id = args.record_id if args.record_id is not None else args.node_name

    sites = list(extract_sites_for_node(args.state_file, args.node_name))
    if not sites:
        print(f"ERROR: no sites found for node {args.node_name!r} in {args.state_file}",
              file=sys.stderr)
        return 1

    sites.sort(key=lambda x: x[0])
    seq = "".join(s for _, s in sites)
    Path(args.output_fasta).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_fasta, "w") as f:
        f.write(f">{record_id}\n{seq}\n")
    print(f"Wrote {len(seq)}-residue ancestral sequence for {args.node_name} "
          f"as {record_id!r} -> {args.output_fasta}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
