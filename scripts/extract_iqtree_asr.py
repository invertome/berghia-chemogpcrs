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
import os
import sys
import tempfile
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



def _write_fasta_atomic(output_fasta: str, payload: str) -> None:
    """Publish the record by rename so a reader never sees a partial file.

    08_structural_analysis.sh globs this directory (`find ... -name '*_asr.fa'`)
    to pool ancestral sequences. Writing the destination directly truncates it
    on open, so a task killed mid-write leaves a short-but-present FASTA that
    the glob picks up and treats as complete -- the failure is silent, because a
    truncated FASTA still parses. Staging to a unique temp in the destination
    directory (same filesystem, or os.replace is not atomic) and renaming makes
    the file appear whole or not at all. The stage is a SLURM array, so the tag
    carries job and task identity as well as the pid: a pid is unique only
    within a node, and this directory is shared across nodes.
    """
    out = Path(output_fasta)
    tag = "{}.{}.{}".format(os.environ.get("SLURM_JOB_ID", "nojob"),
                            os.environ.get("SLURM_ARRAY_TASK_ID", "0"),
                            os.getpid())
    fd, tmp = tempfile.mkstemp(dir=str(out.parent), prefix=out.name + ".",
                               suffix="." + tag + ".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(payload)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, out)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


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
    _write_fasta_atomic(args.output_fasta, f">{record_id}\n{seq}\n")
    print(f"Wrote {len(seq)}-residue ancestral sequence for {args.node_name} "
          f"as {record_id!r} -> {args.output_fasta}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
