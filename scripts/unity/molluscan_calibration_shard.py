#!/usr/bin/env python3
"""molluscan_calibration_shard.py — deterministic FASTA shard extractor.

Each array task writes ONLY its own shard file, so there is no shared temp path
and no truncate/publish race between concurrent tasks (the failure mode bead
`ih5u` records for run_selection_stack.sh). Record i goes to shard i % nshards,
so the shard set is a partition of the input regardless of how many tasks
actually run, and re-running a task reproduces byte-identical output.
"""
from __future__ import annotations

import argparse
import os
import sys


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--nshards", type=int, required=True)
    p.add_argument("--shard", type=int, required=True, help="0-based")
    a = p.parse_args()
    if not 0 <= a.shard < a.nshards:
        print(f"FATAL: shard {a.shard} out of range [0,{a.nshards})", file=sys.stderr)
        return 2

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    n_in = n_out = 0
    keep = False
    with open(a.input) as fh, open(a.out + ".tmp", "w") as out:
        for line in fh:
            if line.startswith(">"):
                keep = (n_in % a.nshards) == a.shard
                n_in += 1
                if keep:
                    n_out += 1
            if keep:
                out.write(line)
    os.replace(a.out + ".tmp", a.out)
    print(f"[shard] {a.shard}/{a.nshards}: {n_out} of {n_in} records -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
