#!/usr/bin/env python3
"""array_plan.py — size a SLURM job array from the orthogroup manifest (bead fxx).

Stages 04/05 carry a hardcoded ``#SBATCH --array=0-999%50``. On a full-557
relaunch the manifest exceeds 1000 orthogroups, so every orthogroup at index
>=1000 falls outside the array. ``assert_array_covers_manifest`` already makes
that fail LOUDLY rather than silently truncating; this planner supplies the
correct range so the submit wrapper can override the directive
(``sbatch`` CLI flags take precedence over ``#SBATCH`` lines).

Pure functions + a thin CLI so the range math is unit-tested. The CLI prints the
array spec on stdout, or exits non-zero with a chunk plan when the manifest
exceeds the cluster's MaxArraySize (Unity: 10001 => highest legal index 10000).
"""
from __future__ import annotations

from typing import List, Optional, Tuple


def array_spec(n_tasks: int, throttle: Optional[int] = None) -> str:
    """``--array`` spec covering task indices ``0..n_tasks-1``.

    ``throttle`` appends the ``%K`` concurrent-task limit when given.
    """
    if n_tasks is None or int(n_tasks) <= 0:
        raise ValueError(f"n_tasks must be > 0 (got {n_tasks!r})")
    n = int(n_tasks)
    spec = f"0-{n - 1}"
    return f"{spec}%{int(throttle)}" if throttle else spec


def fits_single_array(n_tasks: int, max_array_size: int) -> bool:
    """True when indices ``0..n_tasks-1`` fit under the cluster ceiling.

    ``max_array_size`` is Slurm's MaxArraySize (the array SIZE), so the highest
    legal index is ``max_array_size - 1`` and ``n_tasks == max_array_size`` fits.
    """
    return int(n_tasks) <= int(max_array_size)


def chunk_ranges(n_tasks: int, max_array_size: int) -> List[Tuple[int, int]]:
    """Contiguous, non-overlapping global index ranges each within the ceiling.

    Returns ``[(start, end), ...]`` covering ``0..n_tasks-1``. A manifest that
    fits yields a single range. Chunks beyond the first need a per-chunk index
    OFFSET in the stage (the array task id alone cannot exceed the ceiling), so
    the CLI treats multi-chunk plans as an operator decision rather than
    submitting them silently.
    """
    n = int(n_tasks)
    size = int(max_array_size)
    if n <= 0:
        raise ValueError(f"n_tasks must be > 0 (got {n_tasks!r})")
    if size <= 0:
        raise ValueError(f"max_array_size must be > 0 (got {max_array_size!r})")
    out: List[Tuple[int, int]] = []
    start = 0
    while start < n:
        end = min(start + size - 1, n - 1)
        out.append((start, end))
        start = end + 1
    return out


def main(argv=None) -> int:
    import argparse
    import sys

    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--n-tasks", type=int, required=True,
                    help="number of array tasks (orthogroup count)")
    ap.add_argument("--throttle", type=int, default=0,
                    help="max concurrent tasks (%%K); 0 = no throttle")
    ap.add_argument("--max-array-size", type=int, default=10001,
                    help="cluster MaxArraySize (default 10001)")
    args = ap.parse_args(argv)

    if not fits_single_array(args.n_tasks, args.max_array_size):
        plan = chunk_ranges(args.n_tasks, args.max_array_size)
        print(
            f"ERROR: {args.n_tasks} tasks exceed MaxArraySize={args.max_array_size}. "
            f"Chunked submission needs a per-chunk index offset in the stage "
            f"(the array task id cannot exceed {args.max_array_size - 1}). Plan: "
            + ", ".join(f"{s}-{e}" for s, e in plan),
            file=sys.stderr,
        )
        return 3

    print(array_spec(args.n_tasks, args.throttle or None))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
