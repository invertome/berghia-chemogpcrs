#!/usr/bin/env python3
"""parse_busted.py — Parse HyPhy BUSTED / BUSTED-S / BUSTED-MH JSON output.

Bead -urk follow-up. BUSTED is a gene-wide test for episodic positive
selection. BUSTED-S adds synonymous-rate variation (Wisotsky 2020 MBE
37:2430). BUSTED-MH adds multi-nucleotide-substitution correction
(Lucaci 2023 MBE 40:msad150) and is the most stringent / credible
gene-wide signal.

Per-OG output schema (one row per JSON file parsed):
    og_name, model_variant (S | MH | plain), p_value, lrt,
    omega_low, prop_low, omega_mid, prop_mid, omega_high, prop_high,
    is_significant

The unconstrained model fits a 3-class omega distribution; we extract
the high-omega class as the "positive selection" component.

Usage:
    python parse_busted.py --json results/.../OGxx_busted_s.json \\
        --og-name OGxx --variant S --out per_og.csv
    # or batch:
    python parse_busted.py --batch results/selective_pressure/ \\
        --variant S --out busted_s_results.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Iterable, Optional

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
    also keeps staging files out of the consumers' "*_<variant>.csv" globs.
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
    "og_name",
    "model_variant",
    "p_value",
    "lrt",
    "omega_low", "prop_low",
    "omega_mid", "prop_mid",
    "omega_high", "prop_high",
    "is_significant",
]


def _safe_float(x, default: float = float("nan")) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def extract_busted(busted_json: dict) -> dict:
    """Pull the gene-wide test result + omega rate distribution from a BUSTED JSON.

    Returns dict with: p_value, lrt, omega_(low|mid|high), prop_(low|mid|high).
    """
    test = busted_json.get("test results", {}) or busted_json.get("Test results", {})
    p = _safe_float(test.get("p-value", test.get("p value", float("nan"))))
    lrt = _safe_float(test.get("LRT", test.get("lrt", float("nan"))))

    # The "fits" key holds different model fits. The unconstrained model has
    # a 3-class omega distribution under "rate distributions" -> "Test"
    fits = busted_json.get("fits", {})
    unconstr = fits.get("Unconstrained model", {})
    rate_dist = (
        unconstr.get("Rate Distributions", {}).get("Test", [])
        or unconstr.get("rate distributions", {}).get("Test", [])
        or []
    )
    # Sort by omega ascending
    classes = []
    for rc in rate_dist:
        if isinstance(rc, dict):
            o = _safe_float(rc.get("omega", rc.get("Omega", float("nan"))))
            pr = _safe_float(rc.get("proportion", rc.get("weight", float("nan"))))
            classes.append((o, pr))
        elif isinstance(rc, (list, tuple)) and len(rc) >= 2:
            classes.append((_safe_float(rc[0]), _safe_float(rc[1])))
    classes.sort(key=lambda x: x[0])
    if not classes:
        return {
            "p_value": p, "lrt": lrt,
            "omega_low": float("nan"), "prop_low": float("nan"),
            "omega_mid": float("nan"), "prop_mid": float("nan"),
            "omega_high": float("nan"), "prop_high": float("nan"),
        }
    # Always preserve the actual lowest + actual highest classes, even when
    # only 1 or 2 are fitted. The "mid" slot is NaN-padded if missing.
    low = classes[0]
    high = classes[-1]
    if len(classes) == 1:
        mid = (float("nan"), float("nan"))
    elif len(classes) == 2:
        mid = (float("nan"), float("nan"))
    else:
        mid = classes[1]

    return {
        "p_value": p,
        "lrt": lrt,
        "omega_low": low[0], "prop_low": low[1],
        "omega_mid": mid[0], "prop_mid": mid[1],
        "omega_high": high[0], "prop_high": high[1],
    }


def parse_one(json_path: str, og_name: str, variant: str,
              alpha: float = 0.05) -> dict:
    with open(json_path) as f:
        d = json.load(f)
    stats = extract_busted(d)
    is_sig = (not (stats["p_value"] != stats["p_value"])  # not NaN
              and stats["p_value"] < alpha)
    return {
        "og_name": og_name,
        "model_variant": variant,
        **stats,
        "is_significant": int(bool(is_sig)),
    }


def discover_batch(directory: str, variant: str) -> Iterable[tuple[str, str]]:
    """Yield (og_name, json_path) for files like '<og>_busted_<variant>.json'."""
    suffix = f"_busted_{variant.lower()}.json"
    for p in sorted(Path(directory).glob(f"*{suffix}")):
        og = p.name[: -len(suffix)]
        yield og, str(p)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--json", help="Single BUSTED JSON path (with --og-name)")
    src.add_argument("--batch",
                     help="Directory of <og>_busted_{variant}.json files")
    ap.add_argument("--og-name", help="Required when --json is used")
    ap.add_argument("--variant", choices=["S", "MH", "plain"], default="S",
                    help="Model variant tag for the model_variant column")
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="Significance threshold (default 0.05)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    if args.json:
        if not args.og_name:
            print("ERROR: --og-name required with --json", file=sys.stderr)
            return 1
        rows.append(parse_one(args.json, args.og_name, args.variant, args.alpha))
    else:
        for og, path in discover_batch(args.batch, args.variant):
            try:
                rows.append(parse_one(path, og, args.variant, args.alpha))
            except Exception as e:
                print(f"WARN: failed to parse {path}: {e}", file=sys.stderr)

    write_csv_atomic(args.out, FIELDNAMES, rows)
    n_sig = sum(1 for r in rows if r["is_significant"])
    print(f"Wrote {len(rows)} BUSTED-{args.variant} records to {args.out} "
          f"({n_sig} significant at alpha={args.alpha})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
