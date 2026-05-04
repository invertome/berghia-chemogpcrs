#!/usr/bin/env python3
"""parse_meme.py — Parse HyPhy MEME JSON output (per-site episodic selection).

Bead -urk follow-up. MEME (Murrell et al. 2012 PLOS Genet 8:e1002764)
tests each codon site for episodic positive selection on a subset of
branches. It is the site-level companion to BUSTED's gene-wide test
and aBSREL's branch-level test.

Two output modes:
  1. Per-site CSV (one row per site): site, alpha, beta_plus, p_value,
     prop_beta_plus, n_branches_episodic, is_episodic
  2. Per-OG aggregate (one row per JSON parsed): og_name,
     n_sites_total, n_sites_episodic, fraction_episodic, mean_beta_plus

Use --batch to scan a directory of <og>_meme.json files and emit only
the per-OG aggregate (suitable for joining into rank_candidates).

Usage:
    # Single file, both per-site and per-OG outputs:
    python parse_meme.py --json results/.../OGxx_meme.json \\
        --og-name OGxx --out-sites per_site.csv --out-og per_og.csv

    # Batch mode (per-OG aggregate only):
    python parse_meme.py --batch results/selective_pressure/ --out-og per_og.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Iterable

OG_FIELDS = [
    "og_name",
    "n_sites_total",
    "n_sites_episodic",
    "fraction_episodic",
    "mean_beta_plus",
    "median_p_value",
]
SITE_FIELDS = [
    "og_name", "site", "alpha", "beta_plus", "prop_beta_plus",
    "p_value", "n_branches_episodic", "is_episodic",
]


def _safe_float(x, default: float = float("nan")) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def extract_sites(meme_json: dict) -> list[dict]:
    """Pull per-site rows from a MEME JSON."""
    # MEME stores per-site results under "MLE" -> "content" -> "0" -> rows.
    # Headers under "MLE" -> "headers" describe the columns; standard layout:
    #   alpha, beta-, p-, beta+, p+, LRT, p-value, branches, post-mean
    mle = meme_json.get("MLE", {})
    content = mle.get("content", {})
    headers = mle.get("headers", [])
    if not content or not headers:
        return []
    # Index of relevant headers
    def hidx(label_substrings):
        for i, h in enumerate(headers):
            name = (h[0] if isinstance(h, (list, tuple)) else str(h)).lower()
            for s in label_substrings:
                if s in name:
                    return i
        return None

    i_alpha = hidx(["alpha"])
    i_beta_plus = hidx(["&beta;+", "beta+", "beta_plus"])
    i_prop_plus = hidx(["p+"])
    i_pvalue = hidx(["p-value"])
    i_branches = hidx(["# branches", "branches"])

    # Content "0" is partition 0 (most MEME runs have one)
    rows_in = []
    for k, v in content.items():
        if isinstance(v, list):
            rows_in = v
            break
    if not rows_in and isinstance(content, list):
        rows_in = content

    out = []
    for site_idx, row in enumerate(rows_in, start=1):
        if not isinstance(row, list):
            continue
        out.append({
            "site": site_idx,
            "alpha": _safe_float(row[i_alpha]) if i_alpha is not None else float("nan"),
            "beta_plus": _safe_float(row[i_beta_plus]) if i_beta_plus is not None else float("nan"),
            "prop_beta_plus": _safe_float(row[i_prop_plus]) if i_prop_plus is not None else float("nan"),
            "p_value": _safe_float(row[i_pvalue]) if i_pvalue is not None else float("nan"),
            "n_branches_episodic": _safe_float(row[i_branches]) if i_branches is not None else float("nan"),
        })
    return out


def aggregate(sites: list[dict], alpha: float) -> dict:
    """Reduce per-site rows to a single OG-level summary."""
    n_total = len(sites)
    pvals = [s["p_value"] for s in sites if not (s["p_value"] != s["p_value"])]
    episodic = [s for s in sites
                if (s["p_value"] == s["p_value"]) and s["p_value"] < alpha
                and (s["beta_plus"] == s["beta_plus"]) and s["beta_plus"] > 1.0]
    n_ep = len(episodic)
    frac = (n_ep / n_total) if n_total else 0.0
    mean_bp = (sum(s["beta_plus"] for s in episodic) / n_ep) if n_ep else float("nan")
    median_p = float("nan")
    if pvals:
        pvals_sorted = sorted(pvals)
        m = len(pvals_sorted)
        median_p = (pvals_sorted[m // 2] if m % 2 else
                    0.5 * (pvals_sorted[m // 2 - 1] + pvals_sorted[m // 2]))
    return {
        "n_sites_total": n_total,
        "n_sites_episodic": n_ep,
        "fraction_episodic": round(frac, 6),
        "mean_beta_plus": mean_bp,
        "median_p_value": median_p,
    }


def parse_one(json_path: str, og_name: str, alpha: float = 0.05):
    with open(json_path) as f:
        d = json.load(f)
    sites = extract_sites(d)
    for s in sites:
        s["og_name"] = og_name
        s["is_episodic"] = int(bool(
            (s["p_value"] == s["p_value"]) and s["p_value"] < alpha
            and (s["beta_plus"] == s["beta_plus"]) and s["beta_plus"] > 1.0
        ))
    og = {"og_name": og_name, **aggregate(sites, alpha)}
    return sites, og


def discover_batch(directory: str) -> Iterable[tuple[str, str]]:
    suffix = "_meme.json"
    for p in sorted(Path(directory).glob(f"*{suffix}")):
        og = p.name[: -len(suffix)]
        yield og, str(p)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--json")
    src.add_argument("--batch")
    ap.add_argument("--og-name")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--out-sites", help="Per-site CSV (optional)")
    ap.add_argument("--out-og", required=True, help="Per-OG aggregate CSV")
    args = ap.parse_args()

    site_rows = []
    og_rows = []
    if args.json:
        if not args.og_name:
            print("ERROR: --og-name required with --json", file=sys.stderr)
            return 1
        sites, og = parse_one(args.json, args.og_name, args.alpha)
        site_rows.extend(sites)
        og_rows.append(og)
    else:
        for og_name, path in discover_batch(args.batch):
            try:
                sites, og = parse_one(path, og_name, args.alpha)
                site_rows.extend(sites)
                og_rows.append(og)
            except Exception as e:
                print(f"WARN: failed to parse {path}: {e}", file=sys.stderr)

    Path(args.out_og).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_og, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OG_FIELDS)
        w.writeheader()
        w.writerows(og_rows)

    if args.out_sites:
        Path(args.out_sites).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out_sites, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=SITE_FIELDS)
            w.writeheader()
            w.writerows(site_rows)

    print(f"Wrote {len(og_rows)} OG-level + {len(site_rows)} site-level "
          f"MEME records to {args.out_og}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
