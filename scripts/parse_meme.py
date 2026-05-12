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

# bead -7cy: dual-mode (strict + lenient) concordance-stratified output.
# Cross-reference is at the column-vector level — two ClipKit passes on the
# same underlying alignment produce trimmed FASTAs whose columns are subsets
# of the same original; matching residue-vectors identifies "same column"
# without needing ClipKit logs. Per-OG counts drive ranking, per-site rows
# (one per MEME-positive site under either pass) feed sensitivity reports.
OG_FIELDS_DUAL = [
    "og_name",
    "n_strict_positive_sites",
    "n_lenient_positive_sites",
    "high_confidence_sites_n",
    "lenient_only_sites_n",
    "strict_only_sites_n",
    "alignment_robustness_index",
]
SITE_FIELDS_DUAL = [
    "og_name", "lenient_site", "strict_site",
    "strict_pvalue", "lenient_pvalue", "concordance_tier",
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


def _read_fasta(path: str) -> list[tuple[str, str]]:
    """Minimal FASTA reader returning list of (name, sequence) preserving order."""
    out: list[tuple[str, str]] = []
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(seq)))
                name = line[1:].split()[0]
                seq = []
            elif line:
                seq.append(line)
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def _fasta_columns(path: str, codon_aware: bool = False) -> list[tuple[str, ...]]:
    """Return per-column residue-vector tuples (deterministic order = order
    sequences appear in the FASTA).

    When ``codon_aware`` is True, columns are 3-character codons (the input
    sequences must be DNA with length divisible by 3). This is the mode used
    by stage 05's codon-level dual-mode parse: ClipKit -co preserves codon
    boundaries; MEME indexes by codon site; the column-vector comparison
    therefore must operate at codon granularity to align MEME's site space.
    """
    seqs = [s for _, s in _read_fasta(path)]
    if not seqs:
        return []
    L = len(seqs[0])
    if codon_aware:
        n_codons = L // 3
        return [tuple(s[3 * i:3 * (i + 1)] for s in seqs) for i in range(n_codons)]
    return [tuple(s[i] for s in seqs) for i in range(L)]


def map_strict_to_lenient(
    strict_fa: str,
    lenient_fa: str,
    codon_aware: bool = False,
) -> list[int | None]:
    """Map each strict-trimmed column (1-indexed) to its counterpart column
    in the lenient-trimmed alignment (1-indexed), or None if no match.

    Strict is a column-subset of lenient; equality of the full per-column
    residue-vector identifies the corresponding column. Linear-time scan
    with a forward pointer preserves alignment order.

    ``codon_aware=True`` operates at codon granularity (3 characters per
    column-vector entry) — required for matching MEME site indices when
    the trimmed alignments are DNA codon alignments.
    """
    strict_cols = _fasta_columns(strict_fa, codon_aware=codon_aware)
    lenient_cols = _fasta_columns(lenient_fa, codon_aware=codon_aware)
    mapping: list[int | None] = []
    j = 0
    for s_col in strict_cols:
        while j < len(lenient_cols) and lenient_cols[j] != s_col:
            j += 1
        if j < len(lenient_cols):
            mapping.append(j + 1)
            j += 1
        else:
            mapping.append(None)
    return mapping


def parse_one_dual(
    strict_json: str,
    strict_fa: str,
    lenient_json: str,
    lenient_fa: str,
    og_name: str,
    alpha: float = 0.05,
    codon_aware: bool = False,
) -> tuple[list[dict], dict]:
    """Dual-mode parser: read strict + lenient MEME JSONs, map columns
    across the two trimmed alignments, stratify each MEME-positive site
    by concordance tier.

    Returns:
        sites: list of per-site rows (one row per site MEME-positive
               under EITHER pass; tier identifies which)
        og:    OG-level summary with concordance counts + robustness index
    """
    with open(strict_json) as f:
        strict_meme = json.load(f)
    with open(lenient_json) as f:
        lenient_meme = json.load(f)
    strict_sites = extract_sites(strict_meme)
    lenient_sites = extract_sites(lenient_meme)

    def _is_pos(s: dict) -> bool:
        p, b = s["p_value"], s["beta_plus"]
        return (p == p) and p < alpha and (b == b) and b > 1.0

    mapping = map_strict_to_lenient(strict_fa, lenient_fa, codon_aware=codon_aware)

    # Build lookup: lenient-1indexed col → strict-1indexed col (inverse of mapping)
    lenient_to_strict: dict[int, int] = {}
    for s_idx_0, l_idx in enumerate(mapping):
        if l_idx is not None:
            lenient_to_strict[l_idx] = s_idx_0 + 1

    strict_positive_indices = {i + 1 for i, s in enumerate(strict_sites) if _is_pos(s)}
    lenient_positive_indices = {j + 1 for j, s in enumerate(lenient_sites) if _is_pos(s)}

    # Strict positives translated into lenient column space (drops strict
    # sites that have no lenient counterpart — should be empty in practice).
    strict_pos_in_lenient_coords = {
        mapping[i - 1] for i in strict_positive_indices if mapping[i - 1] is not None
    }

    high_confidence = strict_pos_in_lenient_coords & lenient_positive_indices
    lenient_only = lenient_positive_indices - strict_pos_in_lenient_coords
    strict_only = strict_positive_indices - {
        i for i in strict_positive_indices
        if mapping[i - 1] is not None and mapping[i - 1] in lenient_positive_indices
    }

    n_strict_pos = len(strict_positive_indices)
    n_lenient_pos = len(lenient_positive_indices)
    robustness = (len(high_confidence) / n_lenient_pos) if n_lenient_pos else 0.0

    og = {
        "og_name": og_name,
        "n_strict_positive_sites": n_strict_pos,
        "n_lenient_positive_sites": n_lenient_pos,
        "high_confidence_sites_n": len(high_confidence),
        "lenient_only_sites_n": len(lenient_only),
        "strict_only_sites_n": len(strict_only),
        "alignment_robustness_index": round(robustness, 6),
    }

    # Per-site rows: one row per MEME-positive site (under either pass).
    # Indexed by lenient column when possible (since lenient is the superset);
    # strict-only sites carry their strict column and None for lenient.
    rows: list[dict] = []
    # high_confidence + lenient_only — keyed on lenient col
    for l_idx in sorted(high_confidence | lenient_only):
        s_idx = lenient_to_strict.get(l_idx)
        tier = "high_confidence" if l_idx in high_confidence else "lenient_only"
        rows.append({
            "og_name": og_name,
            "lenient_site": l_idx,
            "strict_site": s_idx,
            "strict_pvalue": (strict_sites[s_idx - 1]["p_value"]
                              if s_idx is not None else None),
            "lenient_pvalue": lenient_sites[l_idx - 1]["p_value"],
            "concordance_tier": tier,
        })
    # strict_only — strict sites with no matching positive lenient
    for s_idx in sorted(strict_only):
        l_idx = mapping[s_idx - 1]
        rows.append({
            "og_name": og_name,
            "lenient_site": l_idx,
            "strict_site": s_idx,
            "strict_pvalue": strict_sites[s_idx - 1]["p_value"],
            "lenient_pvalue": (lenient_sites[l_idx - 1]["p_value"]
                               if l_idx is not None else None),
            "concordance_tier": "strict_only",
        })

    return rows, og


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

    # bead -7cy dual-mode: if all four are supplied, also emit concordance CSVs.
    ap.add_argument("--lenient-json", help="Lenient-trim MEME JSON for dual-mode")
    ap.add_argument("--strict-fa", help="Strict-trim protein FASTA (--lenient-json required)")
    ap.add_argument("--lenient-fa", help="Lenient-trim protein FASTA (--lenient-json required)")
    ap.add_argument("--out-og-concordance",
                    help="Per-OG concordance CSV (default: derived from --out-og)")
    ap.add_argument("--out-sites-concordance",
                    help="Per-site concordance CSV (default: derived from --out-sites)")
    ap.add_argument("--codon", action="store_true",
                    help="Treat --strict-fa / --lenient-fa as DNA codon "
                         "alignments (3-char column-vectors); required when "
                         "the dual-mode inputs are codon files (ClipKit -co)")
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

    # Optional dual-mode concordance pass (bead -7cy). Requires the
    # lenient JSON + both trimmed FASTAs. Single-OG only (--json mode):
    # batch concordance is left to a future iteration if useful.
    if args.lenient_json:
        if not (args.json and args.strict_fa and args.lenient_fa and args.og_name):
            print("ERROR: --lenient-json requires --json + --strict-fa + "
                  "--lenient-fa + --og-name", file=sys.stderr)
            return 1
        dual_sites, dual_og = parse_one_dual(
            strict_json=args.json, strict_fa=args.strict_fa,
            lenient_json=args.lenient_json, lenient_fa=args.lenient_fa,
            og_name=args.og_name, alpha=args.alpha,
            codon_aware=args.codon,
        )
        out_og_dual = args.out_og_concordance or args.out_og.replace(
            ".csv", "_concordance.csv")
        Path(out_og_dual).parent.mkdir(parents=True, exist_ok=True)
        with open(out_og_dual, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=OG_FIELDS_DUAL)
            w.writeheader()
            w.writerow(dual_og)
        out_sites_dual = (args.out_sites_concordance
                          or (args.out_sites and args.out_sites.replace(
                              ".csv", "_concordance.csv")))
        if out_sites_dual:
            Path(out_sites_dual).parent.mkdir(parents=True, exist_ok=True)
            with open(out_sites_dual, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=SITE_FIELDS_DUAL)
                w.writeheader()
                w.writerows(dual_sites)
        print(f"Wrote concordance: high={dual_og['high_confidence_sites_n']} "
              f"lenient_only={dual_og['lenient_only_sites_n']} "
              f"strict_only={dual_og['strict_only_sites_n']} "
              f"robustness={dual_og['alignment_robustness_index']:.3f} "
              f"to {out_og_dual}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
