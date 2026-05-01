#!/usr/bin/env python3
"""compute_lrt.py — Likelihood Ratio Test for nested PAML codeml model pairs.

Bead -cio: previous version hardcoded df=1 (only correct for the simplest
nested pair) and the lnL regex required a leading minus, missing rare positive
log-likelihoods. Both fixed here.

Common df values for PAML codeml nested comparisons:
  M1a vs M2a (site-class 'positive selection')        : df=2
  M1a-NULL vs M2a (site-class 'positive', boundary)   : df=1, use boundary mode
  M7 vs M8 (beta vs beta + omega>1)                   : df=2
  M8a vs M8 (omega=1 boundary vs beta + omega>1)      : df=1, use boundary mode
  Branch-site Test1 (CmC vs CmD) common             : df=2
  Branch-site Test2 (alt vs null at boundary)        : df=1, use boundary mode

For boundary tests (M1a-NULL vs M2a, M8a vs M8, branch-site Test2), the
asymptotic null distribution is a 50:50 mixture of chi2(0) and chi2(1)
(Yang 2007). Use --boundary to apply this; effectively halves the p-value.
"""
import argparse
import re
import sys

from scipy.stats import chi2


def extract_lnL(filepath):
    """Extract log-likelihood value from PAML output (sign-optional regex)."""
    with open(filepath) as f:
        for line in f:
            if 'lnL' not in line:
                continue
            # Pattern 1: lnL(ntime: X np: Y): Z = <value>
            m = re.search(r'lnL\([^)]+\):\s*\S+\s*=\s*(-?\d+\.\d+)', line)
            if m:
                return float(m.group(1))
            # Pattern 2: lnL = <value>
            m = re.search(r'lnL\s*=\s*(-?\d+\.\d+)', line)
            if m:
                return float(m.group(1))
            # Pattern 3: any signed float after 'lnL'
            m = re.search(r'lnL.*?(-?\d+\.\d+)', line)
            if m:
                return float(m.group(1))
    raise ValueError(f"No lnL value found in {filepath}")


def lrt_pvalue(chi2_stat: float, df: int, boundary: bool = False) -> float:
    """Return p-value for a likelihood ratio test.

    If ``boundary`` is True, the null distribution is a 50:50 mixture of
    chi2(0) and chi2(df-1) (= chi2(0) point mass at 0 plus chi2(df-1)),
    appropriate for tests where the alternative parameter is at a
    boundary of the parameter space (Self & Liang 1987; Yang 2007).
    """
    if chi2_stat <= 0:
        return 1.0
    if not boundary:
        return float(chi2.sf(chi2_stat, df))
    # Boundary case: 0.5 * 0 + 0.5 * chi2(df) = 0.5 * chi2.sf(stat, df)
    return float(0.5 * chi2.sf(chi2_stat, df))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n', 1)[0])
    ap.add_argument("base", help="Base name for orthogroup")
    ap.add_argument("null_file", help="PAML output for null model")
    ap.add_argument("alt_file", help="PAML output for alternative model")
    ap.add_argument("output_file", help="Output CSV (appended)")
    ap.add_argument("--df", type=int, default=2,
                    help="Degrees of freedom for chi2 (default 2; use 1 for "
                         "single-parameter boundary tests like M8a vs M8 or "
                         "branch-site Test2)")
    ap.add_argument("--boundary", action="store_true",
                    help="Apply 50:50 chi2 mixture for boundary tests "
                         "(Yang 2007). Effectively halves the p-value.")
    args = ap.parse_args()

    try:
        lnL_null = extract_lnL(args.null_file)
        lnL_alt = extract_lnL(args.alt_file)
    except FileNotFoundError as e:
        print(f"Error: File not found for {args.base}: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: Could not parse lnL for {args.base}: {e}", file=sys.stderr)
        sys.exit(1)

    chi2_stat = 2 * (lnL_alt - lnL_null)
    if chi2_stat < 0:
        print(f"Warning: Negative chi2 statistic for {args.base} "
              f"(lnL_null={lnL_null}, lnL_alt={lnL_alt}). "
              f"This may indicate model fitting issues or swapped null/alt files",
              file=sys.stderr)

    p_value = lrt_pvalue(chi2_stat, df=args.df, boundary=args.boundary)

    with open(args.output_file, 'a') as f:
        f.write(f"{args.base},{chi2_stat},{p_value},{args.df},{int(args.boundary)}\n")


if __name__ == "__main__":
    main()
