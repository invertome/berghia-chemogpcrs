#!/usr/bin/env python3
# compute_lrt.py
# Purpose: Compute Likelihood Ratio Test (LRT) statistics for selective pressure analysis from PAML output.
# Inputs: Base name ($1), null model output ($2), alternative model output ($3), output CSV ($4)
# Outputs: Appends LRT results to CSV file
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from scipy.stats import chi2

base = sys.argv[1]  # Base name for orthogroup
null_file = sys.argv[2]  # PAML output for null model
alt_file = sys.argv[3]   # PAML output for alternative model
output_file = sys.argv[4]  # Output CSV file

# Extract log-likelihoods from PAML output
# PAML lnL line format varies but typically contains: lnL(ntime: X np: Y): Z = <value>
# Or: lnL = <value> or lnL: <value>
import re

def extract_lnL(filepath):
    """Extract log-likelihood value from PAML output file."""
    with open(filepath) as f:
        for line in f:
            if 'lnL' in line:
                # Try multiple patterns for PAML output format variations
                # Pattern 1: lnL(ntime: X np: Y): Z = -12345.67
                match = re.search(r'lnL\([^)]+\):\s*\S+\s*=\s*([-\d.]+)', line)
                if match:
                    return float(match.group(1))
                # Pattern 2: lnL = -12345.67
                match = re.search(r'lnL\s*=\s*([-\d.]+)', line)
                if match:
                    return float(match.group(1))
                # Pattern 3: Just find any floating point number after lnL
                match = re.search(r'lnL.*?([-]\d+\.\d+)', line)
                if match:
                    return float(match.group(1))
    raise ValueError(f"No lnL value found in {filepath}")

try:
    lnL_null = extract_lnL(null_file)
    lnL_alt = extract_lnL(alt_file)
except Exception as e:
    print(f"Error processing {base}: {e}", file=sys.stderr)
    lnL_null, lnL_alt = 0, 0

# Compute LRT statistic and p-value
chi2_stat = 2 * (lnL_alt - lnL_null)
p_value = chi2.sf(chi2_stat, 1) if chi2_stat > 0 else 1.0

# Append results to CSV
with open(output_file, 'a') as f:
    f.write(f"{base},{chi2_stat},{p_value}\n")
