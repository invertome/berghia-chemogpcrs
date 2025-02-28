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
try:
    with open(null_file) as f:
        lnL_null = float(next(line for line in f if 'lnL' in line).split()[4])
    with open(alt_file) as f:
        lnL_alt = float(next(line for line in f if 'lnL' in line).split()[4])
except Exception as e:
    print(f"Error processing {base}: {e}", file=sys.stderr)
    lnL_null, lnL_alt = 0, 0

# Compute LRT statistic and p-value
chi2_stat = 2 * (lnL_alt - lnL_null)
p_value = chi2.sf(chi2_stat, 1) if chi2_stat > 0 else 1.0

# Append results to CSV
with open(output_file, 'a') as f:
    f.write(f"{base},{chi2_stat},{p_value}\n")
