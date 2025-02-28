#!/usr/bin/env python3
# test_phyloformer_models.py
# Purpose: Test Phyloformer models (LG, WAG, JTT) and select the best based on AIC.
# Inputs: Input FASTA ($1), output prefix ($2), number of threads ($3)
# Outputs: Best Phyloformer tree (${output_prefix}.tre)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import subprocess
import sys
import os

input_fasta = sys.argv[1]  # Input alignment FASTA
output_prefix = sys.argv[2]  # Output prefix for tree files
threads = sys.argv[3]      # Number of threads for Phyloformer

# Supported Phyloformer models
models = ['LG', 'WAG', 'JTT']
results = {}

# Test each model
for model in models:
    output_tree = f"{output_prefix}_{model}.tre"
    cmd = [
        "phyloformer",
        "-i", input_fasta,
        "-o", output_tree,
        "--model", model,
        "--threads", threads
    ]
    process = subprocess.run(cmd, capture_output=True, text=True)
    if process.returncode != 0:
        print(f"Error running Phyloformer with {model}: {process.stderr}", file=sys.stderr)
        continue
    
    # Extract log-likelihood from output (assumes Phyloformer logs it similarly to IQ-TREE)
    log_likelihood = float(process.stdout.split('Log-likelihood: ')[1].split()[0]) if 'Log-likelihood' in process.stdout else float('-inf')
    # Simplified AIC: 2 * num_params + (-2 * log_likelihood), assuming same params for simplicity
    aic = -2 * log_likelihood  # Placeholder, adjust if Phyloformer provides param count
    results[model] = (aic, output_tree)

# Select best model based on minimum AIC
if results:
    best_model, (best_aic, best_tree) = min(results.items(), key=lambda x: x[1][0])
    print(f"Best model: {best_model} with AIC: {best_aic}", file=sys.stderr)
    os.rename(best_tree, f"{output_prefix}.tre")
    # Clean up other model trees
    for model, (_, tree) in results.items():
        if tree != f"{output_prefix}.tre":
            os.remove(tree)
else:
    print("Error: No models successfully tested.", file=sys.stderr)
    sys.exit(1)
