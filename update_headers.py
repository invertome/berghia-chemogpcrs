#!/usr/bin/env python3
# update_headers.py
# Purpose: Update FASTA headers with short IDs and create an ID mapping file.
# Inputs: FASTA file ($1), output ID map CSV ($2)
# Outputs: Updated FASTA file with short IDs, ID map CSV
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import sys

fasta_file = sys.argv[1]  # Input FASTA file path
id_map_file = sys.argv[2]  # Output CSV file path for ID mapping

# Initialize DataFrame to store ID mappings
id_map = pd.DataFrame(columns=['short_id', 'original_id', 'taxid'])
counter = 1  # Counter for generating short IDs

# Read FASTA and create ID map
with open(fasta_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            original_id = line.split()[0][1:]  # Extract original ID (remove '>')
            short_id = f'ref_{counter}'        # Generate short ID (e.g., ref_1)
            taxid = original_id.split('_')[0]  # Extract TaxID from original ID
            id_map = pd.concat([id_map, pd.DataFrame([{'short_id': short_id, 'original_id': original_id, 'taxid': taxid}])], ignore_index=True)
            counter += 1

# Save ID map to CSV
id_map.to_csv(id_map_file, index=False)

# Write updated FASTA with short IDs
with open(fasta_file, 'r') as f, open(fasta_file.replace('.fa', '_updated.fa'), 'w') as out:
    for line in f:
        if line.startswith('>'):
            original_id = line.split()[0][1:]
            short_id = id_map[id_map['original_id'] == original_id]['short_id'].values[0]
            out.write(f'>{short_id}\n')
        else:
            out.write(line)
