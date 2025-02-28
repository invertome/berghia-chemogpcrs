#!/usr/bin/env python3
# update_headers.py
# Purpose: Update FASTA headers with short IDs from an ID map.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import sys

fasta_file = sys.argv[1]
id_map_file = sys.argv[2]

id_map = pd.DataFrame(columns=['short_id', 'original_id', 'taxid'])
counter = 1
with open(fasta_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            original_id = line.split()[0][1:]
            short_id = f'ref_{counter}'
            taxid = original_id.split('_')[0]
            id_map = pd.concat([id_map, pd.DataFrame([{'short_id': short_id, 'original_id': original_id, 'taxid': taxid}])], ignore_index=True)
            counter += 1
id_map.to_csv(id_map_file, index=False)

with open(fasta_file, 'r') as f, open(fasta_file.replace('.fa', '_updated.fa'), 'w') as out:
    for line in f:
        if line.startswith('>'):
            original_id = line.split()[0][1:]
            short_id = id_map[id_map['original_id'] == original_id]['short_id'].values[0]
            out.write(f'>{short_id}\n')
        else:
            out.write(line)
