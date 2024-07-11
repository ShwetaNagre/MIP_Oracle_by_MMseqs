# -*- coding: utf-8 -*-
"""
@author: Shweta Nagre
"""

import pandas as pd

# Define file paths
mmseqs_result_file = 'mmseqs2_result'
whole_region_file = 'whole_region.txt'
no_hits_file = 'no_hits_ids.txt'

# Read the MMseqs result file to get the IDs with hits
hits_ids = pd.read_csv(mmseqs_result_file, sep='\t', header=None)
hits_ids_set = set(hits_ids[0].tolist())

# Read the whole_region.txt file to get all query IDs
query_ids = []
with open(whole_region_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            query_id = line.strip()[1:]  # Remove '>' for matching purposes
            query_ids.append(query_id)

# Identify the IDs that did not get any hits
no_hits_ids = [query_id for query_id in query_ids if query_id not in hits_ids_set]

# Write the IDs with no hits to a file
with open(no_hits_file, 'w') as f:
    for no_hit_id in no_hits_ids:
        f.write('>' + no_hit_id + '\n')  # Add '>' back when writing to file

print(f"Total sequences with no hits: {len(no_hits_ids)}")
