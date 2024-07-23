# -*- coding: utf-8 -*-
"""
@author: Shweta Nagre
"""

import pandas as pd

# Define file paths
mmseqs_human_hits_result_file = 'Human_Hits_for_QUERY.tab'
mmseqs_nonhuman_hits_result_file = 'NonHuman_Hits_for_QUERY.tab'
whole_region_file = 'whole_region.txt'
no_human_hits_file = 'no_human_hits_ids.txt'
no_nonhuman_hits_file = 'no_nonhuman_hits_ids.txt'
final_no_hits_file = 'final_no_hits_ids.txt'

# Read the MMseqs result file to get the IDs with hits
## First for humans
human_hits_ids = pd.read_csv(mmseqs_human_hits_result_file, sep='\t', header=None)
human_hits_ids_set = set(human_hits_ids[0].tolist())

## Now, for non-humans
nonhuman_hits_ids = pd.read_csv(mmseqs_nonhuman_hits_result_file, sep='\t', header=None)
nonhuman_hits_ids_set = set(nonhuman_hits_ids[0].tolist())

# Read the whole_region.txt file to get all query IDs
query_ids = []
with open(whole_region_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            query_id = line.strip()[1:]  # Remove '>' for matching purposes
            query_ids.append(query_id)

# Identify the IDs that did not get any hits in both human and non-human databases
no_human_hits_ids = [query_id for query_id in query_ids if query_id not in human_hits_ids_set]

no_nonhuman_hits_ids = [query_id for query_id in query_ids if query_id not in nonhuman_hits_ids_set]

# Write the IDs with no human hits to a file
with open(no_human_hits_file, 'w') as f:
    for no_human_hit_id in no_human_hits_ids:
        f.write('>' + no_human_hit_id + '\n')  # Add '>' back when writing to file

print(f"Total sequences with no human hits: {len(no_human_hits_ids)}")

# Write the IDs with no non-human hits to a file
with open(no_nonhuman_hits_file, 'w') as f:
    for no_nonhuman_hits_id in no_nonhuman_hits_ids:
        f.write('>' + no_nonhuman_hits_id + '\n')  # Add '>' back when writing to file

print(f"Total sequences with no non-human hits: {len(no_nonhuman_hits_ids)}")

# Concatenate the two files into a single file
with open(final_no_hits_file, 'w') as outfile:
    with open(no_human_hits_file, 'r') as infile:
        outfile.write(infile.read())
    with open(no_nonhuman_hits_file, 'r') as infile:
        outfile.write(infile.read())

print(f"Concatenated file created: {final_no_hits_file}")
