# -*- coding: utf-8 -*-
"""
@author: Shweta Nagre
"""

# Define file paths
whole_region_file = 'whole_region.txt'
no_hits_file = 'final_no_hits_ids.txt'
final_output_file = 'Final_MIPS.txt'

# Read no hits IDs into a set for fast lookup
with open(no_hits_file, 'r') as f:
    no_hits_ids = {line.strip() for line in f}

# Prepare to collect sequences
no_hits_sequences = []
current_id = None
current_sequence = []

# Read through whole_region.txt and collect sequences for no hits IDs
with open(whole_region_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            # Process the previous sequence if it was a no hit
            if current_id in no_hits_ids:
                no_hits_sequences.append((current_id, ''.join(current_sequence)))
            
            # Start collecting new sequence
            current_id = line.strip()
            current_sequence = []
        else:
            current_sequence.append(line.strip())

# Process the last sequence in the file
if current_id in no_hits_ids:
    no_hits_sequences.append((current_id, ''.join(current_sequence)))

# Write the no hits IDs and their sequences to a new file
with open(final_output_file, 'w') as f:
    for id, sequence in no_hits_sequences:
        f.write(f"{id}\n{''.join(sequence)}\n")

print(f"Total sequences with no hits: {len(no_hits_sequences)}")
print(f"Final output written to {final_output_file}")
