import json
import csv

cofactor_path = 'cofactor_pairs_seed_ids.json'
save_path = 'cofactor_pair_counts_seed_id.csv'

# Read cofactor file
with open(cofactor_path, 'r') as f:
    cofactors = json.load(f)

# Count specific compound pairs per
# pair class / stand-in compounds
# (those used by the rules)
cofactor_counts = {}
for k, v in cofactors.items():
    cofactor_counts[k] = len(v)

# Save counts
with open(save_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(cofactor_counts.items())