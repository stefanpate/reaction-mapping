import csv
import sys
import numpy as np
from mapping_fcns import *

# Run from cmd
# E.g., python map_rxns_to_rules_simple.py minimal1224_all_uniprot.tsv 200_random_metacyc_rxns_seed_1234.json temp_mapped_rxns.csv temp_rdkit_issues_mapping.csv
# rules_path = sys.argv[1]
# rxn_dict_path = sys.argv[2]
# save_to = sys.argv[3]
# rdkit_issues_path = sys.argv[4]

# Set
rules_path = 'test_rule.csv'
save_to = 'test_mapping.csv'
rxn_dict_path = 'test_rxn.json'

# Set and forget
stoich_path = 'stoich_metacyc_rxns_directed_221214.json'

#################################################################################################
# Read in rules
rules = []
with open(rules_path, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        rules.append([row[0], row[2]])

rules = rules[1:] # Remove header

rxn_dict = load_json(rxn_dict_path) # Read in reactions
stoich_dict = load_json(stoich_path) # Read in stoichiometry

n_rxns = len(list(rxn_dict.keys())) # Total no. reactions to map

# Map reactions to rules
rxn_to_rule = []
rxn_ctr = 0
mapped_rxn_binary = np.zeros(shape=(n_rxns,))
for k, rxn in rxn_dict.items():
    row = [k]
    rxn = apply_stoich(k, rxn, stoich_dict)
    for elt in rules:
        rule_name, rule_smarts = elt
        found_match, _, _ = map_rxn2rule(rxn, rule_smarts)

        if found_match:
            print(f"{k} => {rule_name}")
            row.append(rule_name)
            mapped_rxn_binary[rxn_ctr] = 1
        
    rxn_to_rule.append(row)
    rxn_ctr += 1 # Update progress

    # Periodically print progress, save results and empty list
    if (rxn_ctr % 50 == 0) | (rxn_ctr == n_rxns):
        print(f"{mapped_rxn_binary.sum():.0f} / {rxn_ctr} reactions mapped") # Print progress

        # Save results
        with open(save_to, 'w') as f:
            writer = csv.writer(f)
            writer.writerows(rxn_to_rule)

        rxn_to_rule = [] # Empty list of mappings