import csv
import sys
import numpy as np
from mapping_fcns import *

# Run from cmd
# E.g., python map_rxns.py minimal1224_all_uniprot.tsv ../../bottle/mevalonate_pathway_neutralized.json ../../bottle/mapping_mevalonate_pathway_neutralized_min_rules.csv
rules_path = sys.argv[1]
rxn_dict_path = sys.argv[2]
save_to = sys.argv[3]
check_smiles = True # Whether or not to track missing smiles and parse_issues
path_list = save_to.split('/')
missing_smiles_path = '/'.join(path_list[:-1]) + '/' + 'missing_smiles_' + path_list[-1].lstrip('mapping_')
smiles_parse_path = '/'.join(path_list[:-1]) + '/' + 'smiles_parse_issues_' + path_list[-1].lstrip('mapping_')

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
rxns_missing_smiles = []
rxns_w_smiles_parse_issues = []
rxn_ctr = 0
mapped_rxn_binary = np.zeros(shape=(n_rxns,))
for k, rxn in rxn_dict.items():
    row = [k]
    rxn = apply_stoich(k, rxn, stoich_dict)
    for elt in rules:
        rule_name, rule_smarts = elt
        found_match, missing_smiles, smiles_parse_issue = map_rxn2rule(rxn, rule_smarts)

        if found_match:
            print(f"{k} => {rule_name}")
            row.append(rule_name)
            mapped_rxn_binary[rxn_ctr] = 1

        if missing_smiles & ([k] not in rxns_missing_smiles):
            rxns_missing_smiles.append([k])

        if smiles_parse_issue & ([k] not in rxns_w_smiles_parse_issues):
            rxns_w_smiles_parse_issues.append([k])
        
    rxn_to_rule.append(row)
    rxn_ctr += 1 # Update progress

    # Periodically print progress, save results and empty list
    if (rxn_ctr % 50 == 0) | (rxn_ctr == n_rxns):
        print(f"{mapped_rxn_binary.sum():.0f} / {rxn_ctr} reactions mapped") # Print progress

        # Save results
        with open(save_to, 'a') as f:
            writer = csv.writer(f)
            writer.writerows(rxn_to_rule)

        rxn_to_rule = [] # Empty list of mappings

        if check_smiles:

            # Save results
            with open(missing_smiles_path, 'a') as f:
                writer = csv.writer(f)
                writer.writerows(rxns_missing_smiles)

            # Save results
            with open(smiles_parse_path, 'a') as f:
                writer = csv.writer(f)
                writer.writerows(rxns_w_smiles_parse_issues)

            rxns_missing_smiles = []
            rxns_w_smiles_parse_issues = []

