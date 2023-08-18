import csv
import sys
import numpy as np
from mapping_fcns import *

# Run from cmd
# E.g., python map_rxns.py minimal1224_all_uniprot.tsv 10_true_pos_mc_v21_rxns_rnd_seed_1234.json mapping_test_3.csv
rules_path = sys.argv[1]
rxn_dict_path = sys.argv[2]
save_to = sys.argv[3]
check_smiles = True # Whether to save info abt missing smiles & parse_issues
do_template = True # Whether to enforce template matching, ie cofactors
path_list = save_to.split('/')
missing_smiles_path = 'missing_smiles' + path_list[-1].lstrip('mapping')
smiles_parse_path = 'smiles_parse_issues' + path_list[-1].lstrip('mapping')

# Set and forget
stoich_path = 'stoich_metacyc_rxns_directed_221214.json'
paired_cof_path = 'smi2paired_cof.json'
unpaired_cof_path = 'smi2unpaired_cof.json'

#################################################################################################
# Read in rules
rules = []
with open(rules_path, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        rules.append(row[:4])

rules = rules[1:] # Remove header
rxn_dict = load_json(rxn_dict_path) # Read in reactions
stoich_dict = load_json(stoich_path) # Read in stoichiometry
n_rxns = len(list(rxn_dict.keys())) # Total no. reactions to map

# Read in cofactor lookup tables
if do_template:
    smi2paired_cof = load_json(paired_cof_path)
    smi2unpaired_cof = load_json(unpaired_cof_path)
    smi2paired_cof = {tuple(k.split(',')):v[0].split(',') for k,v in smi2paired_cof.items()}

rxn_to_rule = []
rxns_missing_smiles = []
rxns_w_smiles_parse_issues = []
rxn_ctr = 0
mapped_rxn_binary = np.zeros(shape=(n_rxns,))
for k, rxn in rxn_dict.items(): # Iterate over rxns
    # Init data checking flags
    missing_smiles = False
    smiles_parse_issue = False
    row = [k]
    rxn = apply_stoich(k, rxn, stoich_dict) # Returns list of lists of smiles

    # Check that the data is okay
    # Check for missing smiles
    if (None in rxn[0]) or (None in rxn[1]):
        missing_smiles = True
    else:
        # Try to sanitize and remove stereochem
        reactants, reactants_parse_issue = sanitize(rxn[0])
        products, products_parse_issue = sanitize(rxn[1])
        rxn = [reactants, products] # !!
    
        # Catch smiles parse issues
        if reactants_parse_issue or products_parse_issue:
            smiles_parse_issue = True
    
    # If the data is okay, try to map
    if (not missing_smiles) & (not smiles_parse_issue):
        for elt in rules: # Iterate over rules
            did_map = False # Init mapping flag
            rule_name, rule_reactants_template, rule_smarts, rule_products_template = elt

            if do_template:
                matched_idxs = match_template(rxn, rule_reactants_template, rule_products_template, smi2paired_cof, smi2unpaired_cof)

                if len(matched_idxs) > 0:
                    did_map = map_rxn2rule(rxn, rule_smarts, matched_idxs=matched_idxs) # Map if have template matches

            else:
                did_map = map_rxn2rule(rxn, rule_smarts) # Map w/o enforcing templates

            if did_map:
                print(f"{k} => {rule_name}")
                row.append(rule_name)
                mapped_rxn_binary[rxn_ctr] = 1

            if missing_smiles & ([k] not in rxns_missing_smiles):
                rxns_missing_smiles.append([k])

            if smiles_parse_issue & ([k] not in rxns_w_smiles_parse_issues):
                rxns_w_smiles_parse_issues.append([k])
            
    rxn_to_rule.append(row) # Only count parseable smiles to mapping list
    rxn_ctr += 1 # Update progress in any case: missing, unparseable, map, did not map rxn

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
