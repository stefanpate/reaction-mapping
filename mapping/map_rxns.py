import sys
import numpy as np
from mapping_fcns import *
import pandas as pd
import json

# Run from cmd
# E.g., python map_rxns.py JN3604IMT_rules.tsv swissprot_unmapped.json swissprot_unmapped.json
rules_path = sys.argv[1]
rxn_dict_path = sys.argv[2]
save_to = sys.argv[3]
check_smiles = False # Whether to save info abt missing smiles & parse_issues
do_template = True # Whether to enforce template matching, ie cofactors
return_rc = True # Whether to return reaction center while mapping operators
path_list = save_to.split('/')
smiles_parse_path = 'smiles_parse_issues_' + path_list[-1].lstrip('mapping')[:-5] + ".csv"

#################################################################################################
# Read in rules
rules = pd.read_csv(rules_path, sep='\t')
rules.set_index("Name", inplace=True)
rules.drop('Comments', axis=1, inplace=True)

rxn_dict = load_json(rxn_dict_path) # Read in reactions
n_rxns = len(list(rxn_dict.keys())) # Total no. reactions to map

# Read in cofactor lookup tables
if rules_path == 'minimal1224_all_uniprot.tsv':
    paired_cof_path = 'smi2paired_cof_min.json'
    unpaired_cof_path = 'smi2unpaired_cof_min.json'
elif rules_path == 'JN3604IMT_rules.tsv':
    paired_cof_path = 'smi2paired_cof_imt.json'
    unpaired_cof_path = 'smi2unpaired_cof_imt.json'
else:
    paired_cof_path = 'smi2paired_cof.json'
    unpaired_cof_path = 'smi2unpaired_cof.json'

if do_template:
    smi2paired_cof = load_json(paired_cof_path)
    smi2unpaired_cof = load_json(unpaired_cof_path)
    smi2paired_cof = {tuple(k.split(',')):v[0].split(',') for k,v in smi2paired_cof.items()}

processed_reactions = []
mapped_rules = []
reaction_centers = []
rxns_w_smiles_parse_issues = []
for z, (k, rxn) in enumerate(rxn_dict.items()): # Iterate over rxns
    smiles_parse_issue = False # Init data checking flag
    this_mapped_rules = [] # To store mapped ops for this rxn
    this_rcs = [] # Store reaction centers for this rxn

    # Try to sanitize and remove stereochem
    reactants, reactants_parse_issue = sanitize(list(rxn[0].values()))
    products, products_parse_issue = sanitize(list(rxn[1].values()))
    rxn = [reactants, products] # Sanitized reaction substrates

    # Catch smiles parse issues
    if reactants_parse_issue or products_parse_issue:
        smiles_parse_issue = True
        rxns_w_smiles_parse_issues.append(k)
    
    else: # If the data is okay, try to map
        for rule_name, row in rules.iterrows(): # Iterate over rules
            did_map = False # Init mapping flag
            rule_reactants_template, rule_smarts, rule_products_template = row

            if do_template:
                matched_idxs = match_template(rxn, rule_reactants_template, rule_products_template, smi2paired_cof, smi2unpaired_cof)

                if len(matched_idxs) > 0:
                    did_map, rc = map_rxn2rule(rxn, rule_smarts, return_rc=return_rc, matched_idxs=matched_idxs) # Map if have template matches

            else:
                did_map, rc = map_rxn2rule(rxn, rule_smarts, return_rc=return_rc) # Map w/o enforcing templates

            if did_map:
                print(f"{k} => {rule_name}")
                this_mapped_rules.append(rule_name)
                this_rcs.append(rc)

    processed_reactions.append(k)               
    mapped_rules.append(this_mapped_rules)
    reaction_centers.append(this_rcs)

    # Periodically print progress, save results and empty list
    if (z % 50 == 0) | (z == n_rxns):
        print(f"{z} reactions processed") # Print progress

# Save results
all_equal_length = len(processed_reactions) == len(mapped_rules) and len(mapped_rules) == len(reaction_centers)
if not all_equal_length:
    raise Exception("Reactions, rules, rcs lists not of equal length")

res = {
    processed_reactions[i]: {'rules': mapped_rules[i], 'rcs': reaction_centers[i]}
    for i in range(len(processed_reactions))
}

with open(save_to, 'w') as f:
    json.dump(res, f)

parse = pd.DataFrame({'Reactions': rxns_w_smiles_parse_issues})
parse.to_csv(smiles_parse_path, sep='\t', index=False)
