import csv
import json
import sys
import numpy as np
from itertools import combinations, permutations
from rdkit import Chem
from rdkit.Chem import AllChem

# Run from cmd
# E.g., python map_rxns_to_rules_simple.py minimal1224_all_uniprot.tsv 200_random_metacyc_rxns_seed_1234.json temp_mapped_rxns.csv temp_rdkit_issues_mapping.csv
# rules_path = sys.argv[1]
# rxn_dict_path = sys.argv[2]
# save_to = sys.argv[3]
# rdkit_issues_path = sys.argv[4]

# Run from editor
rules_path = 'minimal1224_all_uniprot.tsv'
save_to = 'test_rxn_mapping.csv'
rxn_dict_path = '100_random_metacyc_rxns_rnd_seed_1234.json'
missing_smiles_path = 'filtered_missing_smiles.csv'
parse_issues_path = 'filtered_smiles_parse_issues.csv'
stoich_path = 'stoich_metacyc_rxns_directed_221214.json'

def map_rxn2rule(rxn, rule, max_products=1000):
    '''
    Maps reactions to SMARTS-encoded reaction rule.
    Args:
        - rxn: List of lists, each sublist with smiles of 
        substrates with the correct multiplicity / stoichiometry
        - rule: smarts string
    Returns:
        - did_map (bool)
        - missing_smiles (bool)
        - smiles_parse_issue (bool)
    '''
    pre_sani_reactants, pre_sani_products = rxn # Get lists of smiles

    # Initialize flags
    did_map = False
    missing_smiles = False
    smiles_parse_issue = False

    # Check for missing smiles
    if (None in pre_sani_reactants) or (None in pre_sani_products):
        missing_smiles = True
    
    reactants, products = sanitize(pre_sani_reactants), sanitize(pre_sani_products) # Remove stereochem
    
    if ((len(reactants) != len(pre_sani_reactants)) or (len(products) != len(pre_sani_products))) and (not missing_smiles) :
        smiles_parse_issue = True
    
    products = sorted(products)
    operator = Chem.rdChemReactions.ReactionFromSmarts(rule) # Make reaction object from smarts string
    reactants_mol = tuple([Chem.MolFromSmiles(elt) for elt in reactants]) # Convert reactant smiles to mol obj
    rule_substrate_cts = count_reactants(rule) # [n_reactants, n_products] in a rule
    rxn_substrate_cts = [len(reactants), len(products)]

    # Check if number of reactants / products strictly match
    # rule to reaction
    if rule_substrate_cts != rxn_substrate_cts:
        return did_map, missing_smiles, smiles_parse_issue
        
    # For every permutation of that subset of reactants
    for perm in permutations(reactants_mol):
        outputs = operator.RunReactants(perm, maxProducts=max_products) # Apply rule to that permutation of reactants
        for output in outputs:
            output = [Chem.MolToSmiles(elt) for elt in output] # Convert pred products to smiles
            output = sorted(output)

            # Compare predicted to actual products. If mapped, return
            if output == products: 
                did_map = True
                return did_map, missing_smiles, smiles_parse_issue
    
    return did_map, missing_smiles, smiles_parse_issue

def apply_stoich(rxn_id, rxn, stoich_dict):
    '''
    Take reaction dicts and append incremental smiles
    according to stoichiometry.
    Args:
        - rxn_id: string metacyc identifier
        - rxn: list of two dicts, with id:smi entries for substrates
    Returns:
        - reactants_smi, products_smi: two lists of substrate smiles w/ right multiplicity
    '''
    if rxn_id in stoich_dict.keys():
        output = []
        for i,elt in enumerate(rxn): # Each side of reaction
            temp = []
            for k,v in elt.items(): # Each id:smi
                coeff = stoich_dict[rxn_id][i][k]
                single_substrate = [v for i in range(coeff)]
                temp += single_substrate

            output.append(temp)

    else: # Can't find stoich
        output = [list(rxn[0].values()), list(rxn[1].values())]

    return output


def count_reactants(rule_smarts):
    '''
    Counts number of reactants in a SMARTS-
    encoded operator
    '''
    sides = rule_smarts.split('>>')
    cts = []
    for side in sides:
        dot_split = side.split('.') # Reactants separated by '.'

        # But must catch where pieces of a single compound
        # are split by '.', in which case they'll be surrounded by ()
        left_split_parens = []
        right_split_parens = []
        for i, elt in enumerate(dot_split):
            if (elt[0] == '(') & (elt[-1] != ')'):
                left_split_parens.append(i)
            elif (elt[0] != '(') & (elt[-1] == ')'):
                right_split_parens.append(i)
                
        left_split_parens, right_split_parens = np.array(left_split_parens), np.array(right_split_parens)
        overcount = (right_split_parens - left_split_parens).sum()

        n = len(dot_split) - overcount
        cts.append(int(n))
    return cts

def sanitize(list_of_smiles):
    sanitized_smiles = []
    for elt in list_of_smiles:
        temp_mol = Chem.MolFromSmiles(elt)
        Chem.rdmolops.RemoveStereochemistry(temp_mol)
        sanitized_smiles.append(Chem.MolToSmiles(temp_mol))    
    return sanitized_smiles

# Read in rules
rules = []
with open(rules_path, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        rules.append([row[0], row[2]])

rules = rules[1:] # Remove header

# Read in reactions
with open(rxn_dict_path, 'r') as f:
    rxn_dict = json.load(f)

# Read in stoichiometry
with open(stoich_path, 'r') as f:
    stoich_dict = json.load(f)

n_rxns = len(list(rxn_dict.keys()))

# Map reactions to rules
rxn_to_rule = []
rxns_parse_issue = []
rxns_missing_smiles = []
rxn_ctr = 0
mapped_rxn_binary = np.zeros(shape=(n_rxns,))
for k,v in rxn_dict.items():
    row = [k]
    rxn = apply_stoich(k, v, stoich_dict)
    # print(k)
    for elt in rules:
        rule_name, rule_smarts = elt
        # print(rule_name)
        found_match, missing_smiles, parse_issue = map_rxn2rule(rxn, rule_smarts)

        if found_match:
            print(f"{k} => {rule_name}")
            row.append(rule_name)
            mapped_rxn_binary[rxn_ctr] = 1

        if (parse_issue) and ([k] not in rxns_parse_issue):
            rxns_parse_issue.append([k])

        if (missing_smiles) and ([k] not in rxns_missing_smiles):
            rxns_missing_smiles.append([k])
        
    rxn_to_rule.append(row)

    # Print progress
    if rxn_ctr % 20 == 0:
        print(f"{mapped_rxn_binary.sum()} / {rxn_ctr} reactions mapped")

    rxn_ctr += 1 # Update progress

print(f"{mapped_rxn_binary.sum()} / {n_rxns} reactions mapped") # Final result

# Save results
with open(save_to, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(rxn_to_rule)

with open(parse_issues_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(rxns_parse_issue)

with open(missing_smiles_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(rxns_missing_smiles)