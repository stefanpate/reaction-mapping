import csv
import json
import numpy as np
from itertools import combinations, permutations
from rdkit import Chem
from rdkit.Chem import AllChem

rules_path = 'minimal1224_all_uniprot.tsv'
save_to = 'mapped_reactions_simple.csv'
rxn_dict_path = 'rxn_dict_metacyc_ids.json'
rdkit_issues_path = 'rxns_with_rdkit_issues.csv'

def map_rxn2rule(rxn, rule):
    '''
    Maps reactions to SMARTS-encoded
    reaction rule
    Returns:
        - did_map, did_have_smiles_stereochem_issue
    '''
    reactants, products = list(rxn[0].values()), list(rxn[1].values())
    reactants, products = sanitize(reactants), sanitize(products)
    if (reactants is None) | (products is None):
        return False, True
    else:
        products = sorted(products)
        operator = Chem.rdChemReactions.ReactionFromSmarts(rule)
        reactants_mol = tuple([Chem.MolFromSmiles(elt) for elt in reactants])
        n_rule_reactants = count_reactants(rule)
        for reactant_subset in combinations(reactants_mol, n_rule_reactants):
            for perm in permutations(reactant_subset):
                outputs = operator.RunReactants(perm)
                for output in outputs:
                    output = [Chem.MolToSmiles(elt) for elt in output]
                    output = sorted(output)
                    if output == products:
                        return True, False
    return False, False

def count_reactants(rule_smarts):
    '''
    Counts number of reactants in a SMARTS-
    encoded operator
    '''
    reactants = rule_smarts.split('>>')[0]
    dot_split = reactants.split('.') # Reactants separated by '.'

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

    n_reactants = len(dot_split) - overcount
    return int(n_reactants)

def sanitize(list_of_smiles):
    sanitized_smiles = []
    for elt in list_of_smiles:
        try:
            temp_mol = Chem.MolFromSmiles(elt)
            Chem.rdmolops.RemoveStereochemistry(temp_mol)
            sanitized_smiles.append(Chem.MolToSmiles(temp_mol))
        except:
            try:
                temp_mol = Chem.MolFromSmiles(elt, sanitize=False)
                Chem.rdmolops.RemoveStereochemistry(temp_mol)
                sanitized_smiles.append(Chem.MolToSmiles(temp_mol))
            except:
                pass
    
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

# # Test with this guy
# test_rxn = rxn_dict['|RXN-12920|']
# rxn_dict = {}
# rxn_dict['|RXN-12920|'] = test_rxn

# Map reactions to rules
rxn_to_rule = []
rxns_rdkit_issue = []
mapped_rxns = []
for k,v in rxn_dict.items():
    row = [k]
    for elt in rules:
        rule_name, rule_smarts = elt
        found_match, rdkit_issue = map_rxn2rule(v, rule_smarts)

        if found_match:
            print(f"{k} => {rule_name}")
            row.append(rule_name)
            
            if k not in mapped_rxns:
                mapped_rxns.append(k)

        if (rdkit_issue) and ([k] not in rxns_rdkit_issue):
            rxns_rdkit_issue.append([k])
        
    rxn_to_rule.append(row)

print(f"{len(mapped_rxns)} / {len(rxn_dict)} reactions mapped")

# Save results
with open(save_to, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(rxn_to_rule)

with open(rdkit_issues_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(rxns_rdkit_issue)