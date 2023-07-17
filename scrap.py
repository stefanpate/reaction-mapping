from itertools import combinations, permutations
from rdkit import Chem
import json
from rdkit.Chem import AllChem
import numpy as np

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

    n_start, n_sanitized = len(list_of_smiles), len(sanitized_smiles)
    # if n_start != n_sanitized:
    #     print('smiles issue')
    
    return sanitized_smiles

def map_rxn2rule(rxn, rule):
    '''
    Maps reactions to SMARTS-encoded reaction rule.
    Args:
        - rxn: 
        - rule: smarts string
    Returns:
        - did_map (bool)
        - missing_smiles (bool)
        - smiles_parse_issue (bool)
    '''
    pre_sani_reactants, pre_sani_products = list(rxn[0].values()), list(rxn[1].values()) # Get lists of smiles

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
    n_rule_reactants = count_reactants(rule) # No. reactants in rule

    # For every combo of reaction reactants
    for reactant_subset in combinations(reactants_mol, n_rule_reactants):
        
        # For every permutation of that subset of reactants
        for perm in permutations(reactant_subset):
            outputs = operator.RunReactants(perm, maxProducts=5000) # Apply rule to that permutation of reactants
            for output in outputs:
                output = [Chem.MolToSmiles(elt) for elt in output] # Convert pred products to smiles
                output = sorted(output)

                # Compare predicted to actual products. If mapped, return
                if output == products: 
                    did_map = True
                    return did_map, missing_smiles, smiles_parse_issue
    
    return did_map, missing_smiles, smiles_parse_issue



# Read in reactions
rxn_dict_path = '10_random_metacyc_rxns_filtered_rnd_seed_1234.json'
with open(rxn_dict_path, 'r') as f:
    rxn_dict = json.load(f)

rxn = rxn_dict['|RXN-10132|']
rule_id = 'rule0604'
rule = '([#8:1].[#8:2]-[#15:3]-[#8:4]).[#8:5]>>([#8:2].[#8:1]-[#15:3]-[#8:5]).[#8:4]'
map_rxn2rule(rxn, rule)