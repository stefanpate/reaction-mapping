import json
import numpy as np
from itertools import permutations
from rdkit import Chem
from rdkit.Chem import AllChem, CanonSmiles

def map_rxn2rule(rxn, rule, max_products=10000):
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
            try:
                output = [CanonSmiles(Chem.MolToSmiles(elt)) for elt in output] # Convert pred products to smiles
            except:
                output = [Chem.MolToSmiles(elt) for elt in output]
            
            output = sorted(output)

            # Compare predicted to actual products. If mapped, return
            if output == products: 
                did_map = True
                return did_map, missing_smiles, smiles_parse_issue
    
    return did_map, missing_smiles, smiles_parse_issue

def apply_stoich(rxn_id, rxn, stoich_dict):
    '''
    Take reaction list-of-dicts and append incremental smiles
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

def save_json(data, save_to):
    with open(save_to, 'w') as f:
        json.dump(data, f)
    
def load_json(path):
    with open(path, 'r') as f:
        data = json.load(f)
    return data