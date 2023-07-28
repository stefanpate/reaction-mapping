import csv
import json
import sys
import numpy as np
from itertools import combinations, permutations
from rdkit import Chem
from rdkit.Chem import AllChem, CanonSmiles

# Run from cmd
# E.g., python map_rxns_to_rules_simple.py minimal1224_all_uniprot.tsv 200_random_metacyc_rxns_seed_1234.json temp_mapped_rxns.csv temp_rdkit_issues_mapping.csv
# rules_path = sys.argv[1]
# rxn_dict_path = sys.argv[2]
# save_to = sys.argv[3]
# rdkit_issues_path = sys.argv[4]

# Set
rules_path = 'test_rule.csv'
save_to = 'mapping_w_paired_cofactor_replacement.csv'
rxn_dict_path = 'test_rxn.json'
replace_unpaired_cofactors = False
replace_paired_cofactors = True

# Set and forget
paired_cofactor_path = 'mc_paired_cofactors_lut.json'
unpaired_cofactor_path = 'mc_unpaired_cofactors_lut.json'
cofactor_to_smi_path = 'cofactor_class_to_smi.json'
stoich_path = 'stoich_metacyc_rxns_directed_221214.json'

def map_rxn2rule(rxn, rule, max_products=5000):
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

def substitute_cofactor_smiles(rxn, class2smi, paired_cofactors, unpaired_cofactors):
    '''
    Replaces the smiles values in the rxn dicts with 
    strings from cofactor lookup tables.
    Args:
        - rxn: List of reaction dicts [{'rid':'smi', }, {'pid':'smi', }]
        - paired_cofactors: Dict {'class_pair1,class_pair1':[['id1', 'id2'],],}
        - unpaired_cofactors: {'class_name':['id1',]}
        - cofactor2smi: LUT from cofactor class key to smiles {'class_name':'smi',}
    Returns:
        - rxn: W/ updated smiles
    '''
    rxn_ids_only = [[r for r in rxn[0].keys()], [p for p in rxn[1].keys()]] # To keep track of replacements
    
    # Defensive copy of rxn list-of-dicts...
    # passing the value of a dict actually
    # passes read & write permission, not 
    # just a view
    temp = [elt.copy() for elt in rxn]
    rxn = temp

    # Replace smiles of unpaired cofactors
    if unpaired_cofactors is not None:
        for i, side in enumerate(rxn):
            for id in side.keys():
                for k, v in unpaired_cofactors.items():
                    if id in v:
                        side[id] = class2smi[k]
                        rxn_ids_only[i].remove(id) # Remove ids identified as single cofactors

    # Replace smiles of paired cofactors
    if paired_cofactors is not None:
        for r in rxn_ids_only[0]:
            for k, v in paired_cofactors.items(): # Over 'name1_name2':list of pairs, pairs
                    for elt in v: # Each pair of ids
                        if r in elt:
                            for p in rxn_ids_only[1]:
                                if p in elt:
                                    r_idx, p_idx = elt.index(r), elt.index(p)
                                    r_class, p_class = k.split(',')[r_idx], k.split(',')[p_idx]
                                    r_smi, p_smi = class2smi[r_class], class2smi[p_class]
                                    rxn[0][r], rxn[1][p] = r_smi, p_smi

    return rxn

def save_json(data, save_to):
    with open(save_to, 'w') as f:
        json.dump(data, f)
    
def load_json(path):
    with open(path, 'r') as f:
        data = json.load(f)
    return data

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

# Read in cofactors
cofactor_class_to_smi = load_json(cofactor_to_smi_path)

if replace_unpaired_cofactors:
    unpaired_cofactors = load_json(unpaired_cofactor_path)
else:
    unpaired_cofactors = None

if replace_paired_cofactors:
    paired_cofactors = load_json(paired_cofactor_path)
else:
    paired_cofactors = None

n_rxns = len(list(rxn_dict.keys()))

# Map reactions to rules
rxn_to_rule = []
rxns_parse_issue = []
rxns_missing_smiles = []
rxn_ctr = 0
mapped_rxn_binary = np.zeros(shape=(n_rxns,))
for k, rxn in rxn_dict.items():
    row = [k]
    rxn = substitute_cofactor_smiles(rxn, cofactor_class_to_smi, paired_cofactors, unpaired_cofactors)
    rxn = apply_stoich(k, rxn, stoich_dict)
    for elt in rules:
        rule_name, rule_smarts = elt
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

# # Save results
# with open(save_to, 'w') as f:
#     writer = csv.writer(f)
#     writer.writerows(rxn_to_rule)