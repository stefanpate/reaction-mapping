import csv
import pythoncyc as pc
import json

cpd_id_to_pwy_fn = 'pk_cpd_id_pwy.csv'
save_to = 'rxn_dict_metacyc_ids.json'
starters = ['|PROPIONYL-COA|', '|MALONYL-COA|', '|ACETYL-COA|', '|METHYL-MALONYL-COA|', '|D-METHYL-MALONYL-COA|']
non_smiles = ['an oxidized electron carrier', 'a deaminated amino group donor'] # Catches metacyc substrates with non-smiles strings in their smiles attributes

meta = pc.so('meta') # Pull metacyc

# Read in polyketide pathways
pk_pwys = {}
with open(cpd_id_to_pwy_fn, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if len(row) == 1:
            pk_pwys[row[0]] = []
        else:
            pk_pwys[row[0]] = row[1:]

# Pull rxn ids from pwys
pwy_rxns = [] # Store pwy -> rxn_ids
rxn_id_list = [] # Store list of rxn ids
for cpd, pwy_list in pk_pwys.items():
    for pwy in pwy_list:
        rxns = meta[pwy].reaction_list
        pwy_rxns.append([cpd, pwy] + rxns)
        rxn_id_list = rxn_id_list + rxns

def substrate_list_to_dict(substrates, meta=meta, non_smiles=non_smiles):
    '''
    Given a list of substrates, returns dict of
    substrates with metacyc cpd ids as keys and 
    compound smiles as values. Filters out compounds
    with None or plain words as smiles entries.
    Returns:
        - substrates: dict of substrates or empty dict if 
        problem encountered
    '''

    sub_dict = {}
    for sub in substrates:
        smiles = meta[sub].smiles
        
        # # Catch None SMILES entries
        # if smiles is None:
        #     return {}

        # # Catch 'plain english' SMILES entries
        # for elt in non_smiles:
        #     if elt in smiles:
        #         return {}
    
        
        sub_dict[sub] = smiles

    return sub_dict

print(len(rxn_id_list))
rxn_id_list = list(set(rxn_id_list)) # Filter out repeat reactions
print(len(rxn_id_list))

# Filter out PKS reactions
temp = []
for rxn in rxn_id_list:
    substrates = meta[rxn].substrates
    has_starters = len(set(substrates) & set(starters)) > 0

    if not has_starters:
        temp.append(rxn)

print('Filter out PKS:', len(rxn_id_list), '->', len(temp))
rxn_id_list = sorted(temp)

# Get rxn dict for each rxn id
rxn_dict = {}
for rxn in rxn_id_list:
    reactants, products = meta[rxn].left, meta[rxn].right
    reactant_dict = substrate_list_to_dict(reactants)
    product_dict = substrate_list_to_dict(products)

    # If no issues in making substrate dicts, package into rxn dict
    # if (len(reactant_dict) > 0) & (len(product_dict) > 0):
    rxn_dict[rxn] = [reactant_dict, product_dict]

print(f"Filter out rxns w/out proper SMILES entries: {len(rxn_id_list)} -> {len(rxn_dict)}")

# Save reaction dict
with open(save_to, 'w') as f:
    json.dump(rxn_dict, f)
