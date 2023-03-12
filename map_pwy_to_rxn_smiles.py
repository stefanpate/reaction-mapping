import csv
import pythoncyc as pc
import json

cpd_id_to_pwy_fn = 'pk_cpd_id_pwy.csv'
metacyc_to_seed_id_fn = 'metacyc_to_seed_id.csv'
save_to = 'rxn_dict_for_rule_mapper.json'
rxn_dict_metacyc_path = 'rxn_dict_metacyc_ids.json'
subs_sans_seeds_path = 'substrates_without_model_seed_ids.csv'
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

metacyc_to_seed = {}
# Read in map of metacyc to seed ids
with open(metacyc_to_seed_id_fn, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if len(row) > 1:
            metacyc_to_seed[tuple(row[1:])] = row[0]

# Pull rxn ids from pwys
pwy_rxns = [] # Store pwy -> rxn_ids
rxn_id_list = [] # Store list of rxn ids
for cpd, pwy_list in pk_pwys.items():
    for pwy in pwy_list:
        rxns = meta[pwy].reaction_list
        pwy_rxns.append([cpd, pwy] + rxns)
        rxn_id_list = rxn_id_list + rxns

def check_subs_for_seeds(substrates, meta=meta, metacyc_to_seed=metacyc_to_seed, non_smiles=non_smiles):
    '''
    Given a list of substrates, returns dict of
    substrates of format required by Joseph's rule mapper
    '''
    sub_seeds = []
    subs_sans_seeds = []
    for sub in substrates:
        found_seed = False # Flag whether model SEED id found for substrate
        for k in metacyc_to_seed.keys():

            if sub.upper() in k:
                sub_seeds.append(metacyc_to_seed[k])
                found_seed = True
                break
        
        if not found_seed:
            subs_sans_seeds.append(sub)

    sub_dict = {}
    if len(substrates) == len(sub_seeds):
        
        for i in range(len(substrates)):
            smiles = meta[substrates[i]].smiles
            
            if smiles is None:
                return {}, subs_sans_seeds
            else:
                for elt in non_smiles:
                    if elt in smiles:
                        return {}, subs_sans_seeds
            
                sub_dict[sub_seeds[i]] = smiles

    return sub_dict, subs_sans_seeds

rxn_id_list = list(set(rxn_id_list)) # Filter out repeat reactions

# Filter out PKS reactions
temp = []
for rxn in rxn_id_list:
    substrates = meta[rxn].substrates
    has_starters = len(set(substrates) & set(starters)) > 0

    if not has_starters:
        temp.append(rxn)

print('Filter out PKS:', len(rxn_id_list), '->', len(temp))
rxn_id_list = sorted(temp)

# Convert rxn ids to smiles
rxn_dict = {}
rxn_dict_to_metacyc_id = []
subs_sans_seeds = []
rxn_ct = 0
for rxn in rxn_id_list:
    reactants, products = meta[rxn].left, meta[rxn].right
    reactant_dict, reactants_sans_seeds = check_subs_for_seeds(reactants)
    product_dict, products_sans_seeds = check_subs_for_seeds(products)
    subs_sans_seeds += reactants_sans_seeds + products_sans_seeds

    if (len(reactant_dict) > 0) & (len(product_dict) > 0):
        key = f"rxn{rxn_ct:03}"
        rxn_dict[rxn] = [reactant_dict, product_dict]
        rxn_ct += 1

subs_sans_seeds = sorted(list(set(subs_sans_seeds))) # Remove duplicates

print(f"Filter out rxns w/out SEED id compounds: {len(rxn_id_list)} -> {len(rxn_dict)}")

# Save reaction dict
with open(save_to, 'w') as f:
    json.dump(rxn_dict, f)

# Save substrates without model SEED ids
with open(subs_sans_seeds_path, 'w') as f:
    writer = csv.writer(f)
    for elt in subs_sans_seeds:
        writer.writerow([elt])
