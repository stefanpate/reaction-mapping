from rule_mapping_utils import MapRules
import pandas as pd
import json
import csv

# Designate paths to data
rxn_dict_path = 'rxn_dict_for_rule_mapper.json'
cpd_path = 'brenda_neutralize.tsv'
rules_path = 'minimal1224_all_uniprot.tsv'
cofactor_list_path = 'cofactor_list_seed_ids.csv'
cofactor_pair_path = 'cofactor_pairs_seed_ids.json'
SEED_neutralized = 'SEED_neutralized.tsv'
db_cpd_dict = {k: v['smiles'] for k, v in pd.read_csv(cpd_path, sep='\t', index_col=0).iterrows()}
save_to = 'mapped_reactions.csv'

# Initialize MapRule object
mapper = MapRules(rules_path=rules_path,
                    cofactor_list_path=cofactor_list_path,
                    cofactor_pair_path=cofactor_pair_path,
                    seed_dict=db_cpd_dict)

# # Test on erythronolide monoxygenase step
# # Reactants
# Deoxyerythronolide_SMILES = 'CCC1OC(=O)C(C)C(O)C(C)C(O)C(C)CC(C)C(=O)C(C)C(O)C1C'
# Deoxyerythronolide_SEED_ID = 'cpd02075' # '0:0'

# oxygen_smiles = 'O=O'
# oxygen_SEED_ID = 'cpd00007' # '1:0'

# NADH_smiles = 'NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1'
# NADH_SEED_ID = 'cpd00004' # '2:0' # 

# # Products
# erythronolide_SMILES = 'CCC1OC(=O)C(C)C(O)C(C)C(O)C(C)(O)CC(C)C(=O)C(C)C(O)C1C'
# erythronolide_SEED_ID = 'cpd04039:0' # '3:0'

# water_SMILES = 'O'
# water_SEED_ID = 'cpd00001' # '4:0' #

# NAD_SEED_ID = 'cpd00003' # '5:0' # 
# NAD_SMILES = 'NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1'

# rxn_dict_id = 'rxn127'
# rxn_dict = {rxn_dict_id:
#             [   # Reactants dictionary
#                 {Deoxyerythronolide_SEED_ID:Deoxyerythronolide_SMILES,
#                  NADH_SEED_ID:NADH_smiles,
#                  oxygen_SEED_ID:oxygen_smiles},

#                 # Products dictionary
#                 {erythronolide_SEED_ID:erythronolide_SMILES,
#                  NAD_SEED_ID:NAD_SMILES,
#                  water_SEED_ID:water_SMILES}]}

# Load in reaction dictionary
with open(rxn_dict_path, 'r') as f:
    rxn_dict = json.load(f)

# # Test with erythronolide monoxygenase
# rxn_dict_id = 'rxn041'
# test_rxn = rxn_dict[rxn_dict_id]
# rxn_dict = {rxn_dict_id:test_rxn}

# Test with this guy
test_rxn = rxn_dict['|RXN-12920|']
rxn_dict = {}
rxn_dict['|RXN-12920|'] = test_rxn

# Map rules to reactions
rxn_to_rule = []
for k,v in rxn_dict.items():
    for i in range(1,1225):
        rule = f"rule{i:04}"
        try:
            found_match = mapper.map_pickaxe_rules(v[0], v[1], rule)

            if found_match:
                print(f"{k} => {rule}")
                rxn_to_rule.append([k, rule])
        except ValueError:
            pass

print(f"{len(rxn_to_rule)} / {len(rxn_dict)} reactions mapped")

# Save rule mapping
# with open(save_to, 'w') as f:
#     writer = csv.writer(f)
#     for elt in rxn_to_rule:
#         writer.writerow(elt)