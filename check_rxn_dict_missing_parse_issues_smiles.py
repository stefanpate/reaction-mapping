from rdkit import Chem
import json
import csv

# path = 'metacyc_rxns_directed_221214.json'
# suffix = '_smiles.csv'

path = 'metacyc_directed_rxns_updated_smiles_221214.json'
suffix = '_smiles_post_update.csv'

# Load directed metacyc reactions
with open(path, 'r') as f:
    rxns = json.load(f)

# Break out rxns w/ cpds missing or triggering rdkit parse issues
cpds_missing_smiles = []
rxns_w_cpds_missing_smiles = []
cpds_parse_issues = []
rxns_w_cpds_parse_issues = []
for k,v in rxns.items(): # Over reactions
    for i in range(2): # Over reactant & product dicts
        for k2 in v[i].keys(): # Over cpds in reactant / product dicts
            if (v[i][k2] is None) and ([k2] not in cpds_missing_smiles): # Missing substrates
                cpds_missing_smiles.append([k2])

            if (v[i][k2] is None) and ([k] not in rxns_w_cpds_missing_smiles): # Missing reactions
                rxns_w_cpds_missing_smiles.append([k])
            
            if (v[i][k2] is not None): # Not missing
                temp_mol = Chem.MolFromSmiles(v[i][k2])
                if (temp_mol is None) and ([k2] not in cpds_parse_issues): # Parse issue substrates
                    cpds_parse_issues.append([k2])
                
                if (temp_mol is None) and ([k] not in rxns_w_cpds_parse_issues): # Parse issue reactions
                    rxns_w_cpds_parse_issues.append([k])

labels = ['cpds_missing', 'rxns_missing', 'cpds_parse', 'rxns_parse']
datas = cpds_missing_smiles, rxns_w_cpds_missing_smiles, cpds_parse_issues, rxns_w_cpds_parse_issues

for i in range(len(labels)):
    save_to = f"metacyc_{labels[i]}" + suffix
    with open(save_to, 'w') as f:
        writer = csv.writer(f)
        for elt in datas[i]:
            writer.writerow(elt)