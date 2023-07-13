from rdkit import Chem
import json
import csv

path = 'metacyc_rxns_directed_221214.json'

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
            if (v[i][k2] is None) and ([k2] not in cpds_missing_smiles): # Missing
                cpds_missing_smiles.append([k2])
            elif (v[i][k2] is None) and ([k] not in rxns_w_cpds_missing_smiles): # Missing
                rxns_w_cpds_missing_smiles.append([k])
            elif (v[i][k2] is not None): # Not missing
                temp_mol = Chem.MolFromSmiles(v[i][k2])
                if (temp_mol is None) and ([k2] not in cpds_parse_issues): # Parse issue
                    cpds_parse_issues.append([k2])
                elif (temp_mol is None) and ([k] not in rxns_w_cpds_parse_issues):
                    rxns_w_cpds_parse_issues.append([k])

labels = ['cpds_missing', 'rxns_missing', 'cpds_parse', 'rxns_parse']
datas = cpds_missing_smiles, rxns_w_cpds_missing_smiles, cpds_parse_issues, rxns_w_cpds_parse_issues

for i in range(len(labels)):
    save_to = f"metacyc_{labels[i]}_smiles_post_update.csv"
    with open(save_to, 'w') as f:
        writer = csv.writer(f)
        for elt in datas[i]:
            writer.writerow(elt)