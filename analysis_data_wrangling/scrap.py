from rdkit import Chem
from rdkit.Chem import AllChem
import json

def sanitize(list_of_smiles):
    sanitized_smiles = []
    for elt in list_of_smiles:
        temp_mol = Chem.MolFromSmiles(elt)
        Chem.rdmolops.RemoveStereochemistry(temp_mol)
        sanitized_smiles.append(Chem.MolToSmiles(temp_mol))    
    return sanitized_smiles

# Load final metacyc rxns
mc_rxns_path = 'mc_rxns_final_221214.json'
with open(mc_rxns_path, 'r') as f:
    mc_rxns = json.load(f)

stereoisomerizations = {}
problems = []
for k,v in mc_rxns.items():
    
    destereoed = []
    for i in range(2):
        try:
            destereoed.append(sanitize(v[i].values()))
            intrsct = list(set(destereoed[0]) & set(destereoed[1]))
            if len(intrsct) > 0:
                stereoisomerizations[k] = intrsct
        
        except:
            problems.append(k)

problems = list(set(problems))
print(len(stereoisomerizations.keys()), len(problems))
