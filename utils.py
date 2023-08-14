from rdkit import Chem
from rdkit.Chem import AllChem
import json
from collections import Counter

def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

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

def count_atoms(smi_str):
    '''
    Count atoms in provided smile string and
    return dict of element:count pairs.
    '''
    # Single letter elt dict w lowercase rpts for aromatics
    one_letter_elts = {'B':0, 'C':0, 'N':0, 'O':0,
                        'F':0, 'P':0, 'S':0, 'I':0,
                        'W':0, 'c':0, 'n':0, 'o':0,
                        's':0}
    
    aromatic_elts = ['c', 'n', 'o', 's'] # Keep track of aromatic elts
    
    # Two letter elements
    two_letter_elts = {'Li':0, 'Be':0, 'Ne':0, 'Na':0, 'Mg':0,
                        'Al':0, 'Cl':0, 'Ca':0, 'Fe':0,'Ni':0,
                        'Zn':0, 'Br':0, 'Te':0, 'As':0, 'Sb':0,
                        'Cu':0}

    # Count single letter elements...
    # use Counter for single letter elts
    single_elt_cts = Counter(one_letter_elts)
    single_elt_cts.update(smi_str)
    temp = {} # To replace one_letter_elts
    for k in one_letter_elts.keys():
        temp[k] = single_elt_cts[k]

    # Combine aromatic counts
    for elt in aromatic_elts:
        temp[elt.upper()] = temp[elt.upper()] + temp[elt]
        temp.pop(elt) # Remove lowercase keys after done w em

    # Count two letter elements...
    # use .count() for two letter elts
    for k in two_letter_elts.keys():
        two_letter_elts[k] = smi_str.count(k)

    cts = {**temp, **two_letter_elts} # Combine one, two letter cts
    return cts

def is_balanced(rxn):
    '''
    Returns whether a reaction is balanced.
    Args:
        - rxn: list of dicts w/ [{react1:smi1, }, {prod1:smi1, }]
    '''
    side_cts = [count_atoms(''), count_atoms('')]
    for i, side in enumerate(rxn):
        for cpd in side.values():
            cpd_cts = count_atoms(cpd)
            for k in side_cts[i].keys():
                side_cts[i][k] += cpd_cts[k]

    diff_vec = [side_cts[0][k] - side_cts[1][k] for k in side_cts[0].keys()]
    balanced = all([elt == 0 for elt in diff_vec])
    return balanced