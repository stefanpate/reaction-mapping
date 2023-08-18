import json
import numpy as np
from itertools import permutations, product
from rdkit import Chem
from rdkit.Chem import AllChem, CanonSmiles

def map_rxn2rule(rxn, rule, matched_idxs=None, max_products=10000):
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
    reactants, products = rxn
    
    # If neither missing nor unparseable smiles, continue with mapping
    products = sorted(products)
    operator = Chem.rdChemReactions.ReactionFromSmarts(rule) # Make reaction object from smarts string
    reactants_mol = [Chem.MolFromSmiles(elt) for elt in reactants] # Convert reactant smiles to mol obj
    rule_substrate_cts = count_reactants(rule) # [n_reactants, n_products] in a rule
    rxn_substrate_cts = [len(reactants), len(products)]

    # Check if number of reactants / products strictly match
    # rule to reaction
    if rule_substrate_cts != rxn_substrate_cts:
        return False
    
    # If not enforcing templates,
    # get all permutations of reactant
    # indices
    if matched_idxs is None:
        matched_idxs = list(permutations([i for i in range(len(reactants))]))
        
    # For every permutation of that subset of reactants
    for idx_perm in matched_idxs:
        perm = tuple([reactants_mol[idx] for idx in idx_perm]) # Re-order reactants based on allowable idx perms
        outputs = operator.RunReactants(perm, maxProducts=max_products) # Apply rule to that permutation of reactants

        for output in outputs:
            try:
                output = [CanonSmiles(Chem.MolToSmiles(elt)) for elt in output] # Convert pred products to canonical smiles
            except:
                output = [Chem.MolToSmiles(elt) for elt in output]
            
            output = sorted(output)

            # Compare predicted to actual products. If mapped, return
            if output == products: 
                return True
            
            # Last, try fixing kekulization issues
            postsan_output = postsanitize_smiles(output)
            for elt in postsan_output: # Iterate over sets of outputs w/ diff tautomers
                if sorted(elt) == products:
                    return True
    return False

def match_template(rxn, rule_reactants_template, rule_products_template, smi2paired_cof, smi2unpaired_cof):
    '''
    Returns the permuted indices corresponding to
    a match between reactant and rule templates
    '''
    reactants_smi, products_smi = rxn
    rule_reactants_template = tuple(rule_reactants_template.split(';'))
    rule_products_template = tuple(rule_products_template.split(';'))
    matched_idxs = [] # Return empty if no matches found
    # First check the cardinality of reactants, products matches
    if (len(rule_reactants_template) == len(reactants_smi)) & (len(rule_products_template) == len(products_smi)):

        reactants_template = ['Any' for elt in reactants_smi]
        products_template = ['Any' for elt in products_smi]

        # Search for unpaired cofactors first
        for i, r in enumerate(reactants_smi):
            if r in smi2unpaired_cof:
                reactants_template[i] = smi2unpaired_cof[r]

        for i, p in enumerate(products_smi):
            if p in smi2unpaired_cof:
                products_template[i] = smi2unpaired_cof[p]

        # Search for paired cofactors
        # Only overwriting should be PPi/Pi as phosphate donor/acceptor
        for i, r in enumerate(reactants_smi):
            for j, p in enumerate(products_smi):
                if (r, p) in smi2paired_cof:
                    reactants_template[i] = smi2paired_cof[(r, p)][0]
                    products_template[j] = smi2paired_cof[(r, p)][1]
                elif (p, r) in smi2paired_cof:
                    reactants_template[i] = smi2paired_cof[(p, r)][1]
                    products_template[j] = smi2paired_cof[(p, r)][0]

        reactants_idx_template = [(elt, i) for i, elt in enumerate(reactants_template)]

        # First try to products templates
        product_template_match = False
        for perm in permutations(products_template):
            if perm == rule_products_template:
                product_template_match = True

        # If product templates match
        # find permutations of reactant template that match
        # rule template and keep the indices of those good permutations
        # Else return empty list
        if product_template_match:
            for perm in permutations(reactants_idx_template):
                this_template, this_idx = list(zip(*perm))
                if this_template == rule_reactants_template:
                    matched_idxs.append(this_idx)

    return matched_idxs

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
    parse_issue = False
    for elt in list_of_smiles:
        temp_mol = Chem.MolFromSmiles(elt)
        if temp_mol is None:
            parse_issue = True
            break
        else:
            Chem.rdmolops.RemoveStereochemistry(temp_mol)
            sanitized_smiles.append(Chem.MolToSmiles(temp_mol))
    return sanitized_smiles, parse_issue

def save_json(data, save_to):
    with open(save_to, 'w') as f:
        json.dump(data, f)
    
def load_json(path):
    with open(path, 'r') as f:
        data = json.load(f)
    return data

def postsanitize_smiles(smiles_list):
    """Postsanitize smiles after running SMARTS.
    :returns tautomer list of list of smiles"""

    sanitized_list = []
    tautomer_smarts = '[#7H1X3&a:1]:[#6&a:2]:[#7H0X2&a:3]>>[#7H0X2:1]:[#6:2]:[#7H1X3:3]'

    for s in smiles_list:

        temp_mol = Chem.MolFromSmiles(s, sanitize=False)
        aromatic_bonds = [i.GetIdx() for i in temp_mol.GetBonds() if i.GetBondType() == Chem.rdchem.BondType.AROMATIC]

        for i in temp_mol.GetBonds():
            if i.GetBondType() == Chem.rdchem.BondType.UNSPECIFIED:
                i.SetBondType(Chem.rdchem.BondType.SINGLE)

        try:
            Chem.SanitizeMol(temp_mol)
            Chem.rdmolops.RemoveStereochemistry(temp_mol)
            temp_smiles = Chem.MolToSmiles(temp_mol)

        except Exception as msg:
            if 'Can\'t kekulize mol' in str(msg):
                pyrrole_indices = [i[0] for i in temp_mol.GetSubstructMatches(Chem.MolFromSmarts('n'))]

                # indices to sanitize
                for s_i in pyrrole_indices:
                    temp_mol = Chem.MolFromSmiles(s, sanitize=False)
                    if temp_mol.GetAtomWithIdx(s_i).GetNumExplicitHs() == 0:
                        temp_mol.GetAtomWithIdx(s_i).SetNumExplicitHs(1)
                    elif temp_mol.GetAtomWithIdx(s_i).GetNumExplicitHs() == 1:
                        temp_mol.GetAtomWithIdx(s_i).SetNumExplicitHs(0)
                    try:
                        Chem.SanitizeMol(temp_mol)

                        processed_pyrrole_indices = [i[0] for i in
                                                     temp_mol.GetSubstructMatches(Chem.MolFromSmarts('n'))]
                        processed_aromatic_bonds = [i.GetIdx() for i in
                                                    temp_mol.GetBonds() if i.GetBondType() == Chem.rdchem.BondType.AROMATIC]
                        if processed_pyrrole_indices != pyrrole_indices or aromatic_bonds != processed_aromatic_bonds:
                            continue

                        Chem.rdmolops.RemoveStereochemistry(temp_mol)
                        temp_smiles = Chem.MolToSmiles(temp_mol)
                        break
                    except:
                        continue
                if 'temp_smiles' not in vars():
                    Chem.rdmolops.RemoveStereochemistry(temp_mol)
                    temp_smiles = Chem.MolToSmiles(temp_mol)
                    sanitized_list.append([temp_smiles])
                    continue
            else:
                Chem.rdmolops.RemoveStereochemistry(temp_mol)
                temp_smiles = Chem.MolToSmiles(temp_mol)
                sanitized_list.append([temp_smiles])
                continue
        rxn = AllChem.ReactionFromSmarts(tautomer_smarts)

        try:
            tautomer_mols = rxn.RunReactants((Chem.MolFromSmiles(temp_smiles), ))
        except:
            try:
                tautomer_mols = rxn.RunReactants((Chem.MolFromSmiles(temp_smiles, sanitize=False),))
            except:
                continue

        tautomer_smiles = [Chem.MolToSmiles(m[0]) for m in tautomer_mols]
        sanitized_list.append(sorted(set(tautomer_smiles + [temp_smiles])))

    return list(product(*sanitized_list))