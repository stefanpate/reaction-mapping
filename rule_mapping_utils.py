import pandas as pd
import itertools
from rdkit import Chem
from rdkit.Chem import AllChem
import copy
import json
import os
import csv

class MapRules:
    """Map intended reactants to database reactants described by the same reaction operators."""

    def __init__(self, rules_path=None, molfiles_path=None, seed_dict=None,
                 cofactor_list_path=None, cofactor_pair_path=None, cofactors=True, start_name=None, start_smiles=None):
        self.rules = pd.read_csv(rules_path, sep='\t', index_col=0)
        self.molfiles_path = molfiles_path
        if cofactors:
            self.cofactor_name_dict, self.cofactor_list_dict, self.cofactor_pair_dict = get_cofactors(cofactor_list_path, cofactor_pair_path)
        else:
            self.cofactor_name_dict = {}
            self.cofactor_list_dict = {}
            self.cofactor_pair_dict = {}
        self.seed_dict = seed_dict
        if start_name:
            self.start_name = start_name
        if start_smiles:
            self.start_smiles = start_smiles

    def map_pickaxe_rules(self, lhs_dict, rhs_dict, rule_current, return_reaction_center=False):

        rxn_df = self.rules.loc[rule_current.split(';')[0]]
        rule = rxn_df['SMARTS']
        reactants = rxn_df['Reactants']
        products = rxn_df['Products']

        # remove cofactor, sanitize mols
        lhs_list, rhs_list = self._process_substrates(lhs_dict, rhs_dict, rule_current)

        # match index with pickaxe
        found_match = self._map_rules(rule, lhs_list, rhs_list, reactants, products, return_reaction_center)

        return found_match

    def _map_rules(self, rule, lhs, rhs, reactants, products, return_reaction_center):
        """Operator mapping"""

        rxn = Chem.rdChemReactions.ReactionFromSmarts(rule)
        reactants = reactants.split(';')
        cofactor_index_reactants = [i for i, r in enumerate(reactants) if r != 'Any']

        products = products.split(';')
        cofactor_index_products = [i for i, p in enumerate(products) if p != 'Any']

        # if number of reactants does not match reactant template
        if len(lhs) > reactants.count('Any'):
            repetitive_mols = set(lhs).intersection(set(rhs))

            while repetitive_mols:
                lhs.remove(sorted(repetitive_mols)[0])
                rhs.remove(sorted(repetitive_mols)[0])
                repetitive_mols = set(lhs).intersection(set(rhs))

        lhs_set = set()
        for lhs_perm in itertools.permutations(lhs):
            lhs_set.add(lhs_perm)

        for lhs_perm in lhs_set:
            lhs_temp = list(lhs_perm)

            for c in cofactor_index_reactants:
                if self.molfiles_path:
                    lhs_temp[c:c] = [Chem.MolToSmiles(Chem.MolFromMolFile(os.path.sep.join([self.molfiles_path, self.cofactor_name_dict[reactants[c]] + '.mol'])))]
                elif self.seed_dict:
                    lhs_temp[c:c] = [self.seed_dict[self.cofactor_name_dict[reactants[c]]]]

            # pruned MetaCyc
            try:
                lhs_tuple = tuple([Chem.MolFromSmiles(i) for i in lhs_temp])
                outputs = rxn.RunReactants(lhs_tuple)
            except:
                try:
                    lhs_tuple = tuple([Chem.MolFromSmiles(i, sanitize=False) for i in lhs_temp])
                    outputs = rxn.RunReactants(lhs_tuple)
                except:
                    continue

            for rxn_output in outputs:

                rhs_run = [Chem.MolToSmiles(rhs_mols) for rhs_mols in rxn_output]
                rhs_list = copy.deepcopy(rhs_run)

                for c in cofactor_index_products:
                    rhs_list.remove(rhs_run[c])

                # for all tautomer possibilities of clean rhs
                for rhs in postsanitize_smiles(rhs):
                    rhs = list(rhs)

                    for rhs_list in postsanitize_smiles(rhs_list):

                        # pruned MetaCyc
                        if sorted(list(rhs_list)) == sorted(rhs):
                            return True

        return False

    def _post_process(self, enzyme_list):
        """Post processing"""

        enzyme_list_temp = []

        non_orphan_flag = False

        for e in enzyme_list:
            if 'ENZRXN' in e:
                non_orphan_flag = True
                continue

        if non_orphan_flag:
            for e in enzyme_list:
                if 'ENZRXN' in e:
                    enzyme_list_temp.append(e)
            return enzyme_list_temp
        else:
            return enzyme_list

    def _process_substrates(self, lhs_dict_temp, rhs_dict_temp, rule_current):
        """Process substrates"""

        # check cofactor designation
        rule_reactant_names = self.rules.loc[rule_current, 'Reactants']
        rule_product_names = self.rules.loc[rule_current, 'Products']

        reactant_names, product_names = label_cofactor(sorted(lhs_dict_temp), sorted(rhs_dict_temp), self.cofactor_list_dict, self.cofactor_pair_dict)

        if sorted(rule_reactant_names.split(';')) != sorted(reactant_names.split(';')) \
                or sorted(rule_product_names.split(';')) != sorted(product_names.split(';')):
            raise ValueError('Cofactor designation error.')

        # create list from dict
        lhs_list = [lhs_dict_temp[k] for i, k in enumerate(sorted([k for k in lhs_dict_temp]))
                    if reactant_names.split(';')[i] == 'Any']
        rhs_list = [rhs_dict_temp[k] for i, k in enumerate(sorted([k for k in rhs_dict_temp]))
                    if product_names.split(';')[i] == 'Any']

        # sanitize
        for i, m in enumerate(lhs_list):
            try:
                temp_mol = Chem.MolFromSmiles(m)
                Chem.rdmolops.RemoveStereochemistry(temp_mol)
                lhs_list[i] = Chem.MolToSmiles(temp_mol)
            except:
                temp_mol = Chem.MolFromSmiles(m, sanitize=False)
                Chem.rdmolops.RemoveStereochemistry(temp_mol)
                lhs_list[i] = Chem.MolToSmiles(temp_mol)

        for i, m in enumerate(rhs_list):
            try:
                temp_mol = Chem.MolFromSmiles(m)
                Chem.rdmolops.RemoveStereochemistry(temp_mol)
                rhs_list[i] = Chem.MolToSmiles(temp_mol)
            except:
                temp_mol = Chem.MolFromSmiles(m, sanitize=False)
                Chem.rdmolops.RemoveStereochemistry(temp_mol)
                rhs_list[i] = Chem.MolToSmiles(temp_mol)

        return lhs_list, rhs_list
    
def get_cofactors(input_cofactor_list_path, input_cofactor_pair_path):
    """Get cofactor list & pairs"""

    # cofactor to cpd id dict
    cofactor_name_dict = {}

    # cofactor list name designation
    cofactor_list_dict = {}
    with open (input_cofactor_list_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            cofactor_list_dict[row[0].upper()] = row[1]
            if row[1] not in cofactor_name_dict:
                cofactor_name_dict[row[1]] = row[2]
    # for k, v in pd.read_csv(input_cofactor_list_path, sep=',', index_col=0).iterrows():
    #     cofactor_list_dict[k.upper()] = v['replacement']
    #     if v['replacement'] not in cofactor_name_dict:
    #         cofactor_name_dict[v['replacement']] = k

    # cofactor pair name designation
    cofactor_pair_dict = {}
    with open(input_cofactor_pair_path) as f:
        cofactor_pair_read_json = json.loads(f.read())
    for k, v in cofactor_pair_read_json.items():
        for pair in v:
            cofactor_pair_dict[(pair[1].upper(), pair[2].upper())] = k
            cofactor_pair_dict[(pair[2].upper(), pair[1].upper())] = '%s,%s' % (k.split(',')[1], k.split(',')[0])
            if k.split(',')[0] not in cofactor_name_dict:
                cofactor_name_dict[k.split(',')[0]] = pair[1]
                cofactor_name_dict[k.split(',')[1]] = pair[2]

    return cofactor_name_dict, cofactor_list_dict, cofactor_pair_dict

def label_cofactor(reactant_molfile, product_molfile, cofactor_list, cofactor_pair):  # label cofactors
    """Label cofactors & cofactor pairs for product & reactant names"""

    reactant_molfile = [elt.upper() for elt in reactant_molfile]
    product_molfile = [elt.upper() for elt in product_molfile]
    
    # new substrate labels
    reactant_names = ['Any'] * len(reactant_molfile)
    product_names = ['Any'] * len(product_molfile)

    # get cofactor pairs
    for i_lhs, lhs in enumerate(reactant_molfile):
        for i_rhs, rhs in enumerate(product_molfile):

            # skip if already assigned
            if product_names[i_rhs] != 'Any':
                continue

            # assign cofactor pair designation
            try:
                temp_pair = cofactor_pair[(lhs, rhs)]
                reactant_names[i_lhs] = temp_pair.split(',')[0]
                product_names[i_rhs] = temp_pair.split(',')[1]
                break
            except KeyError:
                continue

        # assign cofactor list if no cofactor pair assigned
        if reactant_names[i_lhs] == 'Any':
            try:
                reactant_names[i_lhs] = cofactor_list[lhs]
            except KeyError:
                continue

    # assign cofactor list for rhs
    for i_rhs, rhs in enumerate(product_molfile):

        # assign cofactor list if no cofactor pair assigned
        if product_names[i_rhs] == 'Any':
            try:
                product_names[i_rhs] = cofactor_list[rhs]
            except KeyError:
                continue

    return ';'.join(reactant_names), ';'.join(product_names)

def postsanitize_smiles(smiles_list):
    """Postsanitize smiles after running SMARTS.
    :returns tautomer list of list of smiles"""

    sanitized_list = []
    # tautomer_smarts = '[#6:1]1:[#6:2]:[#7H1X3:3]:[#6:4]:[#7H0X2:5]:1>>[#6:1]1:[#6:2]:[#7H0X2:3]:[#6:4]:[#7H1X3:5]:1'
    tautomer_smarts = '[#7H1X3&a:1]:[#6&a:2]:[#7H0X2&a:3]>>[#7H0X2:1]:[#6:2]:[#7H1X3:3]'

    for s in smiles_list:

        temp_mol = Chem.MolFromSmiles(s, sanitize=False)

        # # pickaxe
        # temp_mol = Chem.rdmolops.RemoveHs(temp_mol)

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
                # unkekulized_indices = [int(i) for i in str(msg).split('Unkekulized atoms: ')[1].split('.')[0].rstrip(' \n').split(' ')]
                pyrrole_indices = [i[0] for i in temp_mol.GetSubstructMatches(Chem.MolFromSmarts('n'))]

                # indices to sanitize
                # for s_i in set(unkekulized_indices).intersection(set(pyrrole_indices)):
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

    return list(itertools.product(*sanitized_list))