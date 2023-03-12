import csv
import requests
import re
import json

n_entries = 33995 # Total number of compounds on seed as of 11/15/2022
seed_cpds_url = f"https://modelseed.org/solr/compounds/select?wt=json&fl=name,id,formula,mass,abbreviation,deltag,deltagerr,charge,aliases,ontology&q=*&rows={n_entries}&sort=id%20asc"
old_cofactor_list_path = 'cofactor_list_alldb.tsv'
old_cofactor_pairs_path = 'cofactor_pair_alldb.json'
new_cofactor_list_path = 'cofactor_list_seed_ids.csv'
new_cofactor_pairs_path = 'cofactor_pairs_seed_ids.json'

# Get JSON object from model SEED website
json_obj = requests.get(seed_cpds_url).json()
seed_cpds = json_obj['response']['docs']

# Read in contents of old files
old_cofactor_list = []
with open(old_cofactor_list_path, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        old_cofactor_list.append(row[0].split('\t'))

with open(old_cofactor_pairs_path, 'r') as f:
    old_cofactor_pairs = json.load(f)

def get_seed_id(key, seed_cpds=seed_cpds):
    '''
    Returns model seed id give key that may be
    model seed id, metacyc id, or compound name.
    '''
    key = key.lower()
    
    if (key[:3] == 'cpd') & (len(key.split('-')) == 1): # Already seed id
        return key
    elif key[:3] == 'cpd': # Metacyc id
        for cpd in seed_cpds:
            
            if 'aliases' in cpd.keys():
                aliases = cpd['aliases']
                for elt in aliases:
                    
                    if 'MetaCyc' in elt:
                        metacyc_str = elt.split(':')[1].lstrip(' ').lower()
                        metacyc_ids = metacyc_str.split(';')
                        metacyc_ids = [elt.lstrip(' ') for elt in metacyc_ids]
                        break
            
                if key in metacyc_ids:
                    return cpd['id']
    
    else: # Search of compound name
        for cpd in seed_cpds:

            if 'aliases' in cpd.keys():
                aliases = cpd['aliases']
                for elt in aliases:

                    if 'Name' in elt:
                        names_str = elt.split(':')[1].lstrip(' ').lower()
                        names = names_str.split(';')
                        names = [name.lstrip(' ') for name in names]
                        break

                if key in names:
                    return cpd['id']
    
    return None # If no seed id found

# Create new cofactor list w/ model seed ids
new_cofactor_seeds = []
new_cofactor_names = []
new_cofactor_keys = []
for elt in old_cofactor_list:
    name, key = elt
    seed_id = get_seed_id(name)
    if (seed_id is not None) & (seed_id not in new_cofactor_seeds):
        new_cofactor_seeds.append(seed_id)
        new_cofactor_names.append(name)
        new_cofactor_keys.append(key)

new_cofactor_list = list(zip(new_cofactor_seeds, new_cofactor_keys, new_cofactor_names))

# Create new cofactor pairs dict w/ model seed ids
new_cofactor_pairs = {}
for k, v in old_cofactor_pairs.items():
    new_cofactor_pairs[k] = [v[0]] # Must keep first pair value for Josephs cofactor_name_dict
    for pair in v:
        reactant_seed_id = get_seed_id(pair[1])
        product_seed_id = get_seed_id(pair[2])
        if (reactant_seed_id is not None) & (product_seed_id is not None):
            new_cofactor_pairs[k].append([pair[0], reactant_seed_id, product_seed_id])

# Save new cofactor files
with open(new_cofactor_list_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(new_cofactor_list)

with open(new_cofactor_pairs_path, 'w') as f:
    json.dump(new_cofactor_pairs, f)