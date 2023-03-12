import csv
import requests
import re

n_entries = 33995 # Total number of compounds on seed as of 11/15/2022
seed_cpds_url = f"https://modelseed.org/solr/compounds/select?wt=json&fl=name,id,formula,mass,abbreviation,deltag,deltagerr,charge,aliases,ontology&q=*&rows={n_entries}&sort=id%20asc"
save_to = 'metacyc_to_seed_id.csv'

# Get JSON object from model SEED
json_obj = requests.get(seed_cpds_url).json()
seed_cpds = json_obj['response']['docs']

# Pull out seed and metacyc ids
metacyc_to_seed = []
for cpd in seed_cpds:
    try:
        aliases = cpd['aliases']
        for elt in aliases:
            if 'MetaCyc' in elt:
                metacyc_str = elt.split(':')[1].lstrip(' ')
                metacyc_ids = metacyc_str.split(';')
                metacyc_ids = ['|' + elt.lstrip(' ').upper() + '|' for elt in metacyc_ids]
                break
            else:
                metacyc_ids = []
    except KeyError:
        metacyc_ids = []
    metacyc_to_seed.append([cpd['id']] + metacyc_ids)

# Save seed to metacyc id map
with open(save_to, 'w') as f:
    writer = csv.writer(f)
    for elt in metacyc_to_seed:
        writer.writerow(elt)