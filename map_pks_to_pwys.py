import csv
import pythoncyc as pc

# File names
pk_cpd_fn = 'combined_pk_cpd_id.csv'
save_to = 'pk_cpd_id_pwy.csv'

# Read in combined (Olano and Metacyc) list of pk cpd ids
pk_cpds = [] # Store pk cpd ids
with open(pk_cpd_fn, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        pk_cpds.append(row[0])

# Read in reactions and pathways
meta = pc.so('meta')
pathways = meta.pathways

# Gather pathway ids that don't contain 'PWY'
special_pathways = []
for pwy in pathways:
    if 'PWY' not in pwy.frameid:
        special_pathways.append(pwy.frameid)

# Search all pathways for those containing polyketide compounds
pk_pwys = [[] for _ in range(len(pk_cpds))]
for pwy in pathways:
    for rxn_id in pwy.reaction_list:
        if ('PWY' in rxn_id) or (rxn_id in special_pathways): # Skip superpathways
            break
        else:
            rxn = meta[rxn_id]
            for elt in rxn.substrates:
                if elt in pk_cpds:
                    idx = pk_cpds.index(elt)
                    pk_pwys[idx].append(pwy.frameid)

# Eliminate repeat pathways in each sublist of pk_pwys
temp = pk_pwys.copy()
pk_pwys = []
for elt in temp:
    pk_pwys.append(list(set(elt)))

# Save polyketide-associated pathways
with open(save_to, 'w') as f:
    writer = csv.writer(f)
    for i in range(len(pk_pwys)):
        writer.writerow([pk_cpds[i]] + pk_pwys[i])


