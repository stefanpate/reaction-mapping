import csv
import pythoncyc as pc

fn = 'polyketide_target_list_olano_2009.txt' # List of pk names

# Save filenames
olano_name_cpd_id_fn = 'olano_name_cpd_id.csv'
metacyc_name_cpd_id_fn = 'metacyc_name_cpd_id.csv'
olano_name_cpd_id = []
metacyc_name_cpd_id = []

# Pull data from Metacyc
meta = pc.so('meta')
compounds = meta.compounds
polyketides = meta['|POLYKETIDE|']

# Create compound lookup tables
cpd_id2name = {}
cpd_name2id = {}
for elt in compounds:
    cpd_id2name[elt.frameid] = elt.common_name
    cpd_name2id[elt.common_name] = elt.frameid

# Search in metacyc for each Olano compound
with open(fn, 'r') as pk_txt:
    for i, olano_name in enumerate(pk_txt):
        olano_name = olano_name.strip('\n')
        olano_name = olano_name.split(' ')
        olano_name[0] = olano_name[0].lower() # Make first word lowercase
        olano_name = ' '.join(olano_name)
        metacyc_id = cpd_name2id.get(olano_name)
        olano_name_cpd_id.append([olano_name, metacyc_id])

# Pull names and ids of metacyc polyketides
for elt in polyketides:
    metacyc_name_cpd_id.append([elt.common_name, elt.frameid])

# Save lists to csv files
with open(olano_name_cpd_id_fn, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(olano_name_cpd_id)

with open(metacyc_name_cpd_id_fn, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(metacyc_name_cpd_id)
