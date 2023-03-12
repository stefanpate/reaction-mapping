import csv
import pythoncyc as pc

cpd_id_to_pwy_fn = 'pk_cpd_id_pwy.csv'
starters = ['|PROPIONYL-COA|', '|MALONYL-COA|', '|ACETYL-COA|', '|METHYL-MALONYL-COA|', '|D-METHYL-MALONYL-COA|']
save_to = 'pwy_cpd_has_starter.csv'

meta = pc.so('meta') # Pull metacyc

# Read in polyketide pathways
pk_pwys = {}
with open(cpd_id_to_pwy_fn, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if len(row) == 1:
            pk_pwys[row[0]] = []
        else:
            pk_pwys[row[0]] = row[1:]

# Label pwys that contain starters
pwys_w_starters = []
yes_ct = 0
pwy_ct = 0
for cpd, pwy_list in pk_pwys.items():
    for pwy in pwy_list:
        pwy_ct += 1
        has_starters = False
        rxns = meta[pwy].reaction_list
        for rxn in rxns:
            substrates = meta[rxn].substrates
            has_starters = len(set(substrates) & set(starters)) > 0
            
            if has_starters:
                yes_ct += 1
                break

        pwys_w_starters.append([pwy, cpd, has_starters])
        print([pwy, cpd, has_starters])

print(f"Fraction pathways w/ starters: {yes_ct} / {pwy_ct}")

# Save
with open(save_to, 'w') as f:
    writer = csv.writer(f)
    for elt in pwys_w_starters:
        writer.writerow(elt)