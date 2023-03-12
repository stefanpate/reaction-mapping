import csv
import pythoncyc as pc

mapped_rxns_path = 'mapped_reactions_simple.csv'
all_coreactants_path = 'coreactant_counts_all_rxns.csv'
mapped_coreactants_path = 'coreactant_counts_mapped_rxns.csv'
unmapped_coreactants_path = 'coreactant_counts_unmapped_rxns.csv'
pk_cpd_id_path = 'combined_pk_cpd_id.csv'
save_paths = [all_coreactants_path, mapped_coreactants_path, unmapped_coreactants_path]
meta = pc.so('meta') # Access metacyc

def count_coreactants(rxns, meta=meta):
    '''
    Given list of metacyc rxn ids, 
    enumerate frequency of every substrate.
    '''
    coreactant_dict = {}
    for rxn in rxns:
        rxn_frame = meta[rxn]
        substrates = rxn_frame.left + rxn_frame.right
        for sub in substrates:
            common_name = meta[sub].common_name
            if ((sub, common_name) in coreactant_dict) & (sub not in pks):
                    coreactant_dict[(sub, common_name)] += 1
            else:
                coreactant_dict[(sub, common_name)] = 1            

    coreactant_freq = [[elt[1][0], elt[1][1], elt[0]] for elt in sorted(zip(coreactant_dict.values(), coreactant_dict.keys()), reverse=True)]
    return coreactant_freq

# Read in pk cpd ids
pks = []
with open(pk_cpd_id_path, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        pks.append(row[0])

# Read in reactions
rxns = []
mapped_rxns = []
unmapped_rxns = []
with open(mapped_rxns_path, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        rxns.append(row[0])
        if len(row) == 1:
            unmapped_rxns.append(row[0])
        else:
            mapped_rxns.append(row[0])

# Pull coreactants out and count them
for i, elt in enumerate([rxns, mapped_rxns, unmapped_rxns]):
    coreactant_count = count_coreactants(elt) # Count and sort

    # Save each group
    with open(save_paths[i], 'w') as f:
        writer = csv.writer(f)
        writer.writerows(coreactant_count)