from mapping_fcns import load_json, save_json

# Set
rules = 'minimal1224_all_uniprot.tsv' # Path to rules
all_rxns_path = '100_random_mc_v21_rxns_rnd_seed_1234.json'
rxns_prefix = 'mc_v21_rxns_batch_' # Prefix of path to save reaction batches
mapping_prefix = 'mapping_mc_v21_rxns_min_cof_rules_batch_' # Prefix of path to save
submit_prefix = 'submit_batch_'
n = 10 # Number of batches

# Load reaction data
mc_rxns = load_json(all_rxns_path)

default_text_1 = ['#!/bin/bash', '#SBATCH -A b1039', '#SBATCH -p b1039',
                    '#SBATCH -N 1', '#SBATCH -n 15', '#SBATCH --mem-per-cpu=8G',
                    '#SBATCH -t 72:00:00'
                ] 

default_text_2 = ['ulimit -c 0', 'module load python/anaconda3.6', 'source activate mine']

# Split and save
batch_size = int(len(mc_rxns) / n)
ids = list(mc_rxns.keys())
for i in range(n):
    # Split into n batches
    if i < n - 1:
        this_ids = ids[i * batch_size: (i + 1) * batch_size]
    else:
        this_ids = ids[i * batch_size:] # Catch remainder at the end

    this_batch = {elt: mc_rxns[elt] for elt in this_ids}

    # Get file names
    submit = submit_prefix + f"{i}.sh"
    rxns = rxns_prefix + f"{i}.json"
    mapping = mapping_prefix + f"{i}.csv"
    job_name = f"rxn_mapping_{i}"
    outlog = f"outlog_{i}"
    errlog = f"errlog_{i}"

    job_name_log = ['#SBATCH --job-name="' + job_name + '"', f"#SBATCH -o {outlog}", f"#SBATCH -e {errlog}"]

    custom_text = ["python map_rxns.py {} {} {}".format(rules, rxns, mapping)] # Last line of .sh

    # Save n reaction batches
    save_json(this_batch, rxns)

    # Write n submit scripts
    with open(submit, 'w') as f:
        for line in default_text_1 + job_name_log + default_text_2 + custom_text:
            f.write(line)
            f.write('\n')
