import csv

combined_fn = 'combined_pk_cpd_id.csv'

# Files to read in
olano_name_cpd_id_fn = 'olano_name_cpd_id.csv'
metacyc_name_cpd_id_fn = 'metacyc_name_cpd_id.csv'

# Lists for compound ids
olano = []
metacyc = []

olano_total_n = 0 # Total number from Olano review

# Read in compound ids
with open(olano_name_cpd_id_fn, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        olano_total_n += 1
        if row[1]:
            olano.append(row[1])

with open(metacyc_name_cpd_id_fn, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[1]:
            metacyc.append(row[1])

# Count overlap
intersection = list(set(olano) & set(metacyc))

print(f"{len(olano)} / {olano_total_n} of Olano compounds found in Metacyc")
print(f"Metacyc & Olano agree on {len(intersection)} polyketides...")
print(f"...out of {len(olano)} in Olano, and {len(metacyc)} in Metacyc")

# Save union of lists
union = list(set(olano) | set(metacyc))
temp = [[elt,] for elt in union]
union = temp
with open(combined_fn, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(union)