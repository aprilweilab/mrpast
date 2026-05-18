# Given a .ped file and a sample list (one sample ID per line), choose only one sample from
# each family cluster (shared FamilyID) and print the number of families found.
import sys
import csv

if len(sys.argv) < 3:
    print("Usage: ped_no_families.py <ped file> <sample file>", file=sys.stderr)
    exit(1)

ped_filename = sys.argv[1]
sample_filename = sys.argv[2]

indiv2family = {}
with open(ped_filename, newline='') as ped:
    reader = csv.DictReader(ped, delimiter="\t")
    for row in reader:
        indiv2family[row["Individual ID"].strip()] = row["Family ID"].strip()

seen_families = set()
with open(sample_filename) as f:
    for line in f:
        sample_id = line.strip()
        family_id = indiv2family[sample_id]
        if family_id in seen_families:
            print(f"Dropping sample {sample_id} since family {family_id} already seen", file=sys.stderr)
        else:
            seen_families.add(family_id)
            print(sample_id)
