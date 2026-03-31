import os
import glob
import csv
import sys

# User input
if len(sys.argv) != 2:
    print("Usage: python energy_parse.py <METHOD>")
    sys.exit(1)

TARGET_METHOD = sys.argv[1]

# Load SMILES (original & cleaned) for easier calling later
original_map = {}
with open("species.txt") as f:
    next(f) #skip header
    for i, line in enumerate(f):
        ori = line.strip().split()
        original_map[str(i)] = ori[0]

smiles_map = {}
with open("cleaned_smiles.txt") as f:
    for i, line in enumerate(f):
        smiles_map[str(i)] = line.strip()

rows = []

# Iterate through subdirectories
for d in sorted(os.listdir("."), key=lambda x: int(x) if x.isdigit() else float("inf")):
    if not os.path.isdir(d):
        continue

    original = original_map.get(d, "")
    smiles = smiles_map.get(d, "")

    files = glob.glob(os.path.join(d, "std.out"))
    if not files:
        rows.append([d,original, smiles, None, None, None])
        continue

    file = files[0]

    with open(file) as f:
        lines = f.readlines()

    selected = None
    in_table = False

    for line in lines:
        if "Best candidate per method:" in line:
            in_table = True
            continue

        if in_table:
            # skip headers / separators
            if "Conf" in line or "-----" in line:
                continue

            # stop at end of table
            if line.strip() == "":
                break

            parts = line.strip().split()

            # expected format:
            # conf_22 CHGNet -212.988922 False
            if len(parts) >= 4:
                conf, method, energy, reactive = parts[:4]

                if method == TARGET_METHOD:
                    selected = [conf, method, energy, reactive]
                    break

    if selected:
        rows.append([d, original, smiles] + selected)
    else:
        rows.append([d, original, smiles, None, None, None, None])


# Write csv
with open(f"output_{TARGET_METHOD}.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "dir", "original", "smiles", "conf", "method", "energy", "reactive"
    ])
    writer.writerows(rows)
