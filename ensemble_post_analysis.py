import os
import json
import numpy as np
from ase.io import read
from ase.data import covalent_radii
from ase.neighborlist import NeighborList
from scipy.cluster.hierarchy import fclusterdata

# Overall parameters
METHODS = ["GFN1-xTB", "MACE-MP", "CHGNet"]

MAX_ENERGY_WINDOW = 0.8   # eV (per method)
CLUSTER_CUTOFF = 1.0      # clustering threshold in feature space
BOND_SCALE = 1.2          # covalent radii scaling for bond detection
ENERGY_WEIGHT = 0.1       # scaling factor for energy differences in features


# Bond graph builder for later reactivity detection
def build_bond_graph(atoms, scale=BOND_SCALE):
    cutoffs = [covalent_radii[atoms[i].number] * scale for i in range(len(atoms))]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    bonds = set()
    for i in range(len(atoms)):
        neighbors, _ = nl.get_neighbors(i)
        for j in neighbors:
            if j > i:
                bonds.add((i, j))
    return bonds


# Load slab metadata
with open("../surface_atoms.json") as f:
    slab_info = json.load(f)
n_slab_atoms = slab_info["total_atoms"]


# Build reference bond graph from gas-phase adsorbate
gas_ads = read("../gas/POSCAR")
n_ads_ref = len(gas_ads)
initial_bonds = build_bond_graph(gas_ads)

print(f"Reference adsorbate atoms: {n_ads_ref}")
print(f"Reference bond count: {len(initial_bonds)}")


# Reactivity detection helper
def is_reactive(atoms):
    ads = atoms[n_slab_atoms:]
    if len(ads) != n_ads_ref:
        return True
    current_bonds = build_bond_graph(ads)
    for bond in initial_bonds:
        if bond not in current_bonds:
            return True
    return False


# Load structures and energies
confs = []

for d in sorted(x for x in os.listdir(".") if x.startswith("conf_")):
    json_file = os.path.join(d, "ensemble_screen.json")
    if not os.path.exists(json_file):
        continue

    with open(json_file) as jf:
        results = json.load(jf)

    energy_map = {r["method"]: r["energy"] for r in results}

    for method in METHODS:
        xyz_file = os.path.join(d, f"relaxed_{method}.xyz")
        if not os.path.exists(xyz_file):
            continue

        energy = energy_map.get(method, None)
        if energy is None:
            continue

        atoms = read(xyz_file)
        confs.append({
            "conf": d,
            "method": method,
            "energy": float(energy),
            "atoms": atoms
        })

print(f"\nTotal loaded structures: {len(confs)}")


# Method-specific energy filtering
filtered = []

for method in METHODS:
    method_confs = [c for c in confs if c["method"] == method]
    if not method_confs:
        continue

    energies = np.array([c["energy"] for c in method_confs])
    emin = energies.min()

    for c in method_confs:
        if c["energy"] - emin <= MAX_ENERGY_WINDOW:
            filtered.append(c)

print(f"Number of filtered structures: {len(filtered)}")
if len(filtered) == 0:
    print("No structures survived energy filtering.")
    exit()


# Feature construction
# COM + first 2 heavy atoms vector (with reactive/single-atom fallbacks)
features = []
names = []

for c in filtered:
    atoms = c["atoms"]
    ads = atoms[n_slab_atoms:]
    reactive_flag = is_reactive(atoms)

    # Center of mass for translational
    com = ads.get_center_of_mass()

    # Orientation vector for rotational
    symbols = ads.get_chemical_symbols()
    heavy = [i for i, s in enumerate(symbols) if s != "H"]

    if reactive_flag or len(heavy) == 0:
        v = np.zeros(3)
    elif len(heavy) == 1:
        v = ads.positions[heavy[0]] - com
        norm = np.linalg.norm(v)
        if norm > 1e-8:
            v /= norm
        else:
            v = np.zeros(3)
    else:
        p1 = ads.positions[heavy[0]]
        p2 = ads.positions[heavy[1]]
        v = p2 - p1
        norm = np.linalg.norm(v)
        if norm > 1e-8:
            v /= norm
        else:
            v = np.zeros(3)

    # Energy weighting
    method_confs = [fc for fc in filtered if fc["method"] == c["method"]]
    emin = min(fc["energy"] for fc in method_confs)
    energy_offset = (c["energy"] - emin) * ENERGY_WEIGHT

    feature = np.concatenate([com + energy_offset, v])
    features.append(feature)
    names.append(f"{c['conf']}|{c['method']}")

features = np.vstack(features)


# Clustering
if len(features) == 1:
    labels = np.array([1])
else:
    labels = fclusterdata(features, t=CLUSTER_CUTOFF, criterion="distance")

clusters = {}
for name, label in zip(names, labels):
    clusters.setdefault(int(label), []).append(name)


# Representative selection
rep_data = []

for cluster_id, members in clusters.items():
    best = min(
        members,
        key=lambda m: next(
            c["energy"] for c in filtered
            if f"{c['conf']}|{c['method']}" == m
        )
    )

    conf_name, method = best.split("|")
    entry = next(
        c for c in filtered
        if c["conf"] == conf_name and c["method"] == method
    )

    reactive_flag = is_reactive(entry["atoms"])

    rep_data.append({
        "cluster": cluster_id,
        "conf": conf_name,
        "method": method,
        "energy": entry["energy"],
        "reactive": reactive_flag
    })


# Save metadata
summary = {
    "n_loaded": len(confs),
    "n_filtered": len(filtered),
    "n_clusters": len(clusters),
    "representatives": rep_data,
    "clusters": clusters
}

with open("ensemble_screening_summary.json", "w") as f:
    json.dump(summary, f, indent=2)


# Print summary table
print("\nTop representative candidates:")
print("-" * 75)
print(f"{'Cluster':<8} {'Conf':<10} {'Method':<10} "
      f"{'Energy (eV)':<15} {'Reactive':<8}")
print("-" * 75)

for r in sorted(rep_data, key=lambda x: x["energy"]):
    print(f"{r['cluster']:<8} {r['conf']:<10} {r['method']:<10} "
          f"{r['energy']:<15.6f} {str(r['reactive']):<8}")

print("-" * 75)
print(f"Clusters found: {len(clusters)}")

