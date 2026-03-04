#!/bin/bash


# AutoFlow by Adhika Setiawan
# ---------------------------------------------- 
## Auxiliary information

# Strict mode
set -euo pipefail

# Parallelization
MAXPARJOB=4

# Execution timer
start_time=$(date +%s)

trap '
  end_time=$(date +%s)
  elapsed=$((end_time - start_time))
  printf "\nTotal execution time: %02d:%02d:%02d\n" \
    $((elapsed/3600)) $((elapsed%3600/60)) $((elapsed%60))
' EXIT

# Atomic number hash lookup table. On average, O(1)
declare -A Z=(
  [H]=1  [He]=2  [Li]=3  [Be]=4  [B]=5   [C]=6   [N]=7   [O]=8   [F]=9   [Ne]=10
  [Na]=11 [Mg]=12 [Al]=13 [Si]=14 [P]=15 [S]=16 [Cl]=17 [Ar]=18
  [K]=19  [Ca]=20 [Sc]=21 [Ti]=22 [V]=23 [Cr]=24 [Mn]=25 [Fe]=26 [Co]=27 [Ni]=28
  [Cu]=29 [Zn]=30 [Ga]=31 [Ge]=32 [As]=33 [Se]=34 [Br]=35 [Kr]=36
  [Rb]=37 [Sr]=38 [Y]=39 [Zr]=40 [Nb]=41 [Mo]=42 [Tc]=43 [Ru]=44 [Rh]=45 [Pd]=46 [Ag]=47
  [Cd]=48 [In]=49 [Sn]=50 [Sb]=51 [Te]=52 [I]=53 [Xe]=54
  [Cs]=55 [Ba]=56 [La]=57 [Ce]=58 [Pr]=59 [Nd]=60 [Pm]=61 [Sm]=62 [Eu]=63 [Gd]=64 [Tb]=65
  [Dy]=66 [Ho]=67 [Er]=68 [Tm]=69 [Yb]=70 [Lu]=71
  [Hf]=72 [Ta]=73 [W]=74 [Re]=75 [Os]=76 [Ir]=77 [Pt]=78 [Au]=79 [Hg]=80
  [Tl]=81 [Pb]=82 [Bi]=83 [Po]=84 [At]=85 [Rn]=86
  [Fr]=87 [Ra]=88 [Ac]=89 [Th]=90 [Pa]=91 [U]=92 [Np]=93 [Pu]=94 [Am]=95 [Cm]=96 [Bk]=97
  [Cf]=98 [Es]=99 [Fm]=100 [Md]=101 [No]=102 [Lr]=103
  [Rf]=104 [Db]=105 [Sg]=106 [Bh]=107 [Hs]=108 [Mt]=109 [Ds]=110 [Rg]=111 [Cn]=112
  [Fl]=114 [Lv]=116
)

# Add --help alias for UX
if [[ "${1:-}" == "--help" ]]; then
    set -- -h
fi

# Activate python environment
module load miniconda3
conda activate autoflow

# ---------------------------------------------- 
## Argument processing
# Initialize script arguments
slabspec=""
miller=""
adsspec=""
lattconst="None"
packing="fcc"

# getopts to parse options
while getopts ":s:m:a:l:p:h" opt; do
	case $opt in
	h)
		cat << EOF
##################################################
###         AutoFlow Bash Script v.0.7.        ###
##################################################

The AutoFlow script is intended for automating the
generation of initial adsorption structures of the
given adsorbate on the slab of the given element.
At the end of the operation, the script
generates structure files for the gaseous molecule,
the clean slab, and enumerated adsorption
configurations, along with the corresponding VASP
input files to be used for later optimization
using DFT or otherwise. A hybrid screening scheme
based on forces is used to initially filter out
unphysical solutions using a combination of
Grimme's GFN-FF -> GFN1-xTB methods.

Usage: $(basename "$0") -s SLAB -m H,K,L -a SMILES
       [-l LATTCONST] [-p PACKING] [-h/--help]

Required options:
  -s  Slab element (e.g. Cu, Pt).
  -m  Comma-separated Miller indices (e.g. 1,1,1).
  -a  Adsorbate SMILES string (e.g. CH3COOH).

Optional options:
  -l  Lattice constant.
  -p  Packing/crystal structure type. Must be one
      of the following: fcc, hcp, bcc, bct.
      Defaults to fcc when not specified.
  -h  Print this help message and exit. When
      this option is called, all other options
      are ignored, and only the help message is
      printed. Can also be called with --help.

Prepared by Adhika Setiawan
(adhika.n.suryatin@gmail.com)
EOF

    		exit 0
    		;;
	s)
  		slabspec="$OPTARG"
  		slabspec="${slabspec^}" # Foolproof capitalization

		# Ensure that specified element is valid
		if [[ -z "${Z[$slabspec]:-}" ]]; then
			echo "Option Argument Error: -s must be a valid element symbol (e.g. Cu, Pt, Al). You entered: $slabspec"
			exit 1
		fi
		;;
	m)
		miller="$OPTARG"
		
		# Miller index format check
		if ! [[ "$miller" =~ ^-?[0-9]+,-?[0-9]+,-?[0-9]+$ ]]; then
        		echo "Option Argument Error: -m must be three comma-separated integers (e.g. 1,1,1)"
        		exit 1
      		fi

		IFS=',' read -r -a numbers_array <<< "$miller"
  		;;
	a)
		adsspec="$OPTARG"
		;;
	l)
		# Lattice constant format check
		if ! [[ "$OPTARG" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        		echo "Option Argument Error: -l must be a floating-point number"
        		exit 1
    		fi

    		lattconst="$OPTARG"
		;;
	p)
		packing="$OPTARG"
		# Packing type check
		case "$packing" in
        		fcc|hcp|bcc|bct) ;;
        		*)
            			echo "Option Argument Error: -p must be one of: fcc, hcp, bcc, bct"
				# Extend as needed
            			exit 1
            			;;
    		esac
		;;
	\?)
  		echo "Invalid option: -$OPTARG"
		exit 1
  		;;
	:)
		echo "Incomplete argument(s) supplied."
		exit 1
  		;;
	esac
done

# Require all mandatory options
if [[ -z "$slabspec" || -z "$miller" || -z "$adsspec" ]]; then
	echo "Incomplete arguments supplied"
	exit 1
fi

# Print output arguments
echo "$packing Slab Element: $slabspec"
echo "Miller Indices of Simulated Facet: (${numbers_array[0]},${numbers_array[1]},${numbers_array[2]})"
echo "Adsorbate SMILES String: $adsspec"

# ---------------------------------------------- 
## Input processing
# Protective quote stripping
adsspec="${adsspec%\"}"
adsspec="${adsspec#\"}"
adsspec="${adsspec%\'}"
adsspec="${adsspec#\'}"
adsspec="$(echo -e "${adsspec}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

# Extract element symbols from SMILES input
smiles=$adsspec
# RDKit is used to validate SMILES input
read -r -a sorted_elements <<< "$(python3 - << EOF
from rdkit import Chem
smiles = """$adsspec"""
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    raise ValueError(f"Invalid SMILES: {smiles}")
elements = sorted(
    {atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 0},
    key=lambda el: Chem.GetPeriodicTable().GetAtomicNumber(el)
)
print(" ".join(elements))
EOF
)"

# ---------------------------------------------- 
## VASP files for gaseous adsorbate
# Acquire and concatenate gaseous POTCARs
i=0 # Initialize variable to not trip -u
for el in "${sorted_elements[@]}"; do
  cp -f ~/../../share/Apps/vasp/5.4.4.pl2/potentials/PBE.54/$el/POTCAR \
     "$(printf 'POTCAR_%02d' "$i")"
  ((++i))
done
cat POTCAR_*>POTCAR
rm -f POTCAR_*

# Prepare gaseous folder
mkdir ./gas
mv ./POTCAR ./gas
cd ./gas

# Call RDKit & ASE to parse SMILES via Python
cat > gas_gen.py << 'EOF'
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from ase.atoms import Atoms
from ase.build import sort
from ase.io import write
from ase.constraints import FixAtoms
import numpy as np

# SMARTS used to match functional groups
SMARTS_PATTERNS = [
    ("alkene", Chem.MolFromSmarts("C=C")),
    ("alkyne", Chem.MolFromSmarts("C#C")),
    ("carbonyl", Chem.MolFromSmarts("[CX3]=[OX1]")),
    ("alcohol", Chem.MolFromSmarts("[OX2H]")),
    ("amine", Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")),
    ("nitrile", Chem.MolFromSmarts("C#N")),
    ("aromatic_ring", Chem.MolFromSmarts("a1aaaaa1")),
    ("co2", Chem.MolFromSmarts("O=C=O")),
]

def detect_anchor_groups(mol, mol_no_dummy, old_to_new_idx, adsorption_sites):
    """
    Returns anchor groups as lists of atom indices (in mol_no_dummy indexing).

    Automatic selection of point of adsorption is done via filters.
    Three tiers of filters are used:
    1) Open-shell species / intermediates
    2) Closed-shell species with key functional group(s)
    3) Gasteiger charge priority fallback
    """

    anchor_groups = []

    # Dummy atom notation / open-shell species filtering
    if adsorption_sites:
        for site_id, heavy_idx in enumerate(adsorption_sites):
            anchor_groups.append([old_to_new_idx[heavy_idx]])
        return anchor_groups

    # SMARTS pattern matching for closed-shell
    for name, pattern in SMARTS_PATTERNS:
        matches = mol_no_dummy.GetSubstructMatches(pattern)
        for match in matches:
            anchor_groups.append(list(match))

    if anchor_groups:
        return anchor_groups

    # Polar center / Gasteiger charge fallback
    try:
        AllChem.ComputeGasteigerCharges(mol_no_dummy)
        charges = []
        for atom in mol_no_dummy.GetAtoms():
            q = atom.GetProp("_GasteigerCharge")
            charges.append(float(q))
        idx = int(np.argmax(np.abs(charges)))
        anchor_groups.append([idx])
        return anchor_groups
    except Exception:
        return []

    return anchor_groups


# Suppress explicit H warning
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# SMILES string. SMILES injected separately for safety
smiles = """__SMILES__"""

# SMILES -> RDKit molecule
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    raise ValueError(f"Invalid SMILES string: {smiles}")

# Turn on auto add H atoms
mol = Chem.AddHs(mol)

# Identify dummy adsorption sites
dummy_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]

# Map: heavy_atom_idx -> list(dummy_atom_idx)
ads_site_map = {}

for d in dummy_atoms:
    neighbors = [n for n in d.GetNeighbors() if n.GetAtomicNum() != 0]
    if len(neighbors) != 1:
        raise ValueError("Each dummy atom must bond to exactly one heavy atom")
    heavy_idx = neighbors[0].GetIdx()
    ads_site_map.setdefault(heavy_idx, []).append(d.GetIdx())

adsorption_sites = list(ads_site_map.keys())

# Geometry bias to connect representative dummy atoms
if len(adsorption_sites) >= 2:
    rw = Chem.RWMol(mol)
    # Pick one dummy per adsorption site
    rep_dummies = [ads_site_map[h][0] for h in adsorption_sites]

    # Chain-connect them (avoids overbonding)
    for i in range(len(rep_dummies) - 1):
        a1 = rep_dummies[i]
        a2 = rep_dummies[i + 1]
        if not rw.GetBondBetweenAtoms(a1, a2):
            rw.AddBond(a1, a2, Chem.BondType.SINGLE)
    mol = rw.GetMol()

# Remove dummy atoms for further processing
rw = Chem.RWMol(mol)
for idx in sorted([a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 0], reverse=True):
    rw.RemoveAtom(idx)

mol_no_dummy = rw.GetMol()

# Embed molecule
params = AllChem.ETKDGv3()
res = AllChem.EmbedMolecule(mol_no_dummy, params)
if res != 0:
    raise ValueError("3D embedding of molecule failed")

# Optional: pre-optimize geometry
# AllChem.UFFOptimizeMolecule(mol_no_dummy)

# RDKit Mol -> ASE Atoms
conf = mol_no_dummy.GetConformer()
positions = conf.GetPositions()
symbols = [atom.GetSymbol() for atom in mol_no_dummy.GetAtoms()]
atoms = Atoms(symbols=symbols, positions=positions)

# Map heavy atoms to adsorption sites
# old_to_new_idx: mol heavy atom index -> mol_no_dummy index
old_to_new_idx = {}
j = 0
for i, a in enumerate(mol.GetAtoms()):
    if a.GetAtomicNum() != 0:
        old_to_new_idx[i] = j
        j += 1

anchor_groups = detect_anchor_groups(
    mol,
    mol_no_dummy,
    old_to_new_idx,
    adsorption_sites
)

tags = np.full(len(atoms), -1, dtype=int)

for group_id, group in enumerate(anchor_groups):
    for idx in group:
        tags[idx] = group_id

atoms.set_tags(tags)

# Sort atoms (tags automatically permuted)
sorted_atoms = sort(atoms)

# Finalize cell, PBC, centering
sorted_atoms.set_pbc(False)
sorted_atoms.set_cell([20.0, 20.0, 20.0])
sorted_atoms.center()

# Write output
# meta.xyz generated for carrying atom tag metadata
write("meta.xyz", sorted_atoms)
write("POSCAR", sorted_atoms, format="vasp")
EOF

# Safely inject the SMILES string
sed -i "s|__SMILES__|$smiles|" gas_gen.py

# Run Python code
python gas_gen.py || {
    echo "Error: gas-phase POSCAR generation failed"
    exit 1
} # Safety to abort on failure so errors do not silently propagate
rm -f gas_gen.py

# Prepare KPOINTS, INCAR, SLURM submission script
# Generate gas KPOINTS File
cat > KPOINTS <<EOF
Gamma-point only
1
reciprocal
0.0  0.0  0.0  1.0
EOF

# Generate INCAR File for Initial Run
cat > INCAR <<EOF
SYSTEM = Autoflow Gas Phase

# Electronic
PREC   = Accurate
ENCUT  = 500
EDIFF  = 1E-6
ISMEAR = 0
SIGMA  = 0.05
ISPIN  = 1
ISYM   = 0
ALGO   = Normal
LREAL  = .FALSE.

# Ionic relaxation
IBRION = 2
NSW    = 50
ISIF   = 0
EDIFFG = -0.01
IVDW   = 12

# Output
LWAVE  = .FALSE.
LCHARG = .FALSE.
EOF

# Generate submission script in case needed
SLURM_SUBMIT_DIR="" # Initialize variable to not trip -u
cat > submit.sh << EOF
#!/bin/bash

#SBATCH -p haswell
#SBATCH -t 30:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=30
#SBATCH -J vasp
#SBATCH --qos=nogpu
#SBATCH --job-name="Gas"

#cd "$SLURM_SUBMIT_DIR"

source /etc/profile.d/zlmod.sh
module load arch/haswell24v2
module load intel-oneapi-mkl/2024.2.2
module load intel-oneapi-mpi/2021.12.1
module load intel/2025.0.0
module load vasp/5.4.4.pl2

ulimit -s unlimited

mpirun vasp_std > out

exit

EOF

# Return to working directory
cd ../


# ---------------------------------------------- 
## VASP files for surface slab
# Acquire slab species POTCAR
cp -f ~/../../share/Apps/vasp/5.4.4.pl2/potentials/PBE.54/$slabspec/POTCAR ./

# Prepare slab folder
mkdir ./slab
mv ./POTCAR ./slab
cd ./slab

# Create slab with ASE
cat > slab_gen.py << EOF
from ase.build import bulk, surface
from ase.io import write
from ase.constraints import FixAtoms
from ase.neighborlist import NeighborList
from scipy.spatial import Delaunay
from itertools import combinations
import numpy as np
import spglib
import json

# User-supplied variables
element = "__ELEMENT__"
packing = "__PACKING__"
miller = (__MILLER__)
target_thickness = 6.0  # Å
min_lateral_size = 8.0  # Å
layers = 1
vacuum = 15.0
a = __LATTCONST__  # Can be None
freeze_fraction = 0.4  # freeze bottom 40%

# Build bulk and slab
if a is None:
    bulk_met = bulk(element, packing, cubic=True)
else:
    bulk_met = bulk(element, packing, a=a, cubic=True)

# Dynamically determine number of layers
layers = 1

while True:
    slab = surface(bulk_met, miller, layers=layers)
    slab.center(vacuum=vacuum, axis=2)

    positions = slab.get_positions()
    cell = slab.get_cell()

    # Surface normal
    normal = np.cross(cell[0], cell[1])
    normal /= np.linalg.norm(normal)

    heights = positions @ normal
    thickness = heights.max() - heights.min()

    if thickness >= target_thickness:
        break

    layers += 1

# Dynamically determine lateral slab repetition
cell = slab.get_cell()
lx = np.linalg.norm(cell[0])
ly = np.linalg.norm(cell[1])

rx = max(1, int(np.ceil(min_lateral_size / lx)))
ry = max(1, int(np.ceil(min_lateral_size / ly)))

slab = slab.repeat((rx, ry, 1))
slab.set_pbc((True, True, True))

positions = slab.get_positions()
cell = slab.get_cell()
natoms = len(slab)

# Projection of positions to surface normal
heights = positions @ normal
max_height = np.max(heights)

layer_tol = 0.3
sorted_heights = np.sort(heights)
height_layers = []
for h in sorted_heights:
    if not height_layers:
        height_layers.append([h])
    elif abs(h - height_layers[-1][0]) < layer_tol:
        height_layers[-1].append(h)
    else:
        height_layers.append([h])

layer_means = np.array([np.mean(l) for l in height_layers])
layer_means.sort()

# Freeze bottom 2 layers to simulate bulk
sorted_layer_means = np.sort(layer_means)
n_layers = len(sorted_layer_means)
n_freeze = int(np.ceil(freeze_fraction * n_layers))

bottom_cut = sorted_layer_means[n_freeze - 1]
frozen_mask = heights < (bottom_cut + 1e-3) #1e-3 tolerance
slab.set_constraint(FixAtoms(mask=frozen_mask))


# Construct neighbor list
bulk_positions = bulk_met.get_positions()
bulk_cell = bulk_met.get_cell()

bulk_natoms = len(bulk_met)
bulk_dists = []

for i in range(bulk_natoms):
    for j in range(i+1, bulk_natoms):
        d = bulk_met.get_distance(i, j, mic=True)
        bulk_dists.append(d)

nn_dist = np.min(bulk_dists)
cutoff = 1.25 * nn_dist

nl = NeighborList([cutoff] * natoms, self_interaction=False, bothways=True)
nl.update(slab)

# Define undercoordinated sites as surface sites
coordination = np.array([len(nl.get_neighbors(i)[0]) for i in range(natoms)])
# Undercoordination defined relative to 90th percentile bulk coordination
bulk_coord = int(np.percentile(coordination, 90))

surface_atoms_all = [
    i for i in range(natoms)
    if coordination[i] < bulk_coord
    and abs(heights[i] - max_height) < 2.0
]


# Symmetry pruning preparation
lattice = slab.get_cell()
positions_frac = slab.get_scaled_positions()
numbers = slab.get_atomic_numbers()
cell_tuple = (lattice, positions_frac, numbers)
sym_data = spglib.get_symmetry_dataset(cell_tuple, symprec=1e-3)
equiv = sym_data.equivalent_atoms


# Top-layer unique representatives
unique_top = {}
for i in surface_atoms_all:
    cls = equiv[i]
    if cls not in unique_top:
        unique_top[cls] = i
pruned_surface_atoms = list(unique_top.values())

# PBC-aware helper for XY distance
def crosses_pbc(indices, frac_coords, threshold=0.8):
    if isinstance(indices, int):
        return False  # single atom cannot cross PBC
    
    coords = frac_coords[np.array(indices), :2]
    for dim in range(2):  # X and Y
        delta = coords[:, dim][:, None] - coords[:, dim][None, :]
        delta -= np.round(delta)  # wrap into [-0.5, 0.5]
        if np.max(np.abs(delta)) > threshold:
            return True
    return False


# Define a reasonable bridge cutoff based on slab geometry
all_distances = []
for i in range(natoms):
    for j in nl.get_neighbors(i)[0]:
        if i < j:
            d = np.linalg.norm(positions[i] - positions[j])
            all_distances.append(d)

bridge_cutoff = 1.2 * np.min(all_distances)  # slightly larger than shortest NN

# Bridge sites (PBC-aware + symmetry pruning)
bridge_sites = []
# Constructed bridge must also contain pruned top site
for i in surface_atoms_all:
    neighbors, _ = nl.get_neighbors(i)
    for j in neighbors:
        if j in surface_atoms_all and i < j:
            if not crosses_pbc([i,j],positions_frac,threshold=0.5):
                dist = np.linalg.norm(positions[i] - positions[j])
                if dist <= bridge_cutoff:
                    bridge_sites.append((i,j))


# Symmetry + PBC pruning
bridge_sites_pruned = []
seen_keys = set()
for i,j in bridge_sites:
    key = tuple(sorted([equiv[i], equiv[j]]))
    if key not in seen_keys:
        seen_keys.add(key)
        bridge_sites_pruned.append((i,j))


# Threefold candidates via Delaunay triangulation
# Fourflod candidates via four cycle graph search
filtered_threefolds = []
threefold_types = []
max_edge = 1.15 * np.min(all_distances)  # slightly larger than shortest NN

# Local projection function
def project_onto_plane(points, normal):
    return points - np.outer(points @ normal, normal)

# Neighbor connectivity check
def is_connected_triangle(tri_atoms, nl):
    for i in tri_atoms:
        neighbors = set(nl.get_neighbors(i)[0])
        if len(neighbors & set(tri_atoms)) < 2:
            return False
    return True

# Maximum triangle edge length check
def triangle_max_edge(tri_atoms, positions, max_edge=None):
    coords = positions[list(tri_atoms)]
    dists = [
        np.linalg.norm(coords[i] - coords[j])
        for i in range(3) for j in range(i+1, 3)
    ]
    if max_edge is None:
        return True  # no cutoff
    return max(dists) <= max_edge

# Classify between hcp or fcc threefolds
def classify_fcc_hcp(tri_atoms, positions, sub_positions, xy_tol=0.25):
    tri_xy = positions[list(tri_atoms), :2]
    centroid_xy = tri_xy.mean(axis=0)
    sub_xy = sub_positions[:, :2]
    distances = np.linalg.norm(sub_xy - centroid_xy, axis=1)
    return "hcp" if np.min(distances) < xy_tol else "fcc"

# Cluster threefold sites
def cluster_sites(candidates, positions, tol=0.5):
    clustered = []
    for cand in candidates:
        coords = positions[list(cand), :2]
        if not clustered:
            clustered.append(cand)
        else:
            dists = [np.linalg.norm(coords.mean(axis=0) - positions[list(c), :2].mean(axis=0)) for c in clustered]
            if all(dd > tol for dd in dists):
                clustered.append(cand)
    return clustered

# Fourfold terrace adjacency
def build_terrace_adjacency(terrace):

    terrace_set = set(terrace)
    coords = positions[terrace]

    # Project to XY plane (since terrace already planar)
    coords_2d = coords[:, :2]

    # Compute all pair distances
    dists = []
    for i in range(len(coords_2d)):
        for j in range(i+1, len(coords_2d)):
            d = np.linalg.norm(coords_2d[i] - coords_2d[j])
            dists.append(d)

    dists = np.sort(dists)

    # smallest two unique distances
    unique = np.unique(np.round(dists, 3))
    if len(unique) < 2:
        return {}

    short = unique[0]
    long = unique[1]

    edge_tol = 1.2 * long

    adjacency = {}

    for idx_i, i in enumerate(terrace):
        adjacency[i] = []
        for idx_j, j in enumerate(terrace):
            if j == i:
                continue
            d = np.linalg.norm(coords_2d[idx_i] - coords_2d[idx_j])
            if d <= edge_tol:
                adjacency[i].append(j)

    return adjacency


# Fourfold rectangle check
def is_rectangular_ring(quad, angle_tol=15):
    coords = positions[list(quad)]

    # Ensure same height layer
    zs = [heights[i] for i in quad]
    if max(zs) - min(zs) > layer_tol:
        return False

    # Cyclic ordering via centroid-angle
    center = coords.mean(axis=0)
    vecs = coords - center
    angles = np.arctan2(vecs[:, 1], vecs[:, 0])
    order = np.argsort(angles)
    coords = coords[order]

    edges = [coords[(i+1) % 4] - coords[i] for i in range(4)]

    for i in range(4):
        v1 = edges[i] / np.linalg.norm(edges[i])
        v2 = edges[(i+1) % 4] / np.linalg.norm(edges[(i+1) % 4])
        ang = np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1, 1)))

        if abs(ang - 90) > angle_tol:
            return False

    return True


# Fourfold inner atom check
from matplotlib.path import Path

def contains_internal_atom(quad, terrace):

    coords2d = positions[list(quad)][:, :2]

    # order cyclically
    center = coords2d.mean(axis=0)
    angles = np.arctan2(coords2d[:,1]-center[1],
                        coords2d[:,0]-center[0])
    order = np.argsort(angles)
    polygon = coords2d[order]

    path = Path(polygon)

    quad_set = set(quad)

    for atom in terrace:
        if atom in quad_set:
            continue

        point = positions[atom][:2]

        if path.contains_point(point):
            return True

    return False


# Top-layer atoms with bulk-like coordination
top_layer_idx = [i for i in range(natoms) if heights[i] > (layer_means[-1] - layer_tol)]
bulk_top_coord = int(np.median([coordination[i] for i in top_layer_idx]))
valid_top_atoms = [i for i in top_layer_idx if coordination[i] == bulk_top_coord]

top_positions = positions[valid_top_atoms]
top_frac = positions_frac[valid_top_atoms]

# Group top atoms by terrace type
# Initial grouping along global normal
terrace_groups = []
for h, idx in sorted(zip(heights[valid_top_atoms], valid_top_atoms)):
    if not terrace_groups:
        terrace_groups.append([idx])
    elif abs(h - heights[terrace_groups[-1][0]]) < layer_tol:
        terrace_groups[-1].append(idx)
    else:
        terrace_groups.append([idx])

# Local normal and Delaunay triangulation per terrace
threefold_candidates = []
for terrace in terrace_groups:
    terrace_coords = positions[terrace]
    terrace_frac = positions_frac[terrace]

    # Compute local normal via PCA-like fit
    centroid = terrace_coords.mean(axis=0)
    cov = (terrace_coords - centroid).T @ (terrace_coords - centroid)
    eigvals, eigvecs = np.linalg.eigh(cov)
    local_normal = eigvecs[:, np.argmin(eigvals)]
    local_normal /= np.linalg.norm(local_normal)

    # Project to plane
    proj_coords = project_onto_plane(terrace_coords, local_normal)
    coords_2d = proj_coords[:, :2]

    # 2D Delaunay
    try:
        tri = Delaunay(coords_2d)
    except Exception:
        continue  # Skip terrace if triangulation fails

    # Map simplices to global atom indices
    simplices = [tuple(terrace[v] for v in simplex) for simplex in tri.simplices]
    threefold_candidates.extend(simplices)

# Fourfold search per terrace
fourfold_candidates = []

for terrace in terrace_groups:

    adjacency = build_terrace_adjacency(terrace)

    for i in terrace:
        for j in adjacency[i]:
            if j <= i:
                continue
            for k in adjacency[j]:
                if k in (i, j) or k <= i:
                    continue
                for l in adjacency[k]:
                    if l in (i, j, k) or l <= i:
                        continue
                    if i in adjacency[l]:

                        quad = tuple(sorted([i, j, k, l]))

                        if crosses_pbc(quad, positions_frac, threshold=0.8):
                            continue

                        if not is_rectangular_ring(quad):
                            continue

                        if contains_internal_atom(quad, terrace):
                            continue

                        fourfold_candidates.append(quad)


# Filter by connectivity and PBC
filtered_threefolds = [
    tri_atoms
    for tri_atoms in threefold_candidates
    if is_connected_triangle(tri_atoms, nl)
        and not crosses_pbc(tri_atoms, positions_frac, threshold=0.8)
        and triangle_max_edge(tri_atoms, positions, max_edge)
]

# Subsurface probe for FCC/HCP classification (XY-based)
sub_layer_idx = [i for i in range(natoms) if heights[i] > (layer_means[-2] - layer_tol) and heights[i] < (layer_means[-1] - layer_tol)]
sub_positions = positions[sub_layer_idx]

threefold_types = [classify_fcc_hcp(tri, positions, sub_positions) for tri in filtered_threefolds]

# Apply clustering
threefold_sites = cluster_sites(filtered_threefolds, positions)
threefold_types = [threefold_types[filtered_threefolds.index(t)] for t in threefold_sites]

# Threefold symmetry pruning
fcc_candidates = []
hcp_candidates = []

for tri, ttype in zip(threefold_sites, threefold_types):
    key = tuple(sorted([equiv[i] for i in tri]))
    if ttype == "fcc":
        fcc_candidates.append((key, tri))
    else:
        hcp_candidates.append((key, tri))

def prune_unique(candidates):
    unique_dict = {}
    for key, tri in candidates:
        if key not in unique_dict:
            unique_dict[key] = tri
    return list(unique_dict.values())

fcc_sites_pruned = prune_unique(fcc_candidates)[:1]
hcp_sites_pruned = prune_unique(hcp_candidates)[:1]

# Remove fourfold duplicates
fourfold_candidates = list(set(fourfold_candidates))

# Fourfold symmetry pruning
fourfold_unique = {}
for quad in fourfold_candidates:
    key = tuple(sorted([equiv[i] for i in quad]))
    if key not in fourfold_unique:
        fourfold_unique[key] = quad

fourfold_sites_pruned = list(fourfold_unique.values())[:1]


# Write POSCAR
write("POSCAR", slab, vasp5=True, direct=True)

# Export metadata
metadata = {
    "element": element,
    "packing": packing,
    "miller": [int(x) for x in miller],
    "requested_layers": int(layers),
    "detected_layer_count": int(len(layer_means)),
    "cell": slab.get_cell().tolist(),
    "total_atoms": int(natoms),
    "surface_atoms": [int(i) for i in pruned_surface_atoms],
    "bridge_sites": [[int(i), int(j)] for i,j in bridge_sites_pruned],
    "threefold_sites": [[int(i) for i in trip] for trip in hcp_sites_pruned + fcc_sites_pruned],
    "threefold_types": ["hcp"]*len(hcp_sites_pruned) + ["fcc"]*len(fcc_sites_pruned),
    "fourfold_sites": [[int(i) for i in quad] for quad in fourfold_sites_pruned],
    "fourfold_types": ["fourfold"]*len(fourfold_sites_pruned),
    "frozen_atoms": [int(i) for i in np.where(frozen_mask)[0]]
}

with open("surface_atoms.json", "w") as f:
    json.dump(metadata, f, indent=2)
EOF

# Safely inject variable information
sed -i "s|__ELEMENT__|$slabspec|" slab_gen.py
sed -i "s|__PACKING__|$packing|" slab_gen.py
sed -i "s|__MILLER__|$miller|" slab_gen.py
sed -i "s|__LATTCONST__|$lattconst|" slab_gen.py

# Run Python code
python slab_gen.py || {
    echo "Error: slab-phase POSCAR generation failed"
    exit 1
}
rm -f slab_gen.py
mv ./surface_atoms.json ../

# Prepare KPOINTS, INCAR, SLURM submission script
# Generate SLAB KPOINTS File
cat > KPOINTS <<EOF
Slab k-points
0
Gamma
6 6 1
0 0 0
EOF

# Generate INCAR File for Initial Run
cat > INCAR << EOF
SYSTEM = $slabspec slab

# Main
PREC   = Accurate
ENCUT  = 500
EDIFF  = 1E-6
LREAL  = .FALSE.
ALGO   = Normal
ISMEAR = 0
SIGMA  = 0.05
IBRION = 2
NSW    = 150
ISIF   = 2
EDIFFG = -0.02
ISYM   = 0
ISPIN  = 1   # set to 2 only if magnetic
IVDW   = 12
LDIPOL = .TRUE.
IDIPOL = 3

# Output
LWAVE  = .FALSE.
LCHARG = .FALSE.
EOF

# Generate submission script in case needed
cat > submit.sh << EOF
#!/bin/bash

#SBATCH -p haswell
#SBATCH -t 30:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=30
#SBATCH -J vasp
#SBATCH --qos=nogpu
#SBATCH --job-name="Slab"

#cd $SLURM_SUBMIT_DIR

source /etc/profile.d/zlmod.sh
module load arch/haswell24v2
module load intel-oneapi-mkl/2024.2.2
module load intel-oneapi-mpi/2021.12.1
module load intel/2025.0.0
module load vasp/5.4.4.pl2

ulimit -s unlimited

mpirun vasp_std > out

exit

EOF

# Return to working directory
cd ../


# ---------------------------------------------- 
## Adsorption mode enumeration
# Prepare DockonSurf input writer in Python instead of bash
# Streamlined parsing of json metadata
cat > input_gen.py << EOF
import json
from ase.io import read
import numpy as np

# Process slab metadata
with open("surface_atoms.json") as f:
    meta = json.load(f)

# Format pbc cell
cell_str = " ".join(
    f"({v[0]:.6f} {v[1]:.6f} {v[2]:.6f})"
    for v in meta["cell"]
)

# Collect sites
sites = []

top_sites = meta.get("surface_atoms", [])
bridge_sites = meta.get("bridge_sites", [])
threefold_sites = meta.get("threefold_sites", [])

sites += [str(i) for i in top_sites]
sites += [f"({i},{j})" for i, j in bridge_sites]
sites += [f"({i},{j},{k})" for i, j, k in threefold_sites]

sites_str = ", ".join(sites)

# Process gas metadata
gas_atoms = read("gas/meta.xyz")
print("Available arrays:", gas_atoms.arrays.keys())

tags = gas_atoms.get_tags()

unique_ids = sorted(set(tags) - {-1})

molec_sites = []

if len(unique_ids) == 0:
    # Robust fallback: use atom index 0
    molec_sites = ["0"]

else:
    for uid in unique_ids:
        group = np.where(tags == uid)[0].tolist()

        if len(group) == 1:
            molec_sites.append(str(group[0]))
        else:
            molec_sites.append(f"({','.join(map(str, group))})")

molec_ctrs_str = ", ".join(molec_sites)

with open("dockonsurf.inp", "w") as f:
    f.write(f"""
[Global]
run_type = Screening
code = VASP
batch_q_sys = False
project_name = {meta["element"]}_enum
pbc_cell = {cell_str}

[Screening]
screen_inp_file = slab/INCAR slab/KPOINTS POTCAR
surf_file = slab/POSCAR
use_molec_file = gas/POSCAR
molec_ctrs = {molec_ctrs_str}
sites = {sites_str}
adsorption_height = 2.0
set_angles = euler
sample_points_per_angle = 2
surf_normal_vect = z
max_structures = False
""")
EOF

# Concatenate slab and gas POTCAR
cat ./slab/POTCAR ./gas/POTCAR > POTCAR

# Generate DockonSurf input file & run enumeration
python input_gen.py || {
    echo "Error: DockonSurf input generation failed"
    exit 1
}
dockonsurf.py -i dockonsurf.inp -f # DockonSurf run on foreground for clarity

# File cleanup
rm -f input_gen.py POTCAR
mv dockonsurf.inp dockonsurf.log ./screening



# ---------------------------------------------- 
## Coarse optimization + filtering
# Prepare hybrid optimization script
cd ./screening

cat > hybrid_screen.py << 'EOF'
import sys, json, os
import numpy as np
import gc
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import FIRE, LBFGS
from tblite.ase import TBLite
from xtb.ase.calculator import XTB
from mace.calculators import mace_mp
from chgnet.model.dynamics import CHGNetCalculator

# Main optimization function
def run_periodic_relaxation(atoms_in, calculator, label):
    atoms = atoms_in.copy()
    atoms.set_pbc(True)
    atoms.calc = calculator

    opt = LBFGS(atoms, logfile=f"{label}.log")

    try:
        opt.run(fmax=0.2, steps=100)
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        max_force = float(np.linalg.norm(forces, axis=1).max())

        # Store energy inside atoms.info
        atoms.info["method"] = label
        atoms.info["energy"] = float(energy)
        atoms.info["max_force"] = max_force

        # Remove calculator object before writing
        atoms.calc = None

        write(f"relaxed_{label}.xyz", atoms)

        return {
            "method": label,
            "energy": float(energy),
            "max_force": max_force
        }

    except Exception as e:
        with open(f"{label}_FAILED", "w") as f:
            f.write(str(e))
        return None



# Read structure
atoms = read("POSCAR")

# Freeze slab for cheaper optimization
with open(os.path.join("..", "..", "surface_atoms.json")) as f:
    slab_info = json.load(f)
n_slab_atoms = slab_info["total_atoms"]

mask = np.ones(len(atoms), dtype=bool)
mask[:n_slab_atoms] = False  # False -> frozen
atoms.set_constraint(FixAtoms(mask=~mask))  # Fix slab atoms

# Hybrid optimization
# Prerun with GFN-FF and FIRE for cleaning
atoms.set_pbc(False)
atoms.calc = XTB(method="GFN-FF", charge=0, spin=0, maxiter=250, electronic_temperature=3000)
opt = FIRE(atoms, logfile="prerun.log")
opt.run(fmax=2, steps=15)
atoms_after_gfnff = atoms.copy()
atoms_after_gfnff.set_pbc(True)

# Ensemble optimization
results = []

# GFN1-xTB
xtb_calc = TBLite(method="GFN1-xTB", charge=0, spin=0)
results.append(
    run_periodic_relaxation(atoms_after_gfnff, xtb_calc, "GFN1-xTB")
)
del xtb_calc
gc.collect()

# MACE
mace_calc = mace_mp(model="medium", device="cpu")
results.append(
    run_periodic_relaxation(atoms_after_gfnff, mace_calc, "MACE-MP")
)
del mace_calc
gc.collect()

# CHGNet
chgnet_calc = CHGNetCalculator()
results.append(
    run_periodic_relaxation(atoms_after_gfnff, chgnet_calc, "CHGNet")
)
del chgnet_calc
gc.collect()

with open("ensemble_screen.json", "w") as f:
    json.dump([r for r in results if r is not None], f, indent=2)
EOF

# Loop over all enumerated adsorption modes
i=0
for d in conf_*; do
  echo "Hybrid Screening Starting $d"

  (
    set +e
    cd "$d" || exit 1
    python ../hybrid_screen.py > central.log 2>&1 || { [ -f FAILED ] || echo "Screening failed (bash)" > FAILED; }
    [[ -f gfnff_topo ]] && rm gfnff_topo
  ) &

  ((++i))
  if (( i % MAXPARJOB == 0 )); then
    wait || true
  fi
done

wait || true
echo "Hybrid Screening Complete"
