import json
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase import Atoms

def surface_normal(cell):
    a, b = np.array(cell[0]), np.array(cell[1])
    n = np.cross(a, b)
    return n / np.linalg.norm(n)

def centroid(positions, indices):
    return np.array([positions[i] for i in indices]).mean(axis=0)

def generate_monoatomic_configurations(slab_path: Path, gas_path: Path, json_path: Path, output_dir: Path) -> None:
    """Handles monoatomic adsorbate placement directly via ASE."""
    slab = read(slab_path)
    adsorbate = read(gas_path)
    
    if len(adsorbate) != 1:
        raise ValueError("Monoatomic enumeration only supports adsorbates with 1 atom.")

    ads_symbol = adsorbate[0].symbol
    site_data = json.loads(json_path.read_text())
    normal = surface_normal(slab.cell)

    sites = []
    pos = slab.positions
    for i in site_data["surface_atoms"]:
        sites.append(("top", pos[i]))
    for pair in site_data.get("bridge_sites", []):
        sites.append(("bridge", centroid(pos, pair)))
    for tri in site_data.get("threefold_sites", []):
        sites.append(("threefold", centroid(pos, tri)))
    for quad in site_data.get("fourfold_sites", []):
        sites.append(("fourfold", centroid(pos, quad)))

    heights = {"top": 2.0, "bridge": 1.5, "threefold": 1.5, "fourfold": 1.5}

    for idx, (site_type, coord) in enumerate(sites):
        height = heights.get(site_type, 1.1)
        ads_pos = coord + height * normal
        
        struct = slab.copy()
        struct += Atoms(ads_symbol, positions=[ads_pos])

        conf_dir = output_dir / f"conf_{idx}"
        conf_dir.mkdir(parents=True, exist_ok=True)
        write(conf_dir / "POSCAR", struct, format="vasp")

def generate_dockonsurf_input(json_path: Path, gas_meta_path: Path, output_file: Path) -> None:
    """Generates dockonsurf.inp configuration file for polyatomic molecules."""
    meta = json.loads(json_path.read_text())
    cell_str = " ".join(f"({v[0]:.6f} {v[1]:.6f} {v[2]:.6f})" for v in meta["cell"])

    sites = [str(i) for i in meta.get("surface_atoms", [])]
    sites += [f"({i},{j})" for i, j in meta.get("bridge_sites", [])]
    sites += [f"({i},{j},{k})" for i, j, k in meta.get("threefold_sites", [])]
    sites += [f"({i},{j},{k},{l})" for i, j, k, l in meta.get("fourfold_sites", [])]

    gas_atoms = read(gas_meta_path)
    tags = gas_atoms.get_tags()
    unique_ids = sorted(set(tags) - {-1})

    molec_sites = []
    if len(unique_ids) == 0:
        molec_sites = ["0"]
    else:
        for uid in unique_ids:
            group = np.where(tags == uid)[0].tolist()
            molec_sites.append(str(group[0]) if len(group) == 1 else f"({','.join(map(str, group))})")

    output_file.write_text(f"""[Global]
run_type = Screening
code = VASP
batch_q_sys = False
project_name = {meta["element"]}_enum
pbc_cell = {cell_str}

[Screening]
screen_inp_file = INCAR KPOINTS POTCAR
surf_file = slab/POSCAR
use_molec_file = gas/POSCAR
molec_ctrs = {", ".join(molec_sites)}
sites = {", ".join(sites)}
adsorption_height = 2.0
set_angles = euler
sample_points_per_angle = 2
surf_normal_vect = z
max_structures = False
""")
