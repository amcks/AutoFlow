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

## GFN1-xTB (Commented out for now for speedup)
#xtb_calc = TBLite(method="GFN1-xTB", charge=0, spin=0)
#results.append(
#    run_periodic_relaxation(atoms_after_gfnff, xtb_calc, "GFN1-xTB")
#)
#del xtb_calc
#gc.collect()

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
