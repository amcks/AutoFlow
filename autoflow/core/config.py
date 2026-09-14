import os
from pathlib import Path

# VASP Pseudopotential Base Directory
_env_vasp_pp = os.getenv("VASP_PP_PATH")
DEFAULT_VASP_POTENTIAL_PATH = (
    Path(_env_vasp_pp) if _env_vasp_pp 
    else Path("~/../../share/Apps/vasp/5.4.4.pl2/potentials/PBE.54").expanduser()
)

# Pre-trained / Fine-tuned MACE Potential Checkpoint
_env_mace_model = os.getenv("MACE_MODEL_PATH")
DEFAULT_MACE_PATH = (
    Path(_env_mace_model) if _env_mace_model 
    else Path("~/MACE_finetune/mace_adsorption_ft_cpu.model").expanduser()
)

# Atomic Number Hash Table
ELEMENT_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
    "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47,
    "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
    "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65,
    "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
    "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97,
    "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103,
    "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112,
    "Fl": 114, "Lv": 116
}

# INCAR String Templates
GAS_INCAR_TEMPLATE = """
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
"""

SLAB_INCAR_TEMPLATE = """
ENCUT  = 500
PREC   = Accurate
ISMEAR = 1
SIGMA  = 0.2
ISPIN  = 1
IBRION = -1
NSW    = 0
NELM   = 100
EDIFF  = 1E-5
LCHARG = .FALSE.
LWAVE  = .FALSE.
LREAL  = Auto
ISYM   = 0
LDIPOL = .FALSE.
!IDIPOL = 3
NCORE  = 4
"""

SLURM_SCRIPT = """
#!/bin/bash

#SBATCH -p hawkcpu,haswell
#SBATCH -t 07:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -J vasp
#SBATCH --qos=nogpu
#SBATCH --job-name="MPFT_0"

source /etc/profile.d/zlmod.sh

ulimit -s unlimited

./autoflow_v0_10_postFT.sh -s Ag -m 1,1,1 -a "C(=O)C" -l 4.13 -p fcc

exit
"""
