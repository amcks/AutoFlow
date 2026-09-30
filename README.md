# AutoFlow v.1.4

**AutoFlow** is an automated Python pipeline for initial adsorbate structure generation, site enumeration, and multi-tier MLIP screening on metallic surface slabs.

It orchestrates the workflow from gas-phase SMILES processing and surface slab construction down to parallel MLIP relaxations (GFN-FF / MACE), structure clustering, and candidate selection for single-point DFT calculations.

---

## Features

- **Automated Structure Generation**: Generates gas-phase adsorbate geometries from SMILES and metallic surface slabs with defined Miller indices.
- **Site Enumeration**: Direct placement for monoatomic species and seamless integration with DockOnSurf for polyatomic molecules.
- **Parallel MLIP Relaxation**: Multi-threaded MLIP structure relaxations using FIRE prerelaxation and LBFGS main optimization.
- **Ensemble Post-Analysis**: Hierarchical clustering, reactivity/adsorption filtering, and geometric disagreement analysis.
- **DFT Dataset Preparation**: Automatically exports top representative candidates into VASP-ready calculation directories.

---

## Installation & Prerequisites

### 1. Prerequisites
On top of the dependencies installed via pip, AutoFlow also requires `dockonsurf.py` for site enumeration.

```bash
# Create base Conda environment
conda create -n autoflow python=3.10
conda activate autoflow
```

**Note**: [DockOnSurf](https://gitlab.com/lch_interfaces/dockonsurf) must be installed and accessible in your system `$PATH`.

### 2. Installing AutoFlow
After **fulfilling the prerequisites** and **activating the environment**, clone the repository and install with `pip` in editable mode:

```bash
git clone https://github.com/amcks/AutoFlow
cd AutoFlow
pip install -e .
```

## Environment Setup
Configure paths for VASP pseudopotentials and MACE model via environment variables in your `./bashrc` or within your SLURM submission scripts:

```bash
export VASP_PP_PATH="/path/to/vasp/potentials"
export MACE_MODEL_PATH="/path/to/MACE/.model/file"
```

## Usage
AutoFlow provides a command-line interface driven by subcommands.

To view options and available commands, include the `--help` option at either the global or subcommand level:
```bash
autoflow --help
autoflow run --help
```

### 1.Running Screening Pipeline
Execute adsorption mode enumeration, structure generation, screening, and post-analysis:
```bash
autoflow run -s Ag -m 1,1,1 -a "C(=O)C" -l 4.13 -p fcc -j 4
```

To automatically generate VASP DFT inputs of the representative structure(s) for further construction of datasets upon completion of the MLIP screening:
```bash
autoflow run -s Ag -m 1,1,1 -a "C(=O)C" --generate-dft
```

**Available Options**:
- `-s, --slab`: Surface slab element (e.g., Ag, Cu, Pt).
- `-m, --miller`: Comma-separated Miller indices (e.g., 1,1,1).
- `--poscar-slab`: Path to optional slab POSCAR file.
- `--poscar-gas` : Path to optional gas POSCAR file.
- `--site-slab` : Comma-separated atomic indices (e.g. 43,5,10,22) for surface site override in supplied slab POSCAR file.
- `--site-gas` : Comma-separated atomic indices for adsorbate molecule anchor points override in supplied gas POSCAR file.
- `-a, --adsorbate`: Adsorbate SMILES string (e.g., C(=O)C).
- `-l, --lattconst`: Optional lattice constant (Å).
- `-p, --packing`: Crystal structure (fcc, hcp, bcc, bct). Default: fcc.
- `-j, --jobs`: Maximum parallel screening processes. Default: 4.
- `--vasp-pp` : Path to  VASP POTCAR directory if not specified via bash variable.
- `--mace-model` : Path to MACE model file if not specified via bash variable.
- `--no-screen` : \[Flag\] Option to stop after adsorption mode enumeration without performing MLIP screening.
- `--generate-dft`: \[Flag\] Option to prepare DFT input directories in screening/DFT/ after screening.

**Notes**:
- Usage of `-s` & `-m` options are mutually exclusive with the `--poscar-slab` option.
- Usage of the `-a` option is mutually exclusive with the `--poscar-gas` option. Additionally, `--site-gas` must be specified when `--poscar-gas` is in use.


### 2.Standalone DFT Preparation
If MLIP screening was completed previously, VASP DFT input data can be generated separately:
```bash
autoflow prep-dft -d ./screening
```

## Output Directory Structure
Executing `autoflow run` creates the following directory layout:
```
.
├── gas/                     # Gas-phase POSCAR and meta.xyz
├── slab/                    # Clean surface POSCAR
└── screening/               # Enumerated configurations (conf_0, conf_1, ...)
    ├── conf_0               # Relaxed configuration structures
    ├── conf_1 
    ├── surface_atoms.json   # Surface site metadata
    ├── ensemble_screening_summary.json
    ├── dft_selection.json
    └── DFT/                 # Generated if --generate-dft or prep-dft is called
        ├── 0/               # VASP input set (POSCAR, INCAR, KPOINTS, POTCAR, submit.sh)
        └── 1/
```
