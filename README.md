# AutoFlow
The AutoFlow script is intended for automating the generation of initial adsorption structures of the given adsorbate on the slab of the given element.

At the end of the operation, the script generates structure files for the gaseous molecule, the clean slab, and enumerated adsorption configurations, along with the corresponding VASP input files to be used for later optimization using DFT or otherwise.

A hybrid screening scheme based on forces is used to initially filter out unphysical solutions using a combination of Grimme's GFN-FF -> GFN1-xTB methods, as well as other machine learning potential methods like MACE-MP and CHGNet.

## Installation & Dependencies: 

The complete list of dependencies is listed in the included *AF_env.yaml* file, allowing for installation via conda.
```
conda env create -f AF_env.yml
```

At the end of which, an environment called "autoflow" will be created, which is activated automatically in the script, or can be activated manually via the following command.
```
conda activate autoflow
```

## Usage
In bash terms, the usage of this script is as follows:
```
Usage: $(basename "$0") -s SLAB -m H,K,L -a SMILES [-l LATTCONST] [-p PACKING] [-h/--help]
```

The first three options are mandatory, and the script will not execute unless they are supplied:
- **s**:  Slab element (e.g. Cu, Pt).
- **m**:  Comma-separated Miller indices (e.g. 1,1,1).
- **a**:  Adsorbate SMILES string (e.g. CH3COOH).

The remaining options are optional:
- **l**:  Lattice constant. Defaults to ASE's database of lattice constants when not specified.
- **p**:  Packing/crystal structure type. Must be one of the following: fcc, hcp, bcc, bct. Defaults to fcc when not specified.
- **h**:  Print this help message and exit. When this option is called, all other options are ignored, and only the help message is printed. Can also be called with --help.

For reference, a set of the example output of a successful run for the test system of methoxy (CH3O\*) adsorbate on an Ag(211) surface is included in the "example\_output" subdirectory.

## Features to be Implemented
Below are the remaining changes to be implemented in the coming versions.
### High Priority:
- Expand extent of user control in getopts: Calculator control, DockOnSurf parameters, Parallelization.

### Medium Priority:
- Implementation of adsorption *pattern* enumeration for higher coverage cases.
- Implementation of coarse NMA & thermodynamics-based corrections.
- Add functionality to import user-supplied surface slabs instead of having to generate a surface via ASE every time.

### Low Priority:
- Refactoring to full Python script can create considerable speedup, especially in terms of I/O or RAM handling, as well as, possibly, JIT compilation.
