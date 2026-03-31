# AutoFlow
The AutoFlow script is intended for automating the generation and enumeration of initial adsorption structures of the given adsorbate on the surface slab of the given element.

The script generates structure files for the gaseous molecule, the clean slab, and enumerated adsorption configurations, along with the corresponding VASP input files to be used for later optimization using DFT or otherwise. A hybrid screening scheme based on forces is used to initially filter out unphysical solutions using a combination of GFN-FF for initialization, followed by structural optimization with GFN1-xTB, as well as other machine learning potential methods like MACE-MP and CHGNet. Post-analysis is then be performed to cluster the configurations after optimization using all the methods, and representative low-energy structures from each cluster can be selected as initial structures for subsequent, more computationally-exhaustive calculations.

Although the script is primarily written and executed in bash script, the bash script portions act as a wrapper, while the main logic of the operations is written in procedurally-generated Python script blocks. This reflects the nature of the project's evolution, which started as a simpler bash script, with plans on fully refactoring it to Python in later versions.

## Installation & Dependencies: 

The complete list of dependencies is listed in the included *AF_env.yaml* file, allowing for installation via conda.
```
conda env create -f AF_env.yml
```

At the end of which, an environment called 'autoflow' will be created, which is activated automatically in the script, or can be activated manually via the following command.
```
conda activate autoflow
```

The enumeration of the adsorptin sites relies on the usage of [DockOnSurf](https://gitlab.com/lch_interfaces/dockonsurf). Users will need to follow the installation steps outlined in their repository/documentation, but their dependencies have been covered by the provided yaml file. Note: **THIS STEP IS MANDATORY**.

## Usage
The primary script to be executed is 'autoflow\_\<version\>.sh'. In bash terms, the usage of this script is as follows:
```
Usage: autoflow_<version>.sh -s SLAB -m H,K,L -a SMILES [-l LATTCONST] [-p PACKING] [-h/--help]
```

The first three options are mandatory, and the script will not execute unless they are supplied:
- **s**:  Slab element (e.g. Cu, Pt).
- **m**:  Comma-separated Miller indices (e.g. 1,1,1).
- **a**:  Adsorbate SMILES string. Where applicable, bonds formed during adsorption should be specified by a dummy atom \[\*\] (e.g. CO\[\*\] for methoxy, C=C for ethylene, c1ccccc1 for benzene).

The remaining options are optional:
- **l**:  Lattice constant. Defaults to ASE's database of lattice constants when not specified.
- **p**:  Packing/crystal structure type. Must be one of the following: fcc, hcp, bcc, bct. Defaults to fcc when not specified.
- **h**:  Print this help message and exit. When this option is called, all other options are ignored, and only the help message is printed. Can also be called with --help.

For reference, a set of the example output of a successful run for the test system of methoxy (CH3O\*) adsorbate on an Ag(110) surface is included in the 'example\_output' subdirectory. This set of results was generated using the following command:
```
autoflow_v0_8.sh -s Ag -m 1,1,0 -a CO[*] -l 4.13 -p fcc
```

At the end of the initial optimization, a post-analysis procedure is performed to filter and cluster all the solutions from the different methods on the basis of geometry and energies. While previous versions have the post-analysis script as a separate, standalone module, in the current version, the post-analysis module is integrated into the main script, but the python script used is still printed for the sake of debugging.

## Example Applications
Currently, the 'example\_output' subdirectory contains two examples:
- AutoFlow usage on a single adsorption structure optimization
- Integration with the Rule Input Network Generator \([RING](https://doi.org/10.1016/j.compchemeng.2012.06.008)\) workflow to explore feasible mechanisms

More details are included in the `README.md` file within the 'example\_output' subdirectory.

## Features to be Implemented
Below are the remaining changes to be implemented in the coming versions.
### High Priority:
- Expand extent of user control in getopts: Calculator control, DockOnSurf parameters, Parallelization.

### Medium Priority:
- Implementation of adsorption *pattern* enumeration for higher coverage cases.
- Implementation of coarse NMA & thermodynamics-based corrections.
- Add functionality to import user-supplied surface slabs instead of having to generate a surface via ASE every time.
- Integrate post-analysis script into main script operation. 

### Low Priority:
- Refactoring to full Python script can create considerable speedup, especially in terms of I/O or RAM handling, as well as, possibly, JIT compilation.
