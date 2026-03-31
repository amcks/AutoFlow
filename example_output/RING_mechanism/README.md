## Example: Mechanism Sampling
This subdirectory contains an example usage of AutoFlow integrated into the Rule Input Network Generator \([RING](https://doi.org/10.1016/j.compchemeng.2012.06.008)\) workflow. For this example, the system of acetaldehyde combustion on Ag(111) is explored. RING was initially used to generate a pool of relevant species and intermediates, as well as a list of 5,864 reaction mechanisms, with lengths ranging from 6 to 8 elementary steps. It should be noted that RING allows for the use of dummy atoms `[{M}]` to represent sites at which the intermediate is bonded to the catalytic surface, which is first converted to the standardized SMILES dummy atom `[*]` before processing. The former is seen in the native RING species file \(`species.txt`\), while the latter is seen in the pre-processed species file \(`cleaned_smiles.txt`\).

The usage is separated into two stages:
- Energy Mapping: MLP-based AutoFlow usage for species pool
- Mechanism Selection: Combining generated list of mechanisms with energy map to find most thermodynamically feasible reaction pathway

## Energy Mapping
A bash script \(`mass_submit.sh`\) reads the list of relevant species in \(`cleaned_smiles.txt`\), then generates the full set of subdirectories and performs the generation, enumeration, and optimization of each intermediate via AutoFlow. For this example, optimization runs using MACE-MP and CHGNet were both performed on all the structures. For the sake of brevity, instead of including the full optimization result of all 341 sampled intermediate structures, only the example standard output \(`std.out`\) file for one of the structures is included here.

At the end of the optimization, the Python script \(`energy_parse.py`\) can be executed by specifying the method of choice, generating a .csv file which serves as the energy map for the next step.

## Mechanism Selection
With the energy map, the list of mechanisms and their elementary steps contained in `CO2pathways.txt` is read, parsed, and matching with the generated energy mapping from the previous step allows for evaluation of the reaction energies of each step, and therefore, an assessment of which step would be the most thermodynamically feasible.

Considering the relatively sizeable list of reaction mechanisms, an early-exit pruning strategy is adopted. Formulated as a min-max problem, the goal of this step is to find the reaction mechanism in which its most thermodynamically demanding step is the lowest out of any step across the different mechanisms. In this approach, the most endothermic step in a mechanism is recorded, and when evaluating other steps in other mechanisms, encountering a step with an even more endothermic step allows for an early exit to save computing resources, whereas encountering a different mechanism whose most endothermic step is lower will instead cause the code to update the current best mechanism.

As of now, this combination of enumeration and early exit strategy is sufficient, but for future scaling considerations, a graph-based structure combined with a Dijkstra,etc. exploration mechanism can be employed instead for a more efficeint handling of larger reaction networks.
