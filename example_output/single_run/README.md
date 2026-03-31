## Example: Single Run
This subdirectory contains the example output of a successful run for the test system of methoxy (CH3O\*) adsorbate on an Ag(110) surface. This set of results was generated using the following command:
```
autoflow_v0_8.sh -s Ag -m 1,1,0 -a CO[*] -l 4.13 -p fcc
```

At the end of the initial optimization, a post-analysis procedure is performed to filter and cluster all the solutions from the different methods on the basis of geometry and energies. While previous versions have the post-analysis script as a separate, standalone module, in the current version, the post-analysis module is integrated into the main script, but the python script used is still printed for the sake of debugging.

**NOTE**: This set of example results was generated using AutoFlow version 0.8.0, which is deprecated, but accessible through the tags function from the main github repository.
