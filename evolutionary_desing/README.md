# Evolutionary Design Data
This folder contains input files and the results from evolutionary design experiments. 

In particular, the [results](results) folder collects the actual results from several experiments, including the DFT-optimized geometries and the characterization of each candidate ligand visited by a specific run. Foe example, a pathname like `results/run_*/Gen*/*_FIT.sdf` provides the detailed definition of a single candidate, possibly associated with a fitness value (see sdf tag `FITNESS`), while `results/run_*/Gen*/*_DFT.xyz` provides the DFT-optimized geometry used to calculate such fitness. 
All the analysis of the data and the plots reported in the published paper are availanle under [results/analysis](results/analysis/).

In addition, the [test](test) folder contains scripts and data to perform test runs.
