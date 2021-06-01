# Results of Evolutionary Experiments
In this folder you find the collections of candidates evaluated in each experiments. Each `run_#.#` folder contains the data pertaining to one single experiment. In particular, to limit the size of these folders, each candidate is represented by a `*FIT.sdf` file and a `*_DFT.xyz` file.
The `*FIT.sdf` contains the various representations of the candidate, including InChIKey, SMILES, DENOPTIM's graph representation, and an un-unrefined geometry meant only to facilitate visual identification of the molecular structure.
The `*_DFT.xyz` contains the DFT-refined geometry used to calculate the geometrical descriptors, evaluate constraints, and produce the numerical fitness value reported in the `*FIT.sdf` file.


## Note About Pathnames
Note that the these data has been moved from the parent folder into this subfolder. Therefore, relative pathnames should refer to the original location of the folder.
