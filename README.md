# Automated Design of NHC ligands for short Ru=C bond
This repository collects input files used to automatically design N-heterocyclic carbene (NHC) ligands (L) 
meant to generate ruthenium olefin metathesis catalysts with formula (L)Ru(Cl)(Cl)=CH<sub>2</sub> bearing a particularly short Ru=CH<sub>2</sub> bond.

## Prerequisites and Preparation
The experiments were run on an HPC Linux cluster, but amost tasks, i.e., all excluding the highly computational demamding DFT computations can be reproduced on a regular laptop/workstation. The following software in needed to reproduce the experiments:
* <a href="https://github.com/denoptim-project/DENOPTIM">DENOPTIM 2.2.6</a> to run the evolutionary design and fragment space enumeration.
* <a href="https://dasher.wustl.edu/tinker/">Tinker v6.3</a> for performing molecular mechanic calculations like conformational searches. **NB:** Tinker must be compiled with `maxval=12` in Tinker's source file `sizes.i`.
* <a href="https://gaussian.com/">Gaussian 09</a> for running DFT calculations. We expect to run these experiments on a HPC cluster managed by a scheduler (either PBS or SLURM) and Gaussian jobs are submitted to the queue. Submission of the job is heavily machine- and user-specific so we will assume that a job submission BASH script exists and can be called by the [fitness provider script](evolutionary_desing/Ru_14-el_fitness_BndLng.sh).
* <a href="http://openbabel.org/wiki/Main_Page">OpenBabel</a> for chemical file format conversion.
* <a href="">AutoCompChem</a> for performaning molecular modeling tasks automatically.
* Java JDK-8 for running DENOPTIM and AutoCompChem

The file [environment.yml](environment.yml) can be used to recreate a minor part of the environment, namely, install OpenBabel and Java JDK. This can be done with [conda](https://docs.conda.io/en/latest/index.html):
```
conda env create --file environment.yml
```

Next, DENOPTIM and AutoCompChem are both nested into the present repository as submodules. This means that the right version of these tools can be obtained by adding `--recurse-submodules` when cloning the present repository (`git clone --recurse-submodules <URL>`)

Finally, Tinker and Gaussian must be obtained and installed manually. Refer to the respective links above for instructions and license terms.

## Looking for a sandbox dry-run?
You might want to play with the data without engaging in computationally demanding DFT calculations. For this purpose, a super light fitness provider can be chosen and used in evolutionary design. The fitness function will rank the candidates based solely on their molecular weight: Obviously, this does not lead to any chemically meaningful design. Still, it allows you to play with the data.
To run this kind of experiments, you can change the fitness provider configuration in the [evolutionary_desing/evolution.params file](evolutionary_desing/evolution.params) by uncommenting the "Fake provider" and commenting the "Actual provider".
