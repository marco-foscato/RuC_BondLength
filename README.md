# Automated Design of NHC Ligands that Shorten Ru=C Bond
This repository collects input files used to automatically design N-heterocyclic carbene (NHC) ligands (L) 
meant to generate ruthenium olefin metathesis catalysts with formula (L)Ru(Cl)(Cl)=CH<sub>2</sub> bearing a particularly short Ru=CH<sub>2</sub> bond.

## Prerequisites and Preparation
### Overview
The experiments were run on an HPC Linux cluster, but most tasks, i.e., all excluding the highly computational demamding DFT computations can be reproduced on a regular laptop/workstation. The following software in needed to reproduce the experiments, and will be installed in the next section:
* <a href="https://github.com/denoptim-project/DENOPTIM">DENOPTIM 2.2.6</a> to run the evolutionary design and fragment space enumeration.
* <a href="https://dasher.wustl.edu/tinker/">Tinker v6.3</a> for performing molecular mechanic calculations like conformational searches. **NB:** Tinker must be compiled with `maxval=12` in Tinker's source file `sizes.i`.
* <a href="https://gaussian.com/">Gaussian 09</a> for running DFT calculations. We expect to run these experiments on a HPC cluster managed by a scheduler and Gaussian jobs are submitted to the queue. Submission of the job is heavily machine- and user-specific so we will assume that a job submission command exists and can be called by the [fitness provider script](evolutionary_desing/Ru_14-el_fitness_BndLng.sh).
* <a href="http://openbabel.org/wiki/Main_Page">OpenBabel</a> for chemical file format conversion.
* <a href="https://github.com/denoptim-project/AutoCompChem">AutoCompChem</a> for performaning molecular modeling tasks automatically.
* An additional utility for calculating geometrical descriptors, checking geometrical contraints, and cumputing the numerical value of the fitness (`FitnessRuCH2BndLng`).
* Java JDK-8 for running DENOPTIM, AutoCompChem, and FitnessRuCH2BndLng.

In addition, due to a limitation in the length of pathnames that Tinker can handle, you should make sure the path to the root of this folder is short (say, shorter than 30 chracters).

### Procedure
1. The file [environment.yml](environment.yml) can be used to recreate the core of the environment, namely, install OpenBabel and Java JDK. This can be done with [conda](https://docs.conda.io/en/latest/index.html):
```
conda env create --file environment.yml
``` 
After creation you should activate the environment as explained in the log from the above command.

2. The source of DENOPTIM and AutoCompChem is available in the present repository under the [tools folder](tools) as git submodule. This means that the right version of these tools can be obtained by adding `--recurse-submodules` when cloning the present repository (`git clone --recurse-submodules <URL>`). Next, you can build the executable for DENOPTIM
```
cd tools/DENOPTIM/build/
./build-all.sh
cd ../../../
```
and for AutoCompChem
```
./tools/AutoCompChem/scripts/build.sh
```

3. <a href="https://dasher.wustl.edu/tinker/">Tinker v6.3</a> and <a href="https://gaussian.com/">Gaussian 09</a> must be obtained according to the respective license terms. Please refer to the above links for license terms and installation instructions. After installation, make the Tinker excecutable reachable, by exporting an evironmental variable pointing to Tinker's `bin` forlder. For example:
```
export TINKERBIN=/path/to/your/tinker/bin
```
Also, we expect to run Gaussian via a queuing system. Therefore, we assume you have a way to submit Gaussian jobs (send jobs to the queuing system). In our case, this task is performed by the `submit_gaussian_job` command.

4. The `FitnessRuCH2BndLng` utility is compiled with ease:
```
cd tools/FitnessRuCH2BndLng
./build-fitnessruch2bndlng.sh
cd ../..
```

5. Finally, export this variable pointing to the location of the `RuC_BondLength/evolutionary_desing` subfolder:
```
export RUCBONDDESIGN=/your/path/to/RuC_BondLength/evolutionary_desing
```

## Running Evolutionary Design Experiments
After having performed the above preparation, you can run an evolutionary experiment in a machine with sufficient cpu capacity to run very many DFT geometry optimization. This typically means you want to run on an HPC cluster (but to play with the data you can do a [dry run](#dry-runs).
```
cd evolutionary_desing
java -jar ../tools/DENOPTIM/build/DenoptimGA.jar evolution.params
```

### Dry runs
You might want to play with the data without engaging in computationally demanding DFT calculations. For this purpose, a super light fitness provider can be chosen and used in evolutionary design. The fitness function will rank the candidates based solely on their molecular weight: Obviously, this does not lead to any chemically meaningful design. Still, it allows you to play with the data.
To run this kind of experiments, you can change the fitness provider configuration in the [evolutionary_desing/evolution.params file](evolutionary_desing/evolution.params) by uncommenting the "Fake provider" and commenting the "Actual provider".
Notably, this kind of run works without installing Gaussian, Tinker, Openbabel, and AutoCompChem.
