# Automated Design of NHC Ligands that Shorten Ru=C Bond
This repository collects input files used to automatically design N-heterocyclic carbene (NHC) ligands (L) 
meant to generate ruthenium olefin metathesis catalysts with formula (L)Ru(Cl)(Cl)=CH<sub>2</sub> bearing a particularly short Ru=CH<sub>2</sub> bond.

## Prerequisites and Preparation
### Overview
The experiments were run on an HPC Linux cluster, but most tasks, i.e., all excluding the highly computational demamding DFT computations can be reproduced on a regular laptop/workstation. The following software in needed to reproduce the experiments, and will be installed in the next section:
* <a href="https://github.com/denoptim-project/DENOPTIM">DENOPTIM 2.2.6</a> to run the evolutionary design and fragment space enumeration.
* <a href="https://dasher.wustl.edu/tinker/">Tinker v6.3.3</a> for performing molecular mechanic calculations like conformational searches.
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

3. <a href="https://dasher.wustl.edu/tinker/">Tinker v6.3.3</a> and <a href="https://gaussian.com/">Gaussian 09</a> must be obtained according to the respective license terms. Please refer to the above links for license terms and installation instructions. After installation, make the Tinker excecutable reachable, by exporting an evironmental variable pointing to Tinker's `bin` forlder. For example:
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

### Test Your Settings
To test the functionality of the resulting environment, you can run an actual fitness evaluation task using the [runTest.sh](evolutionary_desing/test/runTest.sh) script located in the `evolutionary_desing/test/` folder.
The fitness evaluation task includes: read-in a denoptim graph and produce a 3D molecular model search for a decent conformation, then evaluate atom clashes, prepare the input file for DFT modeling with Gaussian, submit the Gaussian job and post-process the its output, and finally compute the numerical value of the fitness.
The test spends about 10 minutes in the DFT optimization on 40 cpus with shared memory. This can be avoided by using the Gaussian output file provided under `evolutionary_desing/test/` folder. See [runTest.sh](evolutionary_desing/test/runTest.sh) for details.

## Running Evolutionary Design Experiments
The files to run such experiments are collected under `evolutionary_desing`: 
```
cd evolutionary_desing
```
After having performed the above preparation, you can run an evolutionary experiment. A certain cpu capacity is required by the very many DFT geometry optimization tasks. This means you want to run these tasks on an HPC cluster (but to play with the data you can do a [dry run](#dry-runs). With the exclusion of the DFT modeling, the tasks performed by denoptim are sufficiently inexpensive in terms of cpu loading that they can be run on any personal computer or even a laptop.
Therefore, you can run denoptim on a machine that allows you to seamlessly submit the Gaussian jobs to a more powerfull computer center.

To run an evolutionary design experiment in an interactive session that has been configured as explained [above](#prerequisites-and-preparation) (you may want to do this to be able to easily kill the process with `ctrl+C`):
```
java -jar ../tools/DENOPTIM/build/DenoptimGA.jar evolution.params
```
Otherwise, to run the same experiment but in a new, independent session (this allows to log out without killing the denoptim experiment) you can run the following in a terminal that has been configured as explained [above](#prerequisites-and-preparation):
```
nohup bash -c "java -jar ../tools/DENOPTIM/build/DenoptimGA.jar evolution.params &> evolution.log " & echo $! >> runningPIDs.txt ; date >> runningPIDs.txt ; echo "--" >> runningPIDs.txt ; rm nohup.out
```
Note that the PID of the master process is printed into the `runningPIDs.txt` to allow quick identification and managment of the process.


### Dry runs
You might want to play with the data without engaging in computationally demanding DFT calculations. For this purpose, a super light fitness provider can be chosen and used in evolutionary design. The fitness function will rank the candidates based solely on their molecular weight: Obviously, this does not lead to any chemically meaningful design. Still, it allows you to play with the data.
To run this kind of experiments, you can change the fitness provider configuration in the [evolutionary_desing/evolution.params file](evolutionary_desing/evolution.params) by uncommenting the "Fake provider" and commenting the "Actual provider".
Notably, this kind of run works without installing Gaussian, Tinker, Openbabel, and AutoCompChem.
