# Automated Design of NHC ligands for short Ru=C bond
This repository collects input files used to automatically design N-heterocyclic carbene (NHC) ligands (L) 
meant to generate ruthenium olefin metathesis catalysts with formula (L)Ru(Cl)(Cl)=CH<sub>2</sub> bearing a particularly short Ru=CH<sub>2</sub> bond.

## Prerequisites
The following software in needed to reproduce the experiments on a Linux cluster:
* <a href="https://github.com/denoptim-project/DENOPTIM">DENOPTIM 2.2.6</a> to run the evolutionary design and fragment space enumeration.
* <a href="https://dasher.wustl.edu/tinker/">Tinker v6.3</a> for performing molecular mechanic calculations like conformational searches. **NB:** Tinker must be compiled with `maxval=12` in Tinker's source file `sizes.i`.
* <a href="https://gaussian.com/">Gaussian 09</a> for DFT calculations. .
* <a href="http://openbabel.org/wiki/Main_Page">OpenBabel</a> for chemical file format conversion.
* <a href="">AutoCompChem</a> for performaning molecular modeling tasks automatically.
* JAVA JDK-8 for running DENOPTIM and AutoCompChem

The file [environment.yml](environment.yml) can be used to recreate a minor part of the environment (OpenBabel and Java), while DENOPTIM, Tinker, Gaussian, and AutoCompChem must be obtained and installed manually. Refer to the above links for instructions and license.
