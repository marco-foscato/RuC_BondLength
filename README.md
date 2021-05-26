# Automated Design of NHC ligands for short Ru=C bond
This repository collects input files used to automatically design N-heterocyclic carbene (NHC) ligands (L) 
meant to generate ruthenium olefin metathesis catalysts with formula (L)Ru(Cl)(Cl)=CH<sub>2</sub>.

## Prerequisites
To reproduce the experiments the following software tools have to be installed:
* <a href="https://github.com/denoptim-project/DENOPTIM">DENOPTIM 2.2.6</a> to run the evolutionary design and fragment space enumeration.
* <a href="https://dasher.wustl.edu/tinker/">Tinker v6.3</a> for performing conformational searches. **NB:** Tinker must be compiled with `maxval=12` in Tinker's source file `sizes.i`.
* <a href="http://openbabel.org/wiki/Main_Page">OpenBabel</a> for chemical file format conversion.
* <a href="">AutoCompChem</a> for performaning molecular modeling tasks automatically.
* JAVA JDK-8 for running DENOPTIM and AutoCompChem

The file [environment.yml](environment.yml) can be used to recreate part of the environment (OpenBabel and Java). However, DENOPTIM, Tinker and AutoCompChem must be obtained and installed manually.
