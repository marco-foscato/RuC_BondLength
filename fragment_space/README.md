# Fragment Space
This folder collets the files defining the molecular building blocks and the rules used to combine such building blocks to construct molecules.

File [FS.params](FS.params) contains the parameters for definingin the fragment space in DENOPTIM. Note the use of He atoms to replace specific H ones. These He atoms are used to desymmetrize the aromatic substituents so that calculation of InChiKeys becomes able to distinguish some N–A rotamers. The He atoms are to be replaced with H atoms prior to the calculation of the fitness, i.e., the actual element represented by the He labels is H, not He.

File [FS_noHe.params](FS_noHe.params) contains the same information of [FS.params](FS.params) but refers to the list of fragments that do not contain He labels.

The [enumeration](enumeration/README.md) folder contains the files and results from the enumeration of the fragment space.
