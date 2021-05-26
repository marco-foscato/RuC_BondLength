# Enumeration of Fragment Space
Here we want to generate all the molecular constitutions that are encoded in the fragment space. We ignore the issue of hindered torsion along the N–A substituent bond. Therefore, we use fragments that do not have the dummy He atom that de-symmetrizes the Ph substituent. In other words, we use the fragment space defined in `../FS_noHe.params`

## Procedure
To facilitate parallelization of the enumeration tasks, we first have generated the list of all possible root graphs containg the scaffold (i.e., the Ru center) and first layer of fragments (i.e., the NHC core with Ph skeleton on both N atoms). This results in 18 rood graphs. For most of these, the complete set of generated graphs can be obtained by
```
java -jar $DENOPTIMHOME/build/FragSpaceExplorer.jar FSE_$i.params > FSE_$i.log 2>&1
```
where `$DENOPTIMHOME` is the installation folder of denoptim and `$i` is an integer identifying the roog graph (`i = 1–16, 18`). Since one root generates very may graphs, i.e., root_17 generates 189000 graphs, the enumeration of the graphs was further split into independent runs by generating all the 56 combiantions of subtituents in the backbone and using these 56 graphs as roots for as many explorations of the fragment space.

## Results
The `tar.gz` archives contain the unrefineg chemical structure of all 303840 combinations of fragments. These files are named with a unique time stap and reflect the following order of root graphs: 1–16,18,17_1–56.

This is the summary of the total number of fragment combinations generated for each root:
| Root # | Total number of combinations |
| --- | --- |
|1|	27000|
|2|	3375|
|3|	23625|
|4|	27000|
|5|	3375|
|6|	23625|
|7|	3375|
|8|	3375|
|9|	1|
|10|	7|
|11|	8|
|12|	1|
|13|	56|
|14|	1|
|15|	8|
|16|	7|
|17\*|	189000\*|
|18|	1|
| **Sum** |	**303840**|

\*: each of the 56 graphs used to split this task contributes with 3375 graphs.

These numbers, however, do contain duplicates. For instance, many roots, when all attachment points are saturated with the capping group (i.e., an H atom), generate the same molecule. The corresponding InChIKeys can be collected by 
```
find . -type d -name "FSE2*" | while read d ; do grep -A1 "<InChi>" $d/*_inp.sdf | grep -v "<InChi>" | grep -v "^\-\-" >> InChiKeys_with_duplicates ; done
```
To purge the duplicates we collect and count the unique InChI keys, which are generated for each candidate.
```
awk -F".sdf-" '{print $2}' InChiKeys_with_duplicates | sort -u > InChiKeys_unique
```
The result is a list of **154882 unique InChIKeys** generated without distinguishing conformers along the N-Ph bond, i.e., with fragments having no desymmetrizing He atom.
