#!/bin/bash

# Building FitnessRuCH2BndLng
if [ ! -f ../DENOPTIM/build/DENOPTIM-GUI.jar ]; then
	echo "Cannot locate DENOPTIM-GUI.jar file. First build DENOPTIM, then run this script."
    exit -1
fi
cp ../DENOPTIM/build/DENOPTIM-GUI.jar lib

find src/ -name *.java > javafiles.txt
javac -cp lib/cdk-1.4.19.jar:lib/vecmath.jar:lib/DENOPTIM-GUI.jar @javafiles.txt -encoding utf-8 -d .


if [ "$?" != "0" ]; then
    rm javafiles.txt
	echo "Failed to create FitnessRuCH2BndLng.jar."
    exit -1
fi

rm javafiles.txt


echo "Manifest-Version: 1.0" > manifest.mf
echo "Main-Class: fitnessruch2bndlng.FitnessRuCH2BndLng" >> manifest.mf
echo "Class-Path: lib/cdk-1.4.19.jar lib/vecmath.jar lib/DENOPTIM-GUI.jar" >> manifest.mf
echo >> manifest.mf

jar cvfm FitnessRuCH2BndLng.jar manifest.mf fitnessruch2bndlng 


if [ "$?" = "0" ]; then
     rm -rf manifest.mf fitnessruch2bndlng
else
	echo "Failed to create FitnessRuCH2BndLng.jar."
    exit -1
fi

echo "--------------------- Done building FitnessRuCH2BndLng.jar ---------------------"
