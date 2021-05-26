#!/bin/bash

###############################################################
#                                                             #
#      Fitness Calculation for Ru-Catalyst Productivity       #
#                  based on Ru=CH2 bond  length               #
#                                                             #
###############################################################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#                                                             #
# WARNING! Alternatively, make sure that the following        #
#          parameters are in agreement with those given       #
#          to DENOPTIM (<name>.params file)                   # 
#                                                             #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

# Parameters for molecular builder
scaffoldLib="/global/work/marcof/BL-2.2/scaff.sdf"
fragmentLib="/global/work/marcof/BL-2.2/frag.sdf"
cappingLib="/global/work/marcof/BL-2.2/cap.sdf"
cpm="/global/work/marcof/BL-2.2/CPMap.par"

#Setting for the execution of DENOPTIM tools
java="/usr/java/latest/bin/java"
pathToJarFiles="/home/marcof/DENOPTIM/newdenoptim_RuBndLng/build"

#Openbabel
obabel="/home/marcof/OpenBabel/OB-2.3.2_built_Nov14/bin/obabel"

#TINKER - PSSROT
pathToTinkerBin="/home/marcof/Tinker/tinker-6.3.3/bin"
tinkerForceField="/global/work/marcof/BL-2.2/uff_vdw.prm"
tinkerKeyFile="/global/work/marcof/BL-2.2/build_uff.key"
tinkerSubmitFile="/global/work/marcof/BL-2.2/submit_pssrot"
rotSpaceDef="/global/work/marcof/BL-2.2/rotatableBonds-1.0"

#AutoCompChem
ACCpath="/home/marcof/foscato/AutoCompChem/jar"

#Source of 3D models
prev3DSDF="/global/work/marcof/BL-2.2/available3DModels.sdf"
prev3DXYZ="/global/work/marcof/BL-2.2/available3DModelsXYZ"

#DFT with Gaussian09 via AutoCompChem
submitDFTscript="/home/marcof/jobscript/g09_AutoCompChem-PBS-1.2.sh"
jobDetails="/global/work/marcof/BL-2.2/jd_G09_Denoptim_1.1"
nodes=2         # nodes
cpn=16          # processors per node
wt=20           # walltime in hours
mpn=32000       # memory per node
minTime=30      # delay of first checking iteration
minTimeUnit="m" # time unit: s (seconds), m (minutes), h (hours)
step=5          # delay of each checking iteration
stepunit="m"    # time unit: s (seconds), m (minutes), h (hours)
maxwait=350     # 350*5m=29.2h # maximun number of checking iterations

#Exit code for uncomplete evaluation of fitness
# -> set to 0 to return *FIT.sdf file with MOL_ERROR field
# -> set to anything else to make DENOPTIM stop in case of uncomplete evaluation of fitness
E_OPTERROR=0
# Exit code for fatal errors
E_FATAL=1


#############################################
#############################################
##                                         ##
## No need to change paths below this line ##
##                                         ## 
#############################################
#############################################

function cleanup() {
    # t stores $1 argument passed to cleanup()
    FILE=$1
    if [ -f $FILE ];
    then
        rm $FILE
    fi   
}

    
if [ "$#" -lt 4 ]
then
    echo " "
    echo "Usage: `basename $0` required number of arguments not supplied"	
    echo "4 parameters must be supplied (in this order):"
    echo " <inputFileName.sdf>  <outputFileName.sdf> <workingDirectory> <taskID>"
    echo " "
    exit 1
fi

# Graph representation of molecule (in SDF file)
inpSDF=$1
# Will contain the optimized 3D opt. model with FITNESS/MOL_ERROR field
outSDF=$2
# Tasks 'working directory 
# Note tha wrkDir != locDir thae is DENOPTIM's working dir (from which this script is run)
wrkDir=$3
locDir=$(pwd)
# Task ID (useless number in this context, but DENOPTIM will spit it out anyway)
taskId=$4

fname=`basename $inpSDF .sdf`

#Log of this script
log=$wrkDir/$fname"_FProvider.log"
# From here redirect stdout and stderr to log file
exec > $log
exec 2>&1

#Source of 3D models
LookIn3DDBParFile=$wrkDir/$fname"_LookFor3D.par"

#3D builder
DenoptimCG3Dout=$wrkDir/$fname"_3DpreDFT.sdf"
DenoptimCGParFile=$wrkDir/$fname"_DenoptimCG_"$taskId".par"

#AtomClash check
ACCatmclshParFile=$wrkDir/$fname"_AtmClash.par"
ACCatmclshLog=$wrkDir/$fname"_AtmClash.log"

#DFT
queue="PBS"
dftJobName=$wrkDir/$fname"_DFT"
SDFtoDFT=$dftJobName".sdf"
XYZpostDFT=$dftJobName".xyz"
SDFpostDFT=$wrkDir/$fname"_DFTopt.sdf"

#Connectivity check
ACCconnectParFile=$wrkDir/$fname"_ConnectCheck.par"
ACCconnectLog=$wrkDir/$fname"_ConnectCheck.log"
ACCconnectOut=$wrkDir/$fname"_ConnectCheck.out"

#Fitnes calculation
fitParFile=$wrkDir/$fname"_CalcFit.par"


## Search for the molecule in database of previously modelled molecules
echo "Looking for the 3D model in list of previously modelled molecules"
echo "INPSDF=$inpSDF" > $LookIn3DDBParFile
echo "IDFIELD=InChi" >> $LookIn3DDBParFile
echo "DBSDF=$prev3DSDF" >> $LookIn3DDBParFile
echo "DBXYZPATH=$prev3DXYZ" >> $LookIn3DDBParFile
echo "OUTSDF=$SDFpostDFT" >> $LookIn3DDBParFile
echo "OUTXYZ=$XYZpostDFT" >> $LookIn3DDBParFile

rm -f $XYZpostDFT
rm -f $SDFpostDFT
$java -jar $pathToJarFiles/LookInExisting3DMols.jar $LookIn3DDBParFile

if [ $? != 0 ];then
    echo "LookInExisting3DMols failed execution."
    errmsg="#SearchOld3D: non-zero exit status from LookInExisting3DMols"
    exit $E_OPTERROR
fi

if [ -f $SDFpostDFT ]; then
    # No need for modeling 3D structure
    echo "Using previously modelled 3D structure"

    # create copy needed to calculate fitness
    $obabel -isdf $inpSDF -osdf -O $DenoptimCG3Dout --property "Coord3D" "XYZ will be taken from previously modelled 3D"

    ## Change He to H: needed only for special use of fragments bearing dummy He atoms
    echo "Changing He to H"
    sed -i 's/ He / H  /g' $DenoptimCG3Dout

else
    echo "No previous 3D model found for this molecule"
    ## Calculation of 3D conformer using TINKER
    echo "Starting DenoptimCG"
    # prepare param file 
    echo "inpSDF=$inpSDF" > $DenoptimCGParFile
    echo "outSDF=$DenoptimCG3Dout" >> $DenoptimCGParFile
    echo "scaffoldLibFile=$scaffoldLib" >> $DenoptimCGParFile
    echo "fragmentLibFile=$fragmentLib" >> $DenoptimCGParFile
    echo "cappingFragmentLibFile=$cappingLib" >> $DenoptimCGParFile
    echo "compMatrix=$cpm" >> $DenoptimCGParFile
    echo "toolOpenBabel=$obabel" >> $DenoptimCGParFile
    echo "wrkDir=$wrkDir" >> $DenoptimCGParFile
    # definition of rotational space
    echo "ROTBONDSDEF=$rotSpaceDef" >> $DenoptimCGParFile
    # location of the TINKER tools
    echo "PSSROT=$pathToTinkerBin/pssrot" >> $DenoptimCGParFile
    echo "XYZINT=$pathToTinkerBin/xyzint" >> $DenoptimCGParFile
    echo "INTXYZ=$pathToTinkerBin/intxyz" >> $DenoptimCGParFile
    # param file used by Tinker
    echo "PARAM=$tinkerForceField" >> $DenoptimCGParFile
    # key file to be used by tinker with PSSROT
    # this file is copied and edited for every molecule
    echo "KEYFILE=$tinkerKeyFile" >> $DenoptimCGParFile
    # parameters used by PSSROT
    # this file is copied and edited for every molecule
    echo "PSSROTPARAMS=$tinkerSubmitFile" >> $DenoptimCGParFile
    # Atom ordering scheme (1/2)
    echo "atomOrderingScheme=1" >> $DenoptimCGParFile
    
    $java -jar $pathToJarFiles/DenoptimCG.jar $DenoptimCGParFile
    
    if [ ! -f $DenoptimCG3Dout ]; then
        echo "$DenoptimCG3Dout not found."
        errmsg="#DenoptimCG: $DenoptimCG3Dout not found."
        $obabel -isdf $inpSDF -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"    
        exit $E_OPTERROR
    fi
    
    
    ## Change He to H: needed only for special use of fragments bearing dummy He atoms
    echo "Changing He to H"
    sed -i 's/ He / H  /g' $DenoptimCG3Dout
    
    
    ## Check for AtomClashes
    echo "Starting AtomClash detector"
    # prepare param file
    echo "VERBOSITY: 1" > $ACCatmclshParFile
    echo "TASK: AnalyzeVDWClashes" >> $ACCatmclshParFile
    echo "INFILE: $DenoptimCG3Dout" >> $ACCatmclshParFile
    echo "ALLOWANCE13: 1.0" >> $ACCatmclshParFile
    echo "CUTOFF: 0.70" >> $ACCatmclshParFile
    echo "ALLOWANCE: 0.40" >> $ACCatmclshParFile
    
    $java -jar $ACCpath/AutoCompChem.jar $ACCatmclshParFile > $ACCatmclshLog
    
    if [ ! -f $ACCatmclshLog ]; then
        errmsg="#AtomClash Check: $ACCatmclshLog not found."
        $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
        exit $E_OPTERROR
    fi
    if ! grep -q "Termination status: 0" $ACCatmclshLog ; then
        errmsg="#AtomClash Check: non-zero exit status from AutoCompChem"
        $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
        exit $E_OPTERROR
    fi
    if grep -q 'Found 0 mols with one or more atom clashes' $ACCatmclshLog ; then
        echo "No atom clashes"
    else
        echo "Found Atom Clashes"
        errmsg="#AtomClash Check: Found Atom Clashes"
        $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
        exit $E_OPTERROR
    fi
    
    
    ## Submit DFT geometry optimization
    echo "Submitting DFT"
    date
    cp $DenoptimCG3Dout $SDFtoDFT
    cd $wrkDir
    jobID=$(/bin/bash $submitDFTscript $SDFtoDFT $jobDetails $nodes $cpn $wt $mpn)
    #NB: this only works for SLURM. PBS queue may returns a different message
    if [[ "$jobID" == *"Batch job submission failed"* ]] ; then
        errmsg="#SubmitDFT: Batch job submission failed"
        echo $errmsg
        $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
        exit $E_FATAL
    fi
    cd $locDir
    if [ $? != 0 ]; then
        errmsg="#SubmitDFT: non-zero exit status from submission script"
        $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
        exit $E_OPTERROR
    fi
    
    ## Wait for completion of DFT job (= creation of the "task complete label(TCL)")
    if [ $queue == "PBS" ]
    then
        jobID=$(echo $jobID | awk '{print $NF}')
    elif [ $queue == "SLURM" ]
    then
        jobID=$(echo $jobID | grep "Submitted " | awk '{print $NF}')
    fi
    tcl="$wrkDir/tcl_$jobID"
    taskdone=1
    msg="";
    # Wait the minimum time before looping
    echo "Start waiting for completion of DFT: jobID $jobID"
    echo "Name of Task Complete Label: $tcl"
    sleep $minTime$minTimeUnit
    echo "First check for completion of DFT"
    date
    if [ -f $tcl ]; then
        echo "Job done: stop waiting"
        taskdone=0
    else
        echo "Starting job completion checking loop"
        # Start checking for results
        for i in $(seq 1 $maxwait)
        do
           sleep $step$stepunit
           date
           if [ -f $tcl ]
           then
              echo "Job completed: stop waiting"
              taskdone=0
              break
           else
              if [ $i == $maxwait ]
              then
                 errmsg="#WaitingDFT: time limit reached (task abbandoned)"
                 $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
                 exit $E_OPTERROR
              fi
           fi
        done
    fi   
    rm -f $tcl 
    if [ ! -f $XYZpostDFT ]; then
        errmsg="#WaitingDFT: $XYZpostDFT not found"
        $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
        exit $E_OPTERROR
    fi
    
    
    ## XYZ to SDF
    echo "Conversion of XYZ"
    echo "PostDFT XYZ: $XYZpostDFT"
    echo "PostDFT SDF: $SDFpostDFT"
    $obabel -ixyz $XYZpostDFT -osdf -O $SDFpostDFT
    if [ $? != 0 ]; then
        errmsg="#XYZtoSDF failure: non-zero exit status"
        $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
        exit $E_OPTERROR
    fi
    if [ ! -f $SDFpostDFT ]; then
        errmsg="#XYZtoSDF failure: $SDFpostDFT not found"
        $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
        exit $E_OPTERROR
    fi
fi


## Compare connectivity of optimized geometry against initial
echo "Comparing connectivity"
# prepare parameters for connectivity check
echo "VERBOSITY: 2 " > $ACCconnectParFile
echo "TASK: CompareTwoConnectivities" >> $ACCconnectParFile
echo "INFILE: $SDFpostDFT" >> $ACCconnectParFile
echo "REFERENCE: $DenoptimCG3Dout" >> $ACCconnectParFile
echo "OUTFILE: $ACCconnectOut" >> $ACCconnectParFile
# NOTE: the output is removed right after the execution of the task. 
#       $ACCconnectOut not needed but OUTPUT option required.

$java -jar $ACCpath/AutoCompChem.jar $ACCconnectParFile > $ACCconnectLog

if [ ! -f $ACCconnectLog ]; then
    errmsg="#Connectivity Check: $ACCconnectLog not found."
    $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if ! grep -q "Termination status: 0" $ACCconnectLog ; then
    errmsg="#Connectivity Check: non-zero exit status from AutoCompChem"
    $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if grep -q "Inconsistent adjacency" $ACCconnectLog ; then
    errmsg=$(grep "Inconsistent adjacency" $ACCconnectLog)
    errmsg="#Connectivity Check: $errmsg" 
    $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if grep -q "Consistent connectivity" $ACCconnectLog ; then
    echo "Consistent connectivity"
else
    errmsg="#Connectivity Check: Consistent connectivity flag not found"
    $obabel -isdf $DenoptimCG3Dout -osdf -O $outSDF --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi

#removing temporary file
rm -f $ACCconnectOut


## Fitness calculation, this also prepares the final SDF file 
echo "Calculating descriptors and fitness"
# Prepare the parameter file
echo "INPSDF=$DenoptimCG3Dout" > $fitParFile
echo "OUTSDF=$outSDF" >> $fitParFile
echo "OPTSDF=$SDFpostDFT" >> $fitParFile
echo "HPXYZ=$XYZpostDFT" >> $fitParFile
echo "MAXBNDDIST=2.50" >> $fitParFile
echo "MINANGLE=90.0" >> $fitParFile
echo "MAXTORSION=20" >> $fitParFile
echo "MINNBDIST=2.7" >> $fitParFile
echo "WORKDIR=$wrkDir" >> $fitParFile

$java -jar $pathToJarFiles/FitnessRuCH2BndLng.jar $fitParFile

if [ $? != 0 ];then
    echo "FitnessRuCH2BndLng.jar failed execution."
    errmsg="#Fitness Evaluation: non-zero exit status from FitnessRuCH2BndLng"
    exit $E_OPTERROR
fi


# cleanup 
echo "Cleanup" >> $log
cleanup $LookIn3DDBParFile
cleanup $DenoptimCGParFile
#cleanup $DenoptimCG3Dout
cleanup $ACCatmclshParFile
cleanup $ACCatmclshLog
cleanup $SDFtoDFT
#cleanup $XYZpostDFT
cleanup $SDFpostDFT
cleanup $ACCconnectParFile
cleanup $ACCconnectLog
cleanup $fitParFile
#cleanup $log

# Task done
exit 0
