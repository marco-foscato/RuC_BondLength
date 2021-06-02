#!/bin/bash

###############################################################
#                                                             #
#      Fitness Calculation for Ru-Catalyst Productivity       #
#                  based on Ru=CH2 bond  length               #
#                                                             #
###############################################################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#                                                             #
# WARNING! Make sure that the following                       #
#          parameters are in agreement with those given       #
#          to DENOPTIM (<name>.params file)                   # 
#                                                             #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

# Parameters for molecular builder
if [ -z "$RUCBONDDESIGN" ]; then
    echo "ERROR! We expect environmental variable RUCBONDDESIGN to be set."
    exit -1
fi
scaffoldLib="$RUCBONDDESIGN/../fragment_space/scaffold.sdf"
fragmentLib="$RUCBONDDESIGN/../fragment_space/fragments.sdf"
cappingLib="$RUCBONDDESIGN/../fragment_space/capping_groups.sdf"
cpm="$RUCBONDDESIGN/../fragment_space/compatibility_matrix.par"

#Handling of duplicates
allowDuplicates=1 #0 for yes, 1 for no

#Data from previous runs - if a candidate UID has been already evaluated it will be found here
previousResults="$RUCBONDDESIGN/results/run_*"

#Setting for the execution of DENOPTIM tools
java="java"
pathToJarFiles="$RUCBONDDESIGN/../tools/DENOPTIM/build/"

#Openbabel
obabel="obabel"

#Tinker
if [ -z "$TINKERBIN" ]; then
    echo "ERROR! We expect environmental variable TINKERBIN to be set."
    exit -1
fi
pathToTinkerBin="$TINKERBIN"
tinkerForceField="$RUCBONDDESIGN/uff_vdw.prm"
tinkerKeyFile="$RUCBONDDESIGN/build_uff.key"
tinkerSubmitFile="$RUCBONDDESIGN/submit_pssrot"
rotSpaceDef="$RUCBONDDESIGN/rotatableBonds-1.0"

#Sed syntax differs between GNU and BDS sed
sedSyntax=GNU
sed --version >/dev/null 2>&1
if [ $? -ne 0 ]; then
    sedSyntax=Other
fi

#AutoCompChem
ACCpath="$RUCBONDDESIGN/../tools/AutoCompChem/"

#Fitness calculator
fitnessCalculatorPath="$RUCBONDDESIGN/../tools/FitnessRuCH2BndLng"

#Gaussian
submitDFTscript="submit_job_g16-C.01" #TODO back to submitDFTscript
jobDetails="$RUCBONDDESIGN/jd_G09_Denoptim_1.1"
nodes=1         # nodes
cpn=32          # processors per node
wt=20           # walltime in hours
minTime=10      # delay of first checking iteration
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
    FILE="$1"
    if [ -f "$FILE" ];
    then
        rm "$FILE"
    fi   
}

    
if [ "$#" -lt 5 ]
then
    echo " "
    echo "Usage: `basename $0` required number of arguments not supplied"	
    echo "4 parameters must be supplied (in this order):"
    echo " <inputFileName.sdf>  <outputFileName.sdf> <workingDirectory> <taskID> <UIDFile>"
    echo " "
    exit 1
fi

# Graph representation of molecule (in SDF file)
inpSDF="$1"
# Will contain the optimized 3D opt. model with FITNESS/MOL_ERROR field
outSDF="$2"
# Tasks 'working directory'
# Note tha wrkDir != locDir that is DENOPTIM's working dir (from which this script is run)
wrkDir="$3"
locDir="$(pwd)"
# Task ID (useless number in this context, but DENOPTIM will spit it out anyway)
taskId="$4"
UIDFILE="$5"

fname=$(basename "$inpSDF" .sdf)
molName=$(basename "$inpSDF" _I.sdf)

#Log of this script
log="$wrkDir/${fname}_FProvider.log"
# From here redirect stdout and stderr to log file
exec > "$log"
exec 2>&1

#Source of 3D models
LookIn3DDBParFile="$wrkDir/${fname}_LookFor3D.par"

#3D builder
DenoptimCG3Dout="$wrkDir/${fname}_3DpreDFT.sdf"
DenoptimCGParFile="$wrkDir/${fname}_DenoptimCG_"$taskId".par"

#AtomClash check
ACCatmclshParFile="$wrkDir/${fname}_AtmClash.par"
ACCatmclshLog="$wrkDir/${fname}_AtmClash.log"

dftJobName="$wrkDir/${fname}_DFT"
SDFtoDFT="$dftJobName.sdf"
XYZpostDFT="$dftJobName.xyz"
SDFpostDFT="$wrkDir/${fname}_DFTopt.sdf"

#Connectivity check
ACCconnectParFile="$wrkDir/${fname}_ConnectCheck.par"
ACCconnectLog="$wrkDir/${fname}_ConnectCheck.log"
ACCconnectOut="$wrkDir/${fname}_ConnectCheck.out"

#Fitnes calculation
fitParFile="$wrkDir/${fname}_CalcFit.par"

#
# Create UID from input mol
#
if ! grep -q "<UID>" "$inpSDF" 
then
    errmsg="#GrepUID: no UID found"
    "$obabel" -isdf "$inpSDF" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_FATAL
fi
MOLUID="$( grep -A1 "<UID>" "$inpSDF" | tail -n 1)"
uidFile="$wrkDir/${molNum}.uid"
echo "$MOLUID" > "$uidFile"
if [ "$allowDuplicates" == 1 ]
then
    echo "Generating UID for $molNum..."
    "$java" -jar "$pathToJarFiles/UpdateUID.jar" -m "$uidFile" -s "$inpSDF" -k "$UIDFILE"
    if grep -q "Already exists." "$inpSDF"
    then
        echo "UpdateUID: duplicate candidate."
        cp "$inpSDF" "$outSDF"
        echo "All done!"
        exit $E_OPTERROR
    elif grep -q "MOL_ERROR" "$inpSDF"
    then
        echo "UpdateUID: ERROR from checking UID"
        cp "$inpSDF" "$outSDF"
        exit $E_FATAL 
    fi
fi

#
# Use results from previous evaluation of this candidate (if available)
#
echo "Searching for UID=$MOLUID in $previousResults"
fitFilesFromPreviousRun="$( find $previousResults -name "*_FIT.sdf" -exec grep -l "$MOLUID" {} \; | xargs grep -l "FITNESS" )"
 if [[ -n $fitFilesFromPreviousRun ]]
then
    echo "Trying to fetch fitness from ${fitFilesFromPreviousRun}."
    firstFitFile="$( ls $fitFilesFromPreviousRun -rlt | awk '{print $9}' | head -n 1 )"
    preOutSDF="${outSDF}_pre"
    cat "$firstFitFile" > "$preOutSDF"
    firstXYZFile="$(basename $firstFitFile _FIT.sdf)_DFT.xyz"
    if [ -f "$firstXYZFile" ]
    then
        cp "$firstXYZFile" "$XYZpostDFT"
    fi 
    prevMolName="$( head -n 1 $firstFitFile )"
    if [ "$sedSyntax" == "GNU" ]
    then
        sed -i "s/$prevMolName/$molName/g" "$preOutSDF"
    else
        sed -i '' "s/$prevMolName/$molName/g" "$preOutSDF"
    fi
    echo "Fitness collected from previous run: $firstFitFile"
    mv "$preOutSDF" "$outSDF"
    echo "All done!"
    exit 0
fi

## calculation of 3D conformer using TINKER
echo "Starting DenoptimCG"
# prepare param file 
echo "CG-inpSDF=$inpSDF" > "$DenoptimCGParFile"
echo "CG-outSDF=$DenoptimCG3Dout" >> "$DenoptimCGParFile"
echo "CG-workDir=$wrkDir" >> "$DenoptimCGParFile"
echo "FS-ScaffoldLibFile=$scaffoldLib" >> "$DenoptimCGParFile"
echo "FS-FragmentLibFile=$fragmentLib" >> "$DenoptimCGParFile"
echo "FS-CappingFragmentLibFile=$cappingLib" >> "$DenoptimCGParFile"
echo "FS-CompMatrixFile=$cpm" >> "$DenoptimCGParFile"
# definition of rotational space
echo "FS-RotBondsDefFile=$rotSpaceDef" >> "$DenoptimCGParFile"
# location of the TINKER tools
echo "CG-ToolPSSROT=$pathToTinkerBin/pssrot" >> "$DenoptimCGParFile"
echo "CG-ToolXYZINT=$pathToTinkerBin/xyzint" >> "$DenoptimCGParFile"
echo "CG-ToolINTXYZ=$pathToTinkerBin/intxyz" >> "$DenoptimCGParFile"
# param file used by Tinker
echo "CG-ForceFieldFile=$tinkerForceField" >> "$DenoptimCGParFile"
# key file to be used by tinker with PSSROT
# this file is copied and edited for every molecule
echo "CG-KeyFile=$tinkerKeyFile" >> "$DenoptimCGParFile"
# parameters used by PSSROT
# this file is copied and edited for every molecule
echo "CG-PSSROTParams=$tinkerSubmitFile" >> "$DenoptimCGParFile"

"$java" -jar "$pathToJarFiles/DenoptimCG.jar" "$DenoptimCGParFile"

if [ ! -f "$DenoptimCG3Dout" ]; then
    echo "$DenoptimCG3Dout not found."
    errmsg="#DenoptimCG: "$DenoptimCG3Dout" not found."
    "$obabel" -isdf "$inpSDF" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"    
    exit $E_OPTERROR
fi


## Change He to H: needed only for special use of fragments bearing dummy He atoms
echo "Changing He to H"
if [ "$sedSyntax" == "GNU" ]
then
    sed -i 's/ He / H  /g' "$DenoptimCG3Dout"
else
    sed -i '' 's/ He / H  /g' "$DenoptimCG3Dout"
fi

## Check for AtomClashes
echo "Starting AtomClash detector"
# prepare param file
echo "VERBOSITY: 1" > "$ACCatmclshParFile"
echo "TASK: AnalyzeVDWClashes" >> "$ACCatmclshParFile"
echo "INFILE: $DenoptimCG3Dout" >> "$ACCatmclshParFile"
echo "ALLOWANCE13: 1.0" >> "$ACCatmclshParFile"
echo "CUTOFF: 0.70" >> "$ACCatmclshParFile"
echo "ALLOWANCE: 0.40" >> "$ACCatmclshParFile"

"$java" -jar "$ACCpath/AutoCompChem.jar" "$ACCatmclshParFile" > "$ACCatmclshLog"

if [ ! -f "$ACCatmclshLog" ]; then
    errmsg="#AtomClash Check: "$ACCatmclshLog" not found."
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if ! grep -q "Termination status: 0" "$ACCatmclshLog" ; then
    errmsg="#AtomClash Check: non-zero exit status from AutoCompChem"
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if grep -q 'Found 0 mols with one or more atom clashes' "$ACCatmclshLog" ; then
    echo "No atom clashes"
else
    echo "Found Atom Clashes"
    errmsg="#AtomClash Check: Found Atom Clashes"
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi


## Prepare input file for DFT
cp "$DenoptimCG3Dout" "$SDFtoDFT"
ACCmakeInpLog="$wrkDir/${fname}_mkInp.log"
ACCmakeInpParFile="$wrkDir/${fname}_mkInp.par"
dftInpFile="$dftJobName.inp"
dftOutFile="$dftJobName.out"
dftLogFile="$dftJobName.log"
echo "Preparing input file for DFT"
# prepare param file
echo "VERBOSITY: 1" > "$ACCmakeInpParFile"
echo "TASK: PrepareInputGaussian" >> "$ACCmakeInpParFile"
echo "INFILE: $SDFtoDFT" >> "$ACCmakeInpParFile"
echo "JOBDETAILSFILE: $jobDetails" >> "$ACCmakeInpParFile"
echo "OUTNAME: $dftInpFile" >> "$ACCmakeInpParFile"

"$java" -jar "$ACCpath/AutoCompChem.jar" "$ACCmakeInpParFile" > "$ACCmakeInpLog"

if [ ! -f $dftInpFile ]; then
    errmsg="#MakeDFTInp: "$ACCmakeInpLog" not found."
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if [ ! -f "$ACCmakeInpLog" ]; then
    errmsg="#MakeDFTInp: "$ACCmakeInpLog" not found."
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if ! grep -q "Termination status: 0" "$ACCmakeInpLog" ; then
    errmsg="#MakeDFTInp: non-zero exit status from AutoCompChem"
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi

## Submit DFT geometry optimization
echo "Submitting DFT"
date
cd "$wrkDir"
pwd
echo "Running $submitDFTscript $(basename $dftInpFile) -n $nodes -t $wt"
jobID=$("$submitDFTscript" "$(basename $dftInpFile)" -n $nodes -t $wt)
#NB: this only works for SLURM. PBS queue may returns a different message
if [[ "$jobID" == *"Batch job submission failed"* ]] ; then
    errmsg="#SubmitDFT: Batch job submission failed"
    echo $errmsg
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_FATAL
fi
cd "$locDir"
if [ $? != 0 ]; then
    errmsg="#SubmitDFT: non-zero exit status from submission script"
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi

## Wait for completion of DFT job (= creation of the "task complete label(TCL)")
jobID=$(echo "$jobID" | grep "Submitted " | awk '{print $NF}')
taskdone=1
msg="";
# Wait the minimum time before looping
echo "Start waiting for completion of DFT: jobID $jobID"
echo "Waiting for string 'GAUSSIAN JOB ENDED' in '$dftLogFile'"
sleep $minTime$minTimeUnit
echo "First check for completion of DFT"
date
# WARNING: this if.statment depends on the job submission script and its log. You might need to change it!
if [ -f "$dftLogFile" ] && grep -q "GAUSSIAN JOB ENDED" "$dftLogFile" ; then
    echo "Job done: stop waiting"
    taskdone=0
else
    echo "Starting job completion checking loop"
    # Start checking for results
    for i in $(seq 1 $maxwait)
    do
       sleep $step$stepunit
       date
       if [ -f "$dftLogFile" ] && grep -q "GAUSSIAN JOB ENDED" "$dftLogFile"
       then
          echo "Job completed: stop waiting"
          taskdone=0
          break
       else
          if [ $i == $maxwait ]
          then
             errmsg="#WaitingDFT: time limit reached (task abbandoned)"
             "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
             exit $E_OPTERROR
          fi
       fi
    done
fi

## XYZ to SDF
echo "Conversion of XYZ"
echo "DFT output file: $dftOutFile"
echo "Making PostDFT XYZ: $XYZpostDFT"
if [ ! -f "$dftOutFile" ]; then
    errmsg="#DFTResults: file '$dftOutFile' not found"
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR 
fi
"$obabel" -ig09 "$dftOutFile" -oxyz -O "$XYZpostDFT"
if [ $? != 0 ]; then
    errmsg="#G09toXYZ failure: non-zero exit status"
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
echo "Making PostDFT SDF: $SDFpostDFT"
#NB: here we rely on the connectivity perception of OpenBabel!
"$obabel" -ixyz "$XYZpostDFT" -osdf -O "$SDFpostDFT"
if [ $? != 0 ]; then
    errmsg="#XYZtoSDF failure: non-zero exit status"
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if [ ! -f "$SDFpostDFT" ]; then
    errmsg="#XYZtoSDF failure: "$SDFpostDFT" not found"
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi


## Compare connectivity of optimized geometry against initial
echo "Comparing connectivity"
# prepare parameters for connectivity check
echo "VERBOSITY: 2 " > "$ACCconnectParFile"
echo "TASK: CompareTwoConnectivities" >> "$ACCconnectParFile"
echo "INFILE: $SDFpostDFT" >> "$ACCconnectParFile"
echo "REFERENCE: $DenoptimCG3Dout" >> "$ACCconnectParFile"
echo "OUTFILE: $ACCconnectOut" >> "$ACCconnectParFile"
# NOTE: the output is removed right after the execution of the task. 
#       "$ACCconnectOut" not needed but OUTPUT option required.

"$java" -jar "$ACCpath/AutoCompChem.jar" "$ACCconnectParFile" > "$ACCconnectLog"

if [ ! -f "$ACCconnectLog" ]; then
    errmsg="#Connectivity Check: "$ACCconnectLog" not found."
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if ! grep -q "Termination status: 0" "$ACCconnectLog" ; then
    errmsg="#Connectivity Check: non-zero exit status from AutoCompChem"
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if grep -q "Inconsistent adjacency" "$ACCconnectLog" ; then
    errmsg=$(grep "Inconsistent adjacency" "$ACCconnectLog")
    errmsg="#Connectivity Check: $errmsg" 
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi
if grep -q "Consistent connectivity" "$ACCconnectLog" ; then
    echo "Consistent connectivity"
else
    errmsg="#Connectivity Check: Consistent connectivity flag not found"
    "$obabel" -isdf "$DenoptimCG3Dout" -osdf -O "$outSDF" --property "MOL_ERROR" "$errmsg"
    exit $E_OPTERROR
fi

#removing temporary file
rm -f "$ACCconnectOut"


## Fitness calculation, this also prepares the final SDF file 
echo "Calculating descriptors and fitness"
# Prepare the parameter file
echo "INPSDF=$DenoptimCG3Dout" > "$fitParFile"
echo "OUTSDF=$outSDF" >> "$fitParFile"
echo "OPTSDF=$SDFpostDFT" >> "$fitParFile"
echo "HPXYZ=$XYZpostDFT" >> "$fitParFile"
echo "MAXBNDDIST=2.50" >> "$fitParFile"
echo "MINANGLE=90.0" >> "$fitParFile"
echo "MAXTORSION=20" >> "$fitParFile"
echo "MINNBDIST=2.7" >> "$fitParFile"
echo "WORKDIR=$wrkDir" >> "$fitParFile"

"$java" -jar "$fitnessCalculatorPath/FitnessRuCH2BndLng.jar" "$fitParFile"

if [ $? != 0 ];then
    echo "FitnessRuCH2BndLng.jar failed execution."
    errmsg="#Fitness Evaluation: non-zero exit status from FitnessRuCH2BndLng"
    exit $E_OPTERROR
fi


# cleanup 
echo "Cleanup" >> "$log"
cleanup "$uidFile"
cleanup "$LookIn3DDBParFile"
cleanup "$DenoptimCGParFile"
#cleanup "$DenoptimCG3Dout"
cleanup "$ACCatmclshParFile"
cleanup "$ACCatmclshLog"
cleanup "$ACCmakeInpParFile"
cleanup "$ACCmakeInpLog"
cleanup "$SDFtoDFT"
cleanup "$dftJobName.job"
cleanup "$dftJobName.jd"
cleanup "$dftJobName.inp"
#cleanup "$XYZpostDFT"
cleanup "$SDFpostDFT"
cleanup "$ACCconnectParFile"
cleanup "$ACCconnectLog"
cleanup "$fitParFile"
#cleanup $log

# Task done
echo "All done!"
exit 0
