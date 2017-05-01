#!/bin/bash

#REFERENCE
#https://github.com/florez/CONDOR

#------------------------------------------------
#pass the arguments to Analysis/runMe.sh script 
#these arguments will go to the hplusAnalyzer.C
#------------------------------------------------

inNtupleFile=$1
outAnalFile=$2
outAnalDir=$3
date

#------------------------------------------------
#this script runs on some remote condor machine.
#link lxplus to this remote machine, using scram.
#copy the compiled lxplus package to this machine
#------------------------------------------------

cd /afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
echo "CONDOR DIR: $_CONDOR_SCRATCH_DIR"
cd ${_CONDOR_SCRATCH_DIR}


#------------------------------------------------
#copy the lxplus package to the remote machine
#and run the codes at remote machine
#------------------------------------------------

cp -r /afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src/Analysis .
#cp -r /afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src/Analysis .
cd Analysis
./runMe.sh $inNtupleFile $outAnalFile $outAnalDir
echo "OUTPUT: "
ls ${_CONDOR_SCRATCH_DIR}/Analysis/13TeV/$outAnalDir


#------------------------------------------------
#copy the output from remote machine to the lxplus
#or to any other place e.g. Tier-2
#Remove the package, after copying the output
#------------------------------------------------

xrdcp -f -R ${_CONDOR_SCRATCH_DIR}/Analysis/13TeV/$outAnalDir root://se01.indiacms.res.in:1094//cms/store/user/rverma/condor_29April17/
cd ${_CONDOR_SCRATCH_DIR}
rm -rf Analysis
echo "DONE"
date

