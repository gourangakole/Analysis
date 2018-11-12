# Analysis
   
### Download and compile the package  ###  
* source /cvmfs/cms.cern.ch/cmsset_default.sh
* cmsrel CMSSW_9_4_6
* cd CMSSW_9_4_6/src
* cmsenv
* git clone https://github.com/ravindkv/Analysis.git 
* cd Analysis/src
* make clean 
* make
* cd .. 

### Compile and run, in one go ### 
* root -l 'runMe.C("hplusAnalyzer")'

### Submit condor batch jobs  ###

* cd condor
* ./hplusRunCond_TIFR.sh ntupleT2Paths.txt
