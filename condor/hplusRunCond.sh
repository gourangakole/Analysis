#!/bin/bash
#REFERENCE
#https://github.com/florez/CONDOR

#//////////////////////////////////////
#                                     #
# hplusCond.sub runs hplusAnalyzer.sh #
# hplusAnalyzer.sh runs runMe.sh      #
# runMe.sh runs hplusAnalyzer.C       #
#                                     #
#//////////////////////////////////////

#------------------------------------------------
#create a directory where all the outputs will be
#stored, for different merged ntuple input files
#------------------------------------------------
outcond="outcond_29April17"

mkdir $outcond
cp hplusCond.sub $outcond
cp hplusAnalyzer.sh $outcond
cp mergedNtupleT2.txt $outcond
cd $outcond

#------------------------------------------------
#read the file, where paths of ntuples are stored
#do not put empty lines in mergedNtupleT2.txt
#------------------------------------------------
echo "make sure that you have coppied voms \n"
echo "certificate e.g. x509up_u93032 from /tmp\n"
echo "to /afs/cern.ch/user/r/rverma/ \n "

cat mergedNtupleT2.txt | while read ntupleT2Path
do
  #----------------------------------------------
  #print T2Paths of ntuple, on terminal
  #split the T2Paths of ntuples, into an array
  #get the  second last entry of the array
  #remove .root, from the input ntuple
  #----------------------------------------------
  echo " "
  echo -e "\033[01;32m input ntuple: \033[00m" $ntupleT2Path
  IFS='/ ' read -r -a array <<< "$ntupleT2Path"
  len=${#array[@]}
  sec_last=`expr $len - 1`
  #sec_last=`expr $len`
  ntuple=${array[$sec_last]}
  #echo $ntuple
  iFile=${ntuple/.root/""}
 
  #----------------------------------------------
  #copy condor scripts to each input ntuple dir
  #replace hplusCond.sub arguments, as per input
  #submit the condor jobs, for each ntuple
  #----------------------------------------------
  
  mkdir $iFile
  cp hplusCond.sub $iFile
  cp hplusAnalyzer.sh $iFile
  cd $iFile 
  sed -i "s:FNAME:$ntupleT2Path:g" hplusCond.sub
  sed -i "s:OUTPUTFILE:$iFile:g" hplusCond.sub
  sed -i "s:OUTPUTDIR:$iFile:g" hplusCond.sub
  condor_submit hplusCond.sub
  cd ../
done
