# REFERENCE 
# http://research.cs.wisc.edu/htcondor/manual/
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookXrootdService

#---------------------------------------------
# voms-proxy-init --voms cms, to use T2 files
# hplusAnalyzer.sh is run by condor_submit
# Arguments are passed to hplusAnalyzer.sh
#---------------------------------------------

Universe              = vanilla
x509userproxy         = /tmp/x509up_u93032
#use_x509userproxy    = true
Executable            = hplusAnalyzer.sh
Arguments             = FNAME OUTPUTFILE OUTPUTDIR $(Process) $(Cluster)

#---------------------------------------------
# hplusAnalyzer.C output are stored in .stdout
# all kind of errors are stores in .stderr
# CPU usage, RAM etc are stored in .log
#---------------------------------------------

Output                = condor_out_$(Process)_$(Cluster).stdout
Error                 = condor_out_$(Process)_$(Cluster).stderr
Log                   = condor_out_$(Process)_$(Cluster).log
Notification          = Error

#---------------------------------------------
# transfer output from remote machine to T2
# put the condor batch submission in queue 
#---------------------------------------------

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
Queue

