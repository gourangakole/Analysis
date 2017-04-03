universe              = vanilla
Executable            = hplusAnalyzer.sh
Should_Transfer_Files = YES
WhenToTransferOutput  = ON_EXIT_OR_EVICT
Output                = condor_out_$(Process)_$(Cluster).stdout
Error                 = condor_out_$(Process)_$(Cluster).stderr
Log                   = condor_out_$(Process)_$(Cluster).log
Notification          = Error
Arguments             = FNAME OUTPUTFILE OUTPUTDIR $(Process) $(Cluster)
Queue 1

