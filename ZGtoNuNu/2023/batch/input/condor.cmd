universe = vanilla
getenv = True
+JobFlavour = "longlunch"
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
success_exit_code = 0
max_retries = 0
arguments      = $(Process)
executable     = /eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/input/condor.sh
output         = ../../../../user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/logs/condor-$(Cluster)_$(Process).out
error          = ../../../../user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/logs/condor-$(Cluster)_$(Process).err
log            = ../../../../user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/logs/condor-$(Cluster)_$(Process).log
queue 158
