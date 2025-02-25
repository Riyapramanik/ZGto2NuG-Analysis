#!/usr/bin/env bash


function move_files {
 for file in *.root; do
   echo "Moving $file to /eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/output/157/"
   mv $file /eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/output/157/
 done
}

bambooRun --module=/eos/user/r/rpramani/bamboodev/myAnalysis/ZGto2NuG1Jets.py:ZGto2NuGPlotter --distributed=worker --anaConfig=/eos/user/r/rpramani/bamboodev/myAnalysis/config/analysis.yml --plotIt plotIt --samples config/ZGto2NuG1Jets.yml --input=/eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/infiles/ZZ_TuneCP5_13p6TeV_pythia8_in_5.txt --output=ZZ_TuneCP5_13p6TeV_pythia8.root --tree=Events --sample=ZZ_TuneCP5_13p6TeV_pythia8 && move_files