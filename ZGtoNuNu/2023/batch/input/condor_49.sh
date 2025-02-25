#!/usr/bin/env bash


function move_files {
 for file in *.root; do
   echo "Moving $file to /eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/output/49/"
   mv $file /eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/output/49/
 done
}

bambooRun --module=/eos/user/r/rpramani/bamboodev/myAnalysis/ZGto2NuG1Jets.py:ZGto2NuGPlotter --distributed=worker --anaConfig=/eos/user/r/rpramani/bamboodev/myAnalysis/config/analysis.yml --plotIt plotIt --samples config/ZGto2NuG1Jets.yml --input=/eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/infiles/ZGto2NuG-1Jets_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_in_10.txt --output=ZGto2NuG-1Jets_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.root --tree=Events --sample=ZGto2NuG-1Jets_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8 && move_files