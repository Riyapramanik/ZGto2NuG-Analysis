#!/usr/bin/env bash


function move_files {
 for file in *.root; do
   echo "Moving $file to /eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/output/80/"
   mv $file /eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/output/80/
 done
}

bambooRun --module=/eos/user/r/rpramani/bamboodev/myAnalysis/ZGto2NuG1Jets.py:ZGto2NuGPlotter --distributed=worker --anaConfig=/eos/user/r/rpramani/bamboodev/myAnalysis/config/analysis.yml --plotIt plotIt --samples config/ZGto2NuG1Jets.yml --input=/eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/infiles/TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8_in_3.txt --output=TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8.root --tree=Events --sample=TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8 && move_files