#!/usr/bin/env bash


function move_files {
 for file in *.root; do
   echo "Moving $file to /eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/output/1/"
   mv $file /eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/batch/output/1/
 done
}

bambooRun --module=/eos/user/r/rpramani/bamboodev/myAnalysis/ZGto2NuG1Jets.py:ZGto2NuGPlotter --distributed=worker --anaConfig=/eos/user/r/rpramani/bamboodev/myAnalysis/config/analysis.yml --plotIt plotIt --samples config/ZGto2NuG1Jets.yml --input=/eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/infiles/EGamma_Run2023C_in_1.txt --output=EGamma_Run2023C.root --tree=Events --certifiedLumiFile=/eos/user/r/rpramani/bamboodev/myAnalysis/ZGtoNuNu/2023/Cert_Collisions2022_eraD_357538_357900_Golden.json --runRange=357735,357735 --sample=EGamma_Run2023C && move_files