#!/bin/bash -l

current_proxy=$(voms-proxy-info -p)
destination_proxy="$HOME/.x509_proxy"

# Check if the current proxy and destination are different
if [ "$current_proxy" != "$destination_proxy" ]; then
    cp "$current_proxy" "$destination_proxy"
fi

# Export the proxy path as an environment variable
export X509_USER_PROXY=$(realpath "$destination_proxy")
export PYTHONPATH=$PWD:$PYTHONPATH
echo "X509_USER_PROXY is set to $X509_USER_PROXY"

output_dir='ZGtoNuNu'

#samplesYAML='config/ZGto2NuG1Jets_test.yml'
samplesYAML='config/ZGto2NuG1Jets.yml'

run='--maxFiles=1'
#run='--distributed=driver'
#run='--onlyprepare'
#run='--onlypost'
#run='--distributed=finalize'

doSysts=false
era='2023' # choices '2023', '2023PBix', '' the latest will do all

plus_args=''
if [ "$era" != '' ]; then
    output_dir=${output_dir}/$era
    plus_args=' --era='$era
fi


CMD="bambooRun -m ZGto2NuG1Jets.py:Hto4lPlotter config/analysis.yml --samples $samplesYAML $plus_args $run -o $output_dir"
echo "Running command: $CMD"
eval $CMD

        
