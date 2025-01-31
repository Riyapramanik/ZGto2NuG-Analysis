# Bamboo-tutorial of H -> ZZ* -> 4l using Run3 data
You can find out more about bamboo in [the UserGuide](https://bamboo-hep.readthedocs.io/en/latest/index.html). Also feel free to report any issue you encounter in [~bamboo](https://mattermost.web.cern.ch/cms-exp/channels/bamboo) channel on the CERN mattermost, or on [Gitlab](https://gitlab.cern.ch/cp3-cms/bamboo/-/issues).

- Slides on indico: https://indico.cern.ch/event/1445287/contributions/6210844/

## Table of Contents
- [Getting started: Set up your workflow](#getting-started-set-up-your-workflow)
- [Bamboo Installation (only the 1st time)](#bamboo-installation-only-1st-time)
- [Environment setup (Always)](#environment-setup-always)
- [Test your bamboo setup](#test-your-bamboo-setup)
- [How to run the bamboo tutorial ?](#how-to-run-the-bamboo-tutorial-)
- [Useful links](#useful-links)
- [Troule-shooting](#trouble-shooting)

## Getting started: Set up your workflow:

### Bamboo Installation (only 1st time):
- Set up your virtual environment in your ``~/.bashrc`` file:
```bash
alias voms="voms-proxy-init -voms cms -rfc -valid 192:00"

function bambooenv(){
    source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh #( lxplus el9)
    venvdir=$HOME/bamboovenv$LCG_VERSION
    if [ ! -d "$venvdir" ]; then
        python -m venv $venvdir
    fi
    source $venvdir/bin/activate
    export EXTRA_CLING_ARGS=-O2
}
```
- And steup computing environment configuration file in your ``~/.config/bamboorc`` add:
```
[batch]
backend = htcondor  ; or slurm

[htcondor]
jobflavour = "longlunch"

[das]
sitename = T2_CH_CERN
storageroot = /eos/cms
xrootdredirector = cms-xrd-global.cern.ch
checklocalfiles = True
```
- Install bamboo and some other required packages:
```bash
#make a virtualenv
source ~/.bashrc # or in a new shell type bambooenv
bambooenv 

mkdir bamboodev
cd bamboodev

# clone and install bamboo
git clone -o upstream https://gitlab.cern.ch/cp3-cms/bamboo.git
pip install ./bamboo

# clone and install plotIt
git clone -o upstream https://github.com/cp3-llbb/plotIt.git
mkdir build-plotit
cd build-plotit
cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../plotIt
make -j2 install
cd -

# The last two commands above; < pip install ./bamboo > and < make -j2 install >
# need to be run again, every time you upgrade your LCG working version!

# To use scalefactors and weights in the new CMS JSON format, the correctionlib package should be installed with the following cmd
# you can also ignore torch and sphinx pip errors !
pip install --no-binary=correctionlib correctionlib

# To use the calculators modules for jet and MET corrections and systematic variations
pip install git+https://gitlab.cern.ch/cp3-cms/CMSJMECalculators.git

# You will also need the python implementation "pyplotit" of "plotIt", which can be installed with
pip install git+https://gitlab.cern.ch/cp3-cms/pyplotit.git
```

### **Environment Setup (Always)**:
- Every time you want to setup your bamboo enviroment, what you simply need to do:
```bash
voms-proxy-init -voms cms -rfc -valid 192:00
bambooenv
```

### Test your bamboo setup:
```bash
cd bamboodev/bambooo
bambooRun -m examples/nanozmumu.py:NanoZMuMu \
          examples/test1.yml \
          -o test1
```

## How to run the bamboo tutorial ?
To run H-> ZZ* ->4l bamboo tutorial
```bash
bambooRun -m h4l_tutorial.py:Hto4lPlotter \
          config/analysis.yml \
          --samples config/h4l_2023samples.yml \
          --maxFiles=1 \
          -o test \
```
Or you can simply use `run_bamboo.sh` script. 
You can alternatively add: 
- `--maxFiles=1`: Process one root file for each dataset in the analysis configuration.
- `--distributed=driver`: Submit jobs to HTCondor or SLURM depending on the `[batch]` configuration in `~/.config/bamboorc`.
- `--distributed=finalize`: When all initially failed jobs have been rerun, use bambooRun with `finalize` (and the same options as the original submission) to merge outputs and run the postprocessing step. For more details, see [here](https://bamboo-hep.readthedocs.io/en/latest/advanced.html#distributed-rdataframe) and [here](https://bamboo-hep.readthedocs.io/en/latest/recipes.html#dealing-with-failed-batch-jobs).
- For more options, see [the user guide](https://bamboo-hep.readthedocs.io/en/latest/userguide.html#running-bamboorun).

# Useful links: 
- [CERN Batch Service User Guide](https://batchdocs.web.cern.ch/index.html)
- [HTCondor Version 23.8.1 Manual](https://htcondor.readthedocs.io/en/latest/)
- [Slurm workload manager Version 24.05](https://slurm.schedmd.com/quickstart_admin.html)

# Trouble-shooting:
- HTCondor: Using x509 proxy, see [here](https://batchdocs.web.cern.ch/tutorial/exercise2e_proxy.html#troubleshooting). After `voms-proxy-init -voms cms -rfc -valid 192:00` and before you submit your jobs with `--distributed=driver` do:
```bash
cp $(voms-proxy-info -p) ~/.x509_proxy
export X509_USER_PROXY=$(realpath ~/.x509_proxy)
```


