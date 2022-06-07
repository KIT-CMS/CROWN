#!/bin/bash
pathadd() {
    if [[ ":$PATH:" != *":$1:"* ]]; then
        PATH="${PATH:+"$PATH:"}$1"; export PATH
    fi
}
# get the directory of the script
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";
if uname -a | grep -E 'el7' -q
then
    # source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-clang10-opt/setup.sh
    # source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc11-opt/setup.sh
    # source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-clang12-opt/setup.sh
    # source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc11-dbg/setup.sh
    ### Use a more permanent LCG stack, for now LCG 102
    source /cvmfs/sft.cern.ch/lcg/views/LCG_102rc1/x86_64-centos7-gcc11-opt/setup.sh
else
    echo "You are not running on CentOS7, things will propably break..."
fi
# add ~/.local/bin to path if it is not already there
pathadd "${HOME}/.local/bin/"
# set the cmake generator to Ninja
# export CMAKE_GENERATOR="Ninja"
export CMAKE_GENERATOR="Unix Makefiles"

# clone a given analysis if an argument is given
if [ -z "$1" ]
then
    echo "No configuration clone"
else
    if [[ "$1" == "tau" && ! -d "${SCRIPT_DIR}/analysis_configurations/tau" ]]
    then
        echo "Cloning analysis tau into ${SCRIPT_DIR}/analysis_configurations/tau"
        git clone git@github.com:KIT-CMS/TauAnalysis-CROWN.git ${SCRIPT_DIR}/analysis_configurations/tau
    elif [[ "$1" == "earlyrun3" && ! -d "${SCRIPT_DIR}/analysis_configurations/earlyrun3" ]]
    then
        echo "Cloning analysis earlyrun3 into ${SCRIPT_DIR}/analysis_configurations/earlyrun3"
        git clone https://github.com/khaosmos93/CROWN-config-earlyRun3.git ${SCRIPT_DIR}/analysis_configurations/earlyrun3
    fi
fi
