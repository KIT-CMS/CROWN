#!/bin/bash
pathadd() {
    if [[ ":$PATH:" != *":$1:"* ]]; then
        PATH="${PATH:+"$PATH:"}$1"
        export PATH
    fi
}
# get the directory of the script
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" &>/dev/null && pwd 2>/dev/null)"
distro=$(lsb_release -i | cut -f2)
os_version=$(lsb_release -r | cut -f2)
echo "Setting up CROWN for $distro Version $os_version"
# check if the distro is centos
if [[ "$distro" == "CentOS" ]]; then
    # if the first number of os_version is a 7, we are on centOS 7
    if [[ ${os_version:0:1} == "7" ]]; then # if uname -a | grep -E 'el7' -q
        # source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc11-opt/setup.sh
        # source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-clang12-opt/setup.sh
        # source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc11-dbg/setup.sh
        source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
    else
        echo "Unsupported CentOS version, exiting..."
        return 0
    fi
elif [[ "$distro" == "RedHatEnterprise" ]]; then
    if [[ ${os_version:0:1} == "8" ]]; then # elif uname -a | grep -E 'el8' -q
        source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos8-gcc11-opt/setup.sh
    else
        echo "Unsupported CentOS version, exiting..."
        return 0
    fi
elif [[ "$distro" == "Ubuntu" ]]; then
    if [[ ${os_version:0:2} == "20" ]]; then
        source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-ubuntu2004-gcc9-opt/setup.sh
    elif [[ ${os_version:0:2} == "22" ]]; then
        source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-ubuntu2204-gcc11-opt/setup.sh
    else
        echo "Unsupported Ubuntu version, exiting..."
        return 0
    fi
else
    echo "You are not running on CentOS or Ubuntu, exiting..."
    return 0
fi
# add ~/.local/bin to path if it is not already there
pathadd "${HOME}/.local/bin/"
# set the cmake generator to Ninja
# export CMAKE_GENERATOR="Ninja"
export CMAKE_GENERATOR="Unix Makefiles"

# clone a given analysis if an argument is given
if [ -z "$1" ]; then
    echo "No configuration clone"
else
    if [[ "$1" == "tau" && ! -d "${SCRIPT_DIR}/analysis_configurations/tau" ]]; then
        echo "Cloning analysis tau into ${SCRIPT_DIR}/analysis_configurations/tau"
        git clone git@github.com:KIT-CMS/TauAnalysis-CROWN.git "${SCRIPT_DIR}/analysis_configurations/tau"
    elif [[ "$1" == "earlyrun3" && ! -d "${SCRIPT_DIR}/analysis_configurations/earlyrun3" ]]; then
        echo "Cloning analysis earlyrun3 into ${SCRIPT_DIR}/analysis_configurations/earlyrun3"
        git clone https://github.com/khaosmos93/CROWN-config-earlyRun3.git "${SCRIPT_DIR}/analysis_configurations/earlyrun3"
    elif [[ "$1" == "whtautau" && ! -d "${SCRIPT_DIR}/analysis_configurations/whtautau" ]]; then
        echo "Cloning analysis whtautau into ${SCRIPT_DIR}/analysis_configurations/whtautau"
        git clone git@github.com:KIT-CMS/WHTauTauAnalysis-CROWN.git "${SCRIPT_DIR}/analysis_configurations/whtautau"
    elif [[ "$1" == "s" && ! -d "${SCRIPT_DIR}/analysis_configurations/s" ]]; then
        echo "Cloning analysis s-channel into ${SCRIPT_DIR}/analysis_configurations/s"
        git clone git@github.com:nfaltermann/CROWNs.git "${SCRIPT_DIR}/analysis_configurations/s"
    fi
fi
