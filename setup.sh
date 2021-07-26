#!/usr/bin/env bash

############################################################################################
# This script setups all dependencies necessary for making law executable                  #
############################################################################################

action() {
    # determine the directy of this file
    if [ ! -z "$ZSH_VERSION" ]; then
        local this_file="${(%):-%x}"
    else
        local this_file="${BASH_SOURCE[0]}"
    fi

    local base="$( cd "$( dirname "$this_file" )" && pwd )"
    local parent="$( dirname "$base" )"

    _addpy() {
        [ ! -z "$1" ] && export PYTHONPATH="$1:$PYTHONPATH"
    }

    _addbin() {
        [ ! -z "$1" ] && export PATH="$1:$PATH"
    }


    export LAW_HOME="$base/.law"
    export LAW_CONFIG_FILE="$base/law.cfg"
    export LUIGI_CONFIG_PATH="$base/luigi.cfg"
    export ANALYSIS_PATH="$base"
    export ANALYSIS_DATA_PATH="$ANALYSIS_PATH/data"


    # luigi
    _addpy "$base/luigi"
    _addbin "$base/luigi/bin"

    # six
    _addpy "$base/six"

    # law
    _addpy "$base/law"
    _addbin "$base/law/bin"
    source "$( law completion )"

    source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc9-opt/setup.sh
}
action "$@"