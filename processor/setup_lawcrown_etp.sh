#!/bin/sh
action(){

    _addpy() {
        [ ! -z "$1" ] && export PYTHONPATH="$1:$PYTHONPATH"
    }
    _addbin() {
        [ ! -z "$1" ] && export PATH="$1:$PATH"
    }

    SPAWNPOINT=$(pwd)
    conda_env="KingMaker"
    while read $tarball_path; do
        echo "Getting conda tarball from $conda_tarball_path"
    done < "tarball_path*.txt"
    echo "------------------------------------------"
    echo " | conda_env = $conda_env"
    echo " | JOB_TARBALL = $tarball_path"
    echo "------------------------------------------"
    # initiate conda
    source /etc/profile
    conda activate $conda_env
    # first get the conda tarball path from the input textfile:
    (
        # source /cvmfs/grid.cern.ch/umd-c7ui-latest/etc/profile.d/setup-c7-ui-example.sh
        gfal-copy $tarball_path $SPAWNPOINT
    )

    tar -xzf processor*.tar.gz
    rm processor*.tar.gz
    # # add law to path
    # # law
    _addpy "$SPAWNPOINT/law"
    _addbin "$SPAWNPOINT/law/bin"

    # tasks
    _addpy "$SPAWNPOINT/processor"
    _addpy "$SPAWNPOINT/processor/tasks"

    # setup law
    export LAW_HOME="$SPAWNPOINT/.law"
    export LAW_CONFIG_FILE="$SPAWNPOINT/law.cfg"
    export LUIGI_CONFIG_PATH="$SPAWNPOINT/luigi.cfg"
}

action