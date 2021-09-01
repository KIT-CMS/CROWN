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

    echo "------------------------------------------"
    echo " | conda_env = $conda_env"
    echo " | tarball_path = $tarball_path"
    echo "------------------------------------------"

    # initiate conda
    source /etc/profile
    conda activate $conda_env

    # first get the conda tarball path from the input textfile:
    gfal-copy $tarball_path $SPAWNPOINT

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