#!/bin/sh
action(){

    _addpy() {
        [ ! -z "$1" ] && export PYTHONPATH="$1:$PYTHONPATH"
    }
    _addbin() {
        [ ! -z "$1" ] && export PATH="$1:$PATH"
    }
    printenv
    SPAWNPOINT=$(pwd)
    conda_env="KingMaker"
    conda_env_path=$SPAWNPOINT/miniconda/envs/$conda_env
    # untar tarball
    tar -xzf processor*.tar.gz
    ls -la
    rm processor*.tar.gz
    # setup conda
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh
    bash Miniconda3-py39_4.10.3-Linux-x86_64.sh -b -p $SPAWNPOINT/miniconda
    mkdir -p $conda_env_path
    ls -la
    ls -la $conda_env_path
    tar -xzf tarballs/conda.tar.gz -C $conda_env_path
    source $conda_env_path/bin/activate
    conda-unpack
    # source $conda_env_path/bin/deactivate
    # conda activate $conda_env
    # add law to path
    # law
    _addpy "$SPAWNPOINT/law"
    _addbin "$SPAWNPOINT/law/bin"

    # tasks
    _addpy "$SPAWNPOINT/processor"
    _addpy "$SPAWNPOINT/processor/tasks"

    # setup law
    export LAW_HOME="$SPAWNPOINT/.law"
    export LAW_CONFIG_FILE="$SPAWNPOINT/law.cfg"
    export LUIGI_CONFIG_PATH="$SPAWNPOINT/luigi.cfg"

    which law
}

action