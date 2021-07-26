#!/bin/sh

action(){
    SPAWNPOINT=$(pwd)
    # source nightly cvmfs
    source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc9-opt/setup.sh
    # source grid environment
    source /cvmfs/grid.cern.ch/centos7-wn-4.0.5-1_umd4v1/etc/profile.d/setup-c7-wn-example.sh

    # untar tarball
    tar -xzf crown*.tar.gz
    rm crown*.tar.gz

    # setup law
    export LAW_HOME="$PWD/.law"
    export LAW_CONFIG_FILE="$PWD/law.cfg"
    export LUIGI_CONFIG_PATH="$PWD/luigi.cfg"

    export PATH="$PWD/law/bin:$PWD/luigi/bin:$PATH"
    export PYTHONPATH="$PWD/enum34-1.1.10:$PWD/law:$PWD/luigi:$PWD/six:$PWD:$PYTHONPATH"

}

action