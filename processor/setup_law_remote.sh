#!/bin/sh
action(){

    _addpy() {
        [ ! -z "$1" ] && export PYTHONPATH="$1:${PYTHONPATH}" && echo "Add $1 to PYTHONPATH"
    }
    _addbin() {
        [ ! -z "$1" ] && export PATH="$1:${PATH}" && echo "Add $1 to PATH"
    }

    # Set USER as local USER
    export USER={{USER}}
    echo "------------------------------------------"
    echo " | USER = ${USER}"
    echo " | HOSTNAME = $(hostname)"
    echo " | ANA_NAME = {{ANA_NAME}}"
    echo " | ENV_NAME = {{ENV_NAME}}"
    echo " | TAG = {{TAG}}"
    echo " | USE_CVMFS = {{USE_CVMFS}}"
    echo " | TARBALL_PATH = {{TARBALL_PATH}}"
    SPAWNPOINT=$(pwd)
    if [[ "{{USE_CVMFS}}" == "True" ]]; then
        ENV_PATH=/cvmfs/etp.kit.edu/LAW_envs/{{ENV_NAME}}
        echo " | ENV_PATH = ${ENV_PATH}"
    else
        ENV_PATH=${SPAWNPOINT}/miniconda/envs/{{ENV_NAME}}
        echo " | ENV_PATH = $ENV_PATH"
        echo " | TARBALL_ENV_PATH = {{TARBALL_ENV_PATH}}"
    fi
    echo "------------------------------------------"

    # copy and untar process (and environment if necessary)
    if [[ "{{USE_CVMFS}}" == "True" ]]; then
        # Activate environment from cvmfs
        source ${ENV_PATH}/bin/activate
        echo "xrdcp {{TARBALL_PATH}} ${SPAWNPOINT}"
        xrdcp {{TARBALL_PATH}} ${SPAWNPOINT}
    else
        # Copy tarballs 
        (
            source /cvmfs/etp.kit.edu/LAW_envs/KingMaker/bin/activate
            echo "xrdcp {{TARBALL_PATH}} ${SPAWNPOINT}"
            xrdcp {{TARBALL_PATH}} ${SPAWNPOINT}
            echo "xrdcp {{TARBALL_ENV_PATH}} ${SPAWNPOINT}"
            xrdcp {{TARBALL_ENV_PATH}} ${SPAWNPOINT}
        )
        mkdir -p ${ENV_PATH}
        tar -xzf {{ENV_NAME}}.tar.gz -C ${ENV_PATH} && rm {{ENV_NAME}}.tar.gz
        # Activate environment from tarball
        source ${ENV_PATH}/bin/activate
        conda-unpack
    fi

    tar -xzf processor.tar.gz && rm processor.tar.gz

    # # add law to path
    # # law
    _addpy "${SPAWNPOINT}/law"
    _addbin "${SPAWNPOINT}/law/bin"

    # tasks
    _addpy "${SPAWNPOINT}/processor"
    _addpy "${SPAWNPOINT}/processor/tasks"

    # setup law variables
    export LAW_HOME="${SPAWNPOINT}/.law"
    export LAW_CONFIG_FILE="${SPAWNPOINT}/lawluigi_configs/{{ANA_NAME}}_law.cfg"
    export LUIGI_CONFIG_PATH="${SPAWNPOINT}/lawluigi_configs/{{ANA_NAME}}_luigi.cfg"

    export LOCAL_TIMESTAMP="{{LOCAL_TIMESTAMP}}"

}

action
