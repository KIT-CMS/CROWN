#!/bin/sh
action(){

    _addpy() {
        [ ! -z "$1" ] && export PYTHONPATH="$1:${PYTHONPATH}" && echo "Add $1 to PYTHONPATH"
    }
    _addbin() {
        [ ! -z "$1" ] && export PATH="$1:${PATH}" && echo "Add $1 to PATH"
    }

    echo "------------------------------------------"
    echo " | USER = ${USER}"
    echo " | HOSTNAME = $(hostname)"
    echo " | ANA_NAME = ${ANA_NAME}"
    echo " | ENV_NAME = ${ENV_NAME}"
    echo " | TAG = ${TAG}"
    echo " | TARBALL_PATH = ${TARBALL_PATH}"

    # Setup variables
    SPAWNPOINT=$(pwd)

    
    if [[ ! -f "/cvmfs/etp.kit.edu/LAW_envs/${ENV_NAME}/bin/activate" ]]; then
        if [[ -z "${TARBALL_ENV_PATH}" ]]; then
            echo "Environment ${ENV_NAME} is neither available in cvmfs nor as a tarball. Aborting job."
            return 1
        fi
        ENV_FROM_TAR="True"
    fi

    # on ETP ressources, this is not needed, since we can use a docker image containing the conda environment
    # on other systems, these changes here might be needed

    if [[ ! -z "${ENV_FROM_TAR}" ]]; then
        ENV_PATH=${SPAWNPOINT}/miniconda/envs/${ENV_NAME}
        echo " | ENV_PATH = $ENV_PATH"
        echo " | TARBALL_ENV_PATH = ${TARBALL_ENV_PATH}"
    else
        ENV_PATH=/cvmfs/etp.kit.edu/LAW_envs/${ENV_NAME}
        echo " | ENV_PATH = ${ENV_PATH}"
    fi
    echo "------------------------------------------"

    # copy and untar process (and environment if necessary) 
    if [[ -z "${ENV_FROM_TAR}" ]]; then
        # Activate environment from cvmfs
        source ${ENV_PATH}/bin/activate
        echo "xrdcp ${TARBALL_PATH} ${SPAWNPOINT}"
        xrdcp ${TARBALL_PATH} ${SPAWNPOINT}
    else
        (
            source /cvmfs/grid.cern.ch/umd-c7ui-latest/etc/profile.d/setup-c7-ui-example.sh
            echo "xrdcp ${TARBALL_PATH} ${SPAWNPOINT}"
            xrdcp ${TARBALL_PATH} ${SPAWNPOINT}
            echo "xrdcp ${TARBALL_ENV_PATH} ${SPAWNPOINT}"
            xrdcp ${TARBALL_ENV_PATH} ${SPAWNPOINT}
        )
        mkdir -p ${ENV_PATH}
        tar -xzf ${ENV_NAME}_env.tar.gz -C ${ENV_PATH} && rm ${ENV_NAME}_env.tar.gz
        # Activate environment from tarball
        source ${ENV_PATH}/bin/activate
        conda-unpack
    fi

    tar -xzf processor.tar.gz && rm processor.tar.gz

    # (
    #     source /cvmfs/grid.cern.ch/umd-c7ui-latest/etc/profile.d/setup-c7-ui-example.sh
    # )

    # add python modules from conda to pythonpath
    # _addpy "$(find $ENV_PATH/lib/python* -name site-packages)"

    # # add law to path
    # # law
    _addpy "${SPAWNPOINT}/law"
    _addbin "${SPAWNPOINT}/law/bin"

    # tasks
    _addpy "${SPAWNPOINT}/processor"
    _addpy "${SPAWNPOINT}/processor/tasks"

    # setup law variables
    export LAW_HOME="${SPAWNPOINT}/.law"
    export LAW_CONFIG_FILE="${SPAWNPOINT}/lawluigi_configs/${ANA_NAME}_law.cfg"
    export LUIGI_CONFIG_PATH="${SPAWNPOINT}/lawluigi_configs/${ANA_NAME}_luigi.cfg"

    # which law
}

action