############################################################################################
# This script setups all dependencies necessary for making law executable                  #
############################################################################################

action() {

    # Check if current machine is an etp portal machine.
    PORTAL_LIST=("bms1.etp.kit.edu" "bms2.etp.kit.edu" "bms3.etp.kit.edu" "portal1.etp.kit.edu")
    CURRENT_HOST=$(hostname --long)
    if [[ ! " ${PORTAL_LIST[*]} " =~ " ${CURRENT_HOST} " ]]; then  
        echo "Current host (${CURRENT_HOST}) not in list of allowed machines:"
        printf '%s\n' "${PORTAL_LIST[@]}"
        return 1
    else
        echo "Running on ${CURRENT_HOST}."
    fi

    ANA_NAME_GIVEN=$1

    #list of available analyses
    ANA_LIST=("KingMaker" "ML_LAW")
    #Determine analysis to be used. Default is first in list.
    if [[ -z "${ANA_NAME_GIVEN}" ]]; then
        echo "No analysis chosen. Using default analysis ${ANA_LIST[0]}."
        export ANA_NAME=${ANA_LIST[0]}
    else
        #Check if given analysis is in list 
        if [[ ! " ${ANA_LIST[*]} " =~ " ${ANA_NAME_GIVEN} " ]]; then 
            echo "Not a valid name. Allowed choices are:"
            printf '%s\n' "${ANA_LIST[@]}"
            return 1
        else
            echo "Using ${ANA_NAME_GIVEN} analysis." 
            export ANA_NAME="${ANA_NAME_GIVEN}"
        fi
    fi
    
    #Determine which environment to use based on the luigi.cfg file
    export ENV_NAME=$(grep "ENV_NAME" lawluigi_configs/${ANA_NAME}_luigi.cfg | sed 's@ENV_NAME\s*=\s*\([^\s]*\)\s*@\1@g')
    echo "Using ${ENV_NAME} environment."

    # Miniconda version used for all environments
    MINICONDA_VERSION="Miniconda3-py39_4.10.3-Linux-x86_64"
    # determine the directory of this file
    if [ ! -z "${ZSH_VERSION}" ]; then
        local THIS_FILE="${(%):-%x}"
    else
        local THIS_FILE="${BASH_SOURCE[0]}"
    fi

    local BASE_DIR="$( cd "$( dirname "${THIS_FILE}" )" && pwd )"

    _addpy() {
        [ ! -z "$1" ] && export PYTHONPATH="$1:${PYTHONPATH}"
    }

    _addbin() {
        [ ! -z "$1" ] && export PATH="$1:${PATH}"
    }

    #Check if necessary environment is present in cvmfs
    if [[ ! -f "/cvmfs/etp.kit.edu/LAW_envs/${ENV_NAME}/bin/activate" ]]; then
        #If not present
        #install conda in necessary
        echo "${ENV_NAME} environment not found in cvmfs. Using conda."
        if ! command -v conda &> /dev/null
        then
            if [ ! -f "miniconda/bin/activate" ]; then
                echo "conda could not be found, installing conda ..."
                echo "More information can be found in"
                echo "https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html"
                curl -O https://repo.anaconda.com/miniconda/${MINICONDA_VERSION}.sh
                bash ${MINICONDA_VERSION}.sh -b -p miniconda
                rm -f ${MINICONDA_VERSION}.sh
            fi
            source miniconda/bin/activate base
        fi
        # check if correct Conda env is running
        echo "Setting up conda.."
        if [ "${CONDA_DEFAULT_ENV}" != "${ENV_NAME}" ]; then
            if [ -d "miniconda/envs/${ENV_NAME}" ]; then
                echo  "${ENV_NAME} env found, activating..."
                conda activate ${ENV_NAME}
            else
                echo "Creating ${ENV_NAME} env from conda_environments/${ENV_NAME}_env.yml..."
                if [[ ! -f "conda_environments/${ENV_NAME}_env.yml" ]]; then
                    echo "conda_environments/${ENV_NAME}_env.yml not found. Unable to create environment."
                    return 1
                fi
                conda env create -f conda_environments/${ENV_NAME}_env.yml -n ${ENV_NAME}
                echo  "${ENV_NAME} env built, activating..."
                conda activate ${ENV_NAME}
            fi
        fi
        # create conda tarball if env not in cvmfs and it if it doesn't exist already
        if [ ! -f "tarballs/${ENV_NAME}_env.tar.gz" ]; then
            echo "Creating ${ENV_NAME}_env.tar.gz"
            mkdir -p "tarballs"
            conda pack -n ${ENV_NAME} --output tarballs/${ENV_NAME}_env.tar.gz
        fi
        export CVMFS_ENV_PRESENT="False"
    else
        #If present 
        source /cvmfs/etp.kit.edu/LAW_envs/${ENV_NAME}/bin/activate
        export CVMFS_ENV_PRESENT="True"
    fi
    #CVMFS_ENV_PRESENT can be used by the processor/framework.py to determine whether the environment is present in cvmfs

    #Set up other dependencies based on analysis
    ############################################
    case ${ANA_NAME} in
        KingMaker)
            echo "Setup CROWN ..."
            if [ ! -d CROWN ]; then
                git clone git@github.com:KIT-CMS/CROWN
            fi
            ;;
        *)
            ;;
    esac
    ############################################

    # Check is law was cloned, and set it up if not
    if [ -z "$(ls -A law)" ]; then
        git submodule update --init --recursive
    fi

    echo "Setting up Luigi/Law ..."
    export LAW_HOME="${BASE_DIR}/.law"
    export LAW_CONFIG_FILE="${BASE_DIR}/lawluigi_configs/${ANA_NAME}_law.cfg"
    export LUIGI_CONFIG_PATH="${BASE_DIR}/lawluigi_configs/${ANA_NAME}_luigi.cfg"
    export ANALYSIS_PATH="${BASE_DIR}"
    export ANALYSIS_DATA_PATH="${ANALYSIS_PATH}/data"

    # law
    _addpy "${BASE_DIR}/law"
    _addbin "${BASE_DIR}/law/bin"
    source "$( law completion )"

    # tasks
    _addpy "${BASE_DIR}/processor"
    _addpy "${BASE_DIR}/processor/tasks"

    # add voms proxy path
    export X509_USER_PROXY=$(voms-proxy-info -path)

    # start a luidigd scheduler if there is one already running
    if [ -z "$(pgrep -f luigid)" ]; then
        echo "Starting Luigi scheduler..."
        luigid --background --logdir logs --state-path luigid_state.pickle
    fi
}
action "$@"