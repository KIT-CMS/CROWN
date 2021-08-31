############################################################################################
# This script setups all dependencies necessary for making law executable                  #
############################################################################################
action() {
    miniconda="Miniconda3-py39_4.10.3-Linux-x86_64"
    conda_env_name="KingMaker"
    fullhostname=$(hostname -f)
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

    #check if conda is installed
    if ! command -v conda &> /dev/null
    then
        if [ ! -f "miniconda/bin/activate" ]; then
            echo "conda could not be found, installing conda ..."
            echo "More information can be found in"
            echo "https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html"
            curl -O https://repo.anaconda.com/miniconda/$miniconda.sh
            bash $miniconda.sh -b -p miniconda
            rm -f $miniconda.sh
        fi
        source miniconda/bin/activate
    fi

    # check if Conda env is running
    echo "Setting up conda.."
    #
    if [ "$CONDA_DEFAULT_ENV" != "$conda_env_name" ]; then
        if [ -d "miniconda/envs/$conda_env_name" ]; then
            echo  "$conda_env_name env found, activating..."
            conda activate $conda_env_name
        else
            echo "Creating $conda_env_name env from environment.yml..."
            conda env create -f environment.yml
            echo  "KingMaker env found, activating..."
            conda activate $conda_env_name
        fi
    fi
    # since we need a conda tarball for the remote jobs, create it if it doesn't exist
    if [[ $fullhostname != *"etp"* ]]; then
        if [ ! -f "tarballs/conda.tar.gz" ]; then
            echo "Creating conda.tar.gz"
            mkdir -p "tarballs"
            conda pack -n $conda_env_name --output tarballs/conda.tar.gz
        fi
    fi

    echo "Setup CROWN ..."
    if [ ! -d CROWN ]; then
        git clone git@github.com:KIT-CMS/CROWN
    fi

    # Check is law was cloned, and set it up if not
    if [ -z "$(ls -A law)" ]; then
        git submodule update --init --recursive
    fi

    echo "Setting up Luigi/Law ..."
    export LAW_HOME="$base/.law"
    export LAW_CONFIG_FILE="$base/law.cfg"
    export LUIGI_CONFIG_PATH="$base/luigi.cfg"
    export ANALYSIS_PATH="$base"
    export ANALYSIS_DATA_PATH="$ANALYSIS_PATH/data"

    # law
    _addpy "$base/law"
    _addbin "$base/law/bin"
    source "$( law completion )"

    # tasks
    _addpy "$base/processor"
    _addpy "$base/processor/tasks"
}
action "$@"