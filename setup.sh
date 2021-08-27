############################################################################################
# This script setups all dependencies necessary for making law executable                  #
############################################################################################
action() {
    miniconda="Miniconda3-py39_4.10.3-Linux-x86_64"

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
        echo "conda could not be found, installing conda ..."
        echo "More information can be found in"
        echo "https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html"
        echo "Follow the installer instructions .."
        curl -O https://repo.anaconda.com/miniconda/$miniconda.sh
        bash $miniconda.sh
        rm $miniconda.sh
    fi

    # check if Conda env is running
    echo "Setting up conda env"
    #
    if [ $CONDA_DEFAULT_ENV!="KingMaker" ]; then
        if conda env list | grep -q "KingMaker"; then
            echo  "KingMakter env found, activating..."
            conda activate KingMaker
        else
            echo "Creating KingMaker env ..."
            conda env create -f environment.yml
            echo  "KingMakter env found, activating..."
            conda activate KingMaker
            conda pack -n KingMaker --output tarballs/conda.tar.gz
        fi
    fi
    # since we need a conda tarball for the remote jobs, create it if it doesn't exist
    if [ ! -f "tarballs/conda.tar.gz" ]; then
        echo "Creating conda.tar.gz"
        conda pack -n KingMaker --output tarballs/conda.tar.gz
    fi

    if [ ! -d CROWN ]; then
        echo "Setup CROWN ..."
        git clone git@github.com:KIT-CMS/CROWN
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
