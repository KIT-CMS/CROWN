#!/bin/bash

# --- Usage Check ---
if [[ $# -lt 1 ]]; then
    echo "Usage: source $0 <container_path_or_url> [analysis_name]"
    echo "Example: source $0 docker://tvoigtlaender/kingmaker_standalone:V1.3 tau"
    return 1 2>/dev/null || exit 1
fi

CONTAINER_IMG="$1"
ANA_NAME="$2"

# --- Local Path Setup ---
# Get the directory of the script
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" &>/dev/null && pwd 2>/dev/null)"
ANALYSIS_PATH="${SCRIPT_DIR}"

# --- Clone Analysis Logic (if second arg is provided) ---
if [ -n "$2" ]; then
    case "$2" in
        tau) REPO="git@github.com:KIT-CMS/TauAnalysis-CROWN.git" ;;
        earlyrun3) REPO="git@github.com:KIT-CMS/earlyRun3Analysis-CROWN" ;;
        whtautau) REPO="git@github.com:KIT-CMS/WHTauTauAnalysis-CROWN.git" ;;
        boosted_h_tautau) REPO="git@github.com:KIT-CMS/BoostedHiggsTauTauAnalysis-CROWN.git" ;;
        s) REPO="git@github.com:nfaltermann/CROWNs.git" ;;
        xyh_bbtautau) REPO="git@github.com:KIT-CMS/XYHBBTauTauAnalysis-CROWN.git" ;;
        haa) REPO="git@github.com:KIT-CMS/HaaAnalysis-CROWN.git" ;;
    esac

    if [ -n "$REPO" ] && [ ! -d "${SCRIPT_DIR}/analysis_configurations/$2" ]; then
        echo "--> Cloning analysis $2..."
        git clone "$REPO" "${SCRIPT_DIR}/analysis_configurations/$2"
    fi
fi

# --- VOMS Proxy Check ---
if ! voms-proxy-info -exists -file "$X509_USER_PROXY" >/dev/null 2>&1; then
    echo "Warning: No valid VOMS proxy found. You may need to run 'voms-proxy-init' before the container."
fi

# --- Define the Internal Environment ---
# This string is executed once the container starts
INT_CMD="
    echo '--- Initializing Conda Environment ---';
    source /opt/conda/bin/activate env;

    export ANALYSIS_PATH=${ANALYSIS_PATH};
    export ANA_NAME=${ANA_NAME};
    export CMAKE_GENERATOR='Unix Makefiles';
    export EXTRA_CLING_ARGS='-O2';
    export X509_USER_PROXY=${X509_USER_PROXY};

    cd ${ANALYSIS_PATH};
    echo '--- Container Ready. Analysis: ${ANA_NAME} ---';
    exec bash -i
"

# --- Execute Singularity ---
echo "--> Launching Container: ${CONTAINER_IMG}"
singularity exec -e \
    -B /etc/grid-security/certificates \
    -B ${SCRIPT_DIR}:${SCRIPT_DIR} \
    -B ${HOME}:${HOME} \
    "${CONTAINER_IMG}" \
    bash -c "${INT_CMD}"
