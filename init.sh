#!/bin/bash

action() {

    # --- Define list of available analyses ---
    ANALYSIS_LIST=("tau" "earlyrun3" "whtautau" "boosted_h_tautau" "s" "xyh_bbtautau" "haa")

    # --- Define defaults ---
    DEFAULT_CROWN_ANALYSIS=""
    DEFAULT_CONTAINER="docker://tvoigtlaender/kingmaker_standalone:V1.4j"
    DEFAULT_DRY_RUN=false
    CROWN_ANALYSIS=${DEFAULT_CROWN_ANALYSIS}
    CONTAINER=${DEFAULT_CONTAINER}
    DRY_RUN=${DEFAULT_DRY_RUN}

    # --- Parse arguments ---
    while [[ $# -gt 0 ]]; do
        case $1 in
            -a|--analysis)
                if [[ -z "$2" || "$2" == -* ]]; then
                    echo "Error: --analysis requires a non-option argument."
                    return 1
                fi
                CROWN_ANALYSIS="$2"
                shift 2
                ;;
            -c|--container)
                if [[ -z "$2" || "$2" == -* ]]; then
                    echo "Error: --container requires a non-option argument."
                    return 1
                fi
                CONTAINER="$2"
                shift 2
                ;;
            -d|--dry-run)
                DRY_RUN=true
                shift 1
                ;;
            -l|--list)
                echo "Available CROWN analyses:"
                echo "-------------------"
                for workflow in "${ANALYSIS_LIST[@]}"; do
                    echo "* ${workflow}"
                done
                return 1
                ;;
            -h|--help)
                echo "Usage: source setup.sh [options]"
                echo ""
                echo "Options:"
                echo "  -a, --analysis ANALYSIS    The CROWN analysis configs to use"
                echo "  -l, --list                 List available CROWN analyses"
                echo "  -c, --container CONTAINER  The container to use as environment"
                echo "                             [default: ${DEFAULT_CONTAINER}]"
                echo "  -h, --help                 Show this help message"
                echo ""
                return 1
                ;;
            *)
                echo "Error: Unknown option $1"
                echo "Use --help to see available options"
                return 1
                ;;
        esac
    done
    if [[ $? -eq "1" ]]; then
        return 1
    fi

    # --- Local Path Setup ---
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" &>/dev/null && pwd 2>/dev/null)"
    ANALYSIS_PATH="${SCRIPT_DIR}"

    # --- Clone Analysis Logic (if second arg is provided) ---
    if [ -n "${CROWN_ANALYSIS}" ]; then
        case "${CROWN_ANALYSIS}" in
            tau)              REPO="git@github.com:KIT-CMS/TauAnalysis-CROWN.git" ;;
            earlyrun3)        REPO="git@github.com:KIT-CMS/earlyRun3Analysis-CROWN" ;;
            whtautau)         REPO="git@github.com:KIT-CMS/WHTauTauAnalysis-CROWN.git" ;;
            boosted_h_tautau) REPO="git@github.com:KIT-CMS/BoostedHiggsTauTauAnalysis-CROWN.git" ;;
            s)                REPO="git@github.com:nfaltermann/CROWNs.git" ;;
            xyh_bbtautau)     REPO="git@github.com:KIT-CMS/XYHBBTauTauAnalysis-CROWN.git" ;;
            haa)              REPO="git@github.com:KIT-CMS/HaaAnalysis-CROWN.git" ;;
            *)                echo "Error: '${CROWN_ANALYSIS}' is not a valid analysis name."
                              echo "See --list argument for available analyses."
                              return 1
            ;;
        esac

        if [ -n "${REPO}" ] && [ ! -d "${SCRIPT_DIR}/analysis_configurations/${CROWN_ANALYSIS}" ]; then
            echo "--> Cloning analysis ${CROWN_ANALYSIS}..."
            git clone "${REPO}" "${SCRIPT_DIR}/analysis_configurations/${CROWN_ANALYSIS}"
        else
            echo "--> Analysis configuration '${CROWN_ANALYSIS}' already exists."
        fi
    fi

    # --- Dry Run Exit Point ---
    if [ "${DRY_RUN}" = true ]; then
        echo "--> Dry run complete. Analysis configs are ready."
        return 0 2>/dev/null || exit 0
    fi

    # --- VOMS Proxy Check ---
    if ! voms-proxy-info -exists -file "${X509_USER_PROXY}" >/dev/null 2>&1; then
        echo "Warning: No valid VOMS proxy found. Grid storage may be inaccessible."
    fi

    # --- Get the absolute path of the parent git repository ---
    # checks/git-status.sh fails if top level git directory is not mounted in
    GIT_ROOT=$(git -C "${SCRIPT_DIR}" rev-parse --show-superproject-working-tree)

    # --- Define the Internal Environment ---
    # This string is executed once the container starts
    INT_CMD="
        echo '--- Initializing Conda Environment ---';
        export ANALYSIS_PATH=${ANALYSIS_PATH};
        export CROWN_ANALYSIS=${CROWN_ANALYSIS};
        export CCACHE_DIR=${CROWN_ANALYSIS}/.cache/ccache;
        export CMAKE_GENERATOR='Unix Makefiles';
        export EXTRA_CLING_ARGS='-O2';
        export X509_USER_PROXY=${X509_USER_PROXY};
        echo '--- Container Ready. Analysis: ${CROWN_ANALYSIS} ---';
        bash --rcfile /etc/bashrc -i
    "

    # --- Execute Singularity ---
    echo "--> Launching Container: ${CONTAINER}"
    singularity exec -e \
        -B /etc/grid-security/certificates \
        -B "${GIT_ROOT}:${GIT_ROOT}" \
        -B "${HOME}:${HOME}" \
        "${CONTAINER}" \
        bash -c "${INT_CMD}"
}
action "$@"
