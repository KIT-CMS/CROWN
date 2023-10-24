#!/bin/bash

basepath=$1
analysis=$2
# throw exception if args are not provided
if [ -z "$analysis" ] || [ -z "$basepath" ]; then
    echo "ERROR: analysis and basepath must be provided"
    exit 1
fi
# go to the basepath
cd "$basepath"
# get the status of the git repo and store it in a variable
if output=$(git status --porcelain) && [ -z "$output" ]; then
    crown_is_clean="true"
else
    # Uncommitted changes
    crown_is_clean="false"
fi
# get the current commit hash and store it in a variable
maingitcommit=$(git rev-parse HEAD)
# cd into the analysis directory
analysispath="$basepath/analysis_configurations/$analysis"
if [ -d "$analysispath" ]; then
    cd "$analysispath"
    # get the status of the git repo and store it in a variable
    if output=$(git status --porcelain) && [ -z "$output" ]; then
        analysis_is_clean="true"
    else
        # Uncommitted changes
        analysis_is_clean="false"
    fi
    # get the commit hash of the analysis and store it in a variable
    analysiscommit=$(git rev-parse HEAD)
else
    echo "ERROR: analysis $analysis does not exist"
    analysis_is_clean="false"
    analysiscommit="not_found"
fi
# echo the final results
echo "crown_is_clean=$crown_is_clean"
echo "commit_hash=$maingitcommit"
echo "analysis_is_clean=$analysis_is_clean"
echo "analysis_commit_hash=$analysiscommit"
