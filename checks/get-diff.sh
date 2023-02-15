#!/bin/bash

basepath=$1
analysis=$2
resultfolder=$3
# throw exception if args are not provided
if [ -z "$analysis" ] || [ -z "$basepath" ] || [ -z "$resultfolder" ]; then
    echo "ERROR: analysis, basepath and resultfolder must be provided"
    exit 1
fi
# create the result folder if it does not exist
if [ ! -d "$resultfolder" ]; then
    mkdir -p $resultfolder
fi
# go to the basepath
cd $basepath
base_diff=$resultfolder/base_diff.patch
commit_hash=$resultfolder/base_commit_hash.txt
touch $base_diff
for next in $(git ls-files --others --exclude-standard); do
    git --no-pager diff --no-index /dev/null $next >> $base_diff;
done
git rev-parse HEAD > $commit_hash

# cd into the analysis directory
analysispath=$basepath/analysis_configurations/$analysis
if [ -d "$analysispath" ]; then
    cd $analysispath
    analysis_diff=$resultfolder/analysis_diff.patch
    analysis_commit=$resultfolder/analysis_commit_hash.txt
    analysis_name=$resultfolder/analysis_name.txt
    touch $analysis_diff
    for next in $(git ls-files --others --exclude-standard); do
        git --no-pager diff --no-index /dev/null $next >> $analysis_diff;
    done
    git rev-parse HEAD > $analysis_commit
    echo $analysis > $analysis_name
else
    analysiscommit="not_found"
fi
