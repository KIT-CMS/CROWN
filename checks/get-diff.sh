#!/bin/bash

base_path=$1
analysis=$2
result_dir=$3
# throw exception if args are not provided
if [ -z "$analysis" ] || [ -z "$base_path" ] || [ -z "$result_dir" ]; then
    echo "ERROR: analysis, base_path and result_dir must be provided"
    exit 1
fi
# create the result folder if it does not exist
if [ ! -d "$result_dir" ]; then
    mkdir -p "$result_dir"
fi
# go to the base_path
cd "$base_path"
base_diff=$result_dir/base_diff.patch
commit_hash=$result_dir/base_commit_hash.txt
touch "$base_diff"
for next in $(git ls-files --others --exclude-standard); do
    git --no-pager diff --no-index /dev/null "$next" >>"$base_diff"
done
git rev-parse HEAD >"$commit_hash"

# cd into the analysis directory
analysis_path=$base_path/analysis_configurations/$analysis
if [ -d "$analysis_path" ]; then
    cd "$analysis_path"
    analysis_diff=$result_dir/analysis_diff.patch
    analysis_commit=$result_dir/analysis_commit_hash.txt
    analysis_name=$result_dir/analysis_name.txt
    touch "$analysis_diff"
    for next in $(git ls-files --others --exclude-standard); do
        git --no-pager diff --no-index /dev/null "$next" >>"$analysis_diff"
    done
    git rev-parse HEAD >"$analysis_commit"
    echo "$analysis" >"$analysis_name"
else
    echo "ERROR: analysis directory does not exist"
fi
