#!/bin/bash

tarball=$1
checkout_location=$2
# check if tarball is provided
if [ -z "$tarball" ]; then
    echo "ERROR: path to a tarball or an unpacked folder must be provided"
    exit 1
fi
# create the checkout location if it does not exist
if [ ! -d "$checkout_location" ]; then
    mkdir -p "$checkout_location"
fi

if [ -d "$tarball" ]; then
    echo "Using tarball folder $tarball"
    mkdir temp
    cp -r "$tarball" temp/
    tarball=$(realpath "temp/$tarball")
# if tarball is a file, unpack it here
elif [ -f "$tarball" ]; then
    echo "Using tarball file $tarball"
    mkdir temp
    tar -xzf "$tarball" -C temp/
    tarball=$(realpath temp)
else
    echo "ERROR: tarball must be a file or a folder"
    exit 1
fi
# go to the checkout location
cd "$checkout_location"
echo "Using temporary tarball copy in $tarball"
# read the commit hashes from the diff folder in the tarball
base_commit=$(cat "$tarball/diff/base_commit_hash.txt")
analysis_commit=$(cat "$tarball/diff/analysis_commit_hash.txt")
analysis_name=$(cat "$tarball/diff/analysis_name.txt")

# checkout the base repository with the commit hash
git clone --recursive git@github.com:KIT-CMS/CROWN
cd CROWN
git checkout "$base_commit"
# apply the base diff
git apply "$tarball/diff/base_diff.patch"
# setup the analysis with the init script
bash init.sh "$analysis_name"
# checkout the analysis repository with the commit hash
cd "analysis_configurations/$analysis_name"
git checkout "$analysis_commit"
# apply the analysis diff
git apply "$tarball/diff/analysis_diff.patch"
# cleanup the temporary tarball
rm -rf "$tarball"
echo "**************************************************************"
echo "*  Setup from tarball finished. You can now run the analysis *"
echo "**************************************************************"
