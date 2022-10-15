#!/bin/bash

for FILENAME in $(find . -name "*.[h,c]xx");
do
    clang-format -i $FILENAME
done

DIFF=$(git --no-pager diff)
if [[ -z $DIFF ]];
then
    echo "Nothing to format, all good!"
    exit 0
else
    echo "Found following diff with clang-format:"
    git --no-pager diff
    echo "Please format the code with clang-format -i /path/to/files using the .clang-format config from this repository"
    exit 1
fi
