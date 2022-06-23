#!/bin/bash

FOUND_ISSUE=0

for FILENAME in $(find . -path ./law -prune -o -name "*.py");
do
    # only run the check if the filename ends with .py
    if [[ $FILENAME == *.py ]]; then
        echo "Checking $FILENAME"
        black --check $FILENAME
        RETURN_VALUE=$?
        if [ $RETURN_VALUE -ne 0 ]
        then
            black --diff $FILENAME 2> /dev/null
            FOUND_ISSUE=1
        fi
    fi
done

exit $FOUND_ISSUE
