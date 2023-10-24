#!/bin/bash

FOUND_ISSUE=0

for FILENAME in $(find . -name "*.py");
do
    black --check "$FILENAME"
    RETURN_VALUE=$?
    if [ $RETURN_VALUE -ne 0 ]
    then
        black --diff "$FILENAME" 2> /dev/null
        FOUND_ISSUE=1
    fi
done

exit $FOUND_ISSUE
