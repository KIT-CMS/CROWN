#!/bin/bash

EXECUTABLE=$1
INPUTFILE=$2
OUTPUTFILE=$3

# Record data with the massif tool from valgrind.
# Take care, the exectuable runs considerable slower now!
valgrind --time-unit=ms --tool=massif --massif-out-file="massif.out" $EXECUTABLE $INPUTFILE $OUTPUTFILE

# Make the output understandable using ms_print
# Look for the last snapshot!
ms_print massif.out > massif.log
