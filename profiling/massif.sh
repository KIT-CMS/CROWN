#!/bin/bash

EXECUTABLE=$1
OUTPUTFILE=$2
INPUTFILE=$3

# Record data with the massif tool from valgrind.
# Take care, the exectuable runs considerable slower now!
valgrind --time-unit=ms --tool=massif --massif-out-file="massif.out" $EXECUTABLE $OUTPUTFILE $INPUTFILE

# Make the output understandable using ms_print
# Look for the last snapshot!
ms_print massif.out > massif.log
