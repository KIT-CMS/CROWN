#!/bin/bash

EXECUTABLE=$1

# Record data with the massif tool from valgrind.
# Take care, the exectuable runs considerable slower now!
valgrind --tool=massif --massif-out-file="massif.out" $EXECUTABLE

# Make the output understandable using ms_print
# Look for the last snapshot!
ms_print massif.out > massif.log
