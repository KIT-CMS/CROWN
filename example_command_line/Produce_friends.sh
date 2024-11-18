#!/bin/bash

# Command to run friend generation of solution with open data example file
# Run this from the CROWN dir
ANALYSIS=solution
CONFIG=config_friends
SAMPLES=data
ERAS=2012
SCOPE=mm
QUANTITIESMAP=${PWD}/build/bin/out_${SCOPE}.root

source init.sh
mkdir build -p
cd build
# Run cmake to generate c++ code
cmake .. -DANALYSIS=${ANALYSIS} -DCONFIG=${CONFIG} -DSAMPLES=${SAMPLES} -DERAS=${ERAS} -DSCOPES=${SCOPE} -DQUANTITIESMAP=${QUANTITIESMAP}
# Compile generated code
nice -19 make install -j 20
# Use executable on Example file
cd bin
./${CONFIG}_${SAMPLES}_${ERAS}_${SCOPE} out_friend.root ${QUANTITIESMAP}
