#!/bin/bash

# Command to run ntuple generation of hello_world with open data example file
# Run this from the CROWN dir
ANALYSIS=hello_world
CONFIG=config
SAMPLES=dyjets
ERAS=2012
SCOPE=mm
EXAMPLE_FILE="root://eospublic.cern.ch//eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/Run2012BC_DoubleMuParked_Muons.root"

cd build
# Run cmake to generate c++ code
cmake .. -DANALYSIS=${ANALYSIS} -DCONFIG=${CONFIG} -DSAMPLES=${SAMPLES} -DERAS=${ERAS} -DSCOPES=${SCOPE}
# Compile generated code
nice -19 make install -j 20
# Use executable on Example file
cd bin
./${CONFIG}_${SAMPLES}_${ERAS} out.root ${EXAMPLE_FILE}
