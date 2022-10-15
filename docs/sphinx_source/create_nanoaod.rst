How to generate dummy NanoAOD samples
==========================================

We use self-generated NanoAOD samples for testing, which can be produced with the open source CMSSW framework (see https://github.com/cms-sw/cmssw). The following script generates a file from scratch using fast simulation (see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFastSimulationExamples).

The used files are hosted in the repository https://github.com/KIT-CMS/CROWNTestingSamples.

It should be noted that these files must be used solely for testing purposes and do not contain any useful information for physics analysis.

.. code-block:: console

    #!/bin/bash

    set -e

    CMSSW_VERSION=CMSSW_10_6_24
    THIS=$(pwd)

    # Set up CMSSW
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    if [[ -d $CMSSW_VERSION ]]
    then
        # Just source CMSSW
        cd $CMSSW_VERSION/src
        eval `scramv1 runtime -sh`
        cd $THIS
    else
        # Check out CMSSW
        scram project CMSSW $CMSSW_VERSION
        cd $CMSSW_VERSION
        eval `scramv1 runtime -sh`

        # Add generator configurations and build
        git cms-addpkg Configuration/Generator
        scram b

        cd $THIS
    fi

    # Create an output directory
    OUTPUT_DIR=output
    mkdir -p $OUTPUT_DIR
    cd $OUTPUT_DIR

    # Settings
    MINBIAS_CONFIG=MinBias_13TeV_pythia8_TuneCUETP8M1_cfi
    PROCESS_CONFIG=GGToHtautau_13TeV_pythia8_cff
    NUM_EVENTS=100
    CONDITIONS=auto:run2_mc
    ERA=Run2_2016

    # Generate the AOD file (with pileup mixing)
    cmsDriver.py ${MINBIAS_CONFIG} \
        --mc \
        --fast \
        --step GEN,SIM,RECOBEFMIX \
        --conditions $CONDITIONS \
        --era $ERA \
        -n ${NUM_EVENTS} \
        --eventcontent FASTPU \
        --datatier GEN-SIM-RECO

    cmsDriver.py ${PROCESS_CONFIG}  \
        --mc \
        --fast  \
        --step GEN,SIM,RECOBEFMIX,DIGI:pdigi_valid,L1,DIGI2RAW,L1Reco,RECO \
        --conditions $CONDITIONS \
        --era $ERA \
        -n ${NUM_EVENTS} \
        --eventcontent AODSIM \
        --datatier AODSIM \
        --pileup AVE_35_BX_25ns \
        --pileup_input file://${MINBIAS_CONFIG}_GEN_SIM_RECOBEFMIX.root

    # Generate the miniAOD
    cmsDriver.py miniAOD \
        --fast \
        --mc \
        --step PAT \
        --conditions $CONDITIONS \
        --era $ERA \
        --eventcontent MINIAODSIM \
        -n ${NUM_EVENTS} \
        --runUnscheduled \
        --filein file://${PROCESS_CONFIG}_GEN_SIM_RECOBEFMIX_DIGI_L1_DIGI2RAW_L1Reco_RECO_PU.root

    # Generate the nanoAOD
    cmsDriver.py nanoAOD \
        --fast \
        --mc \
        --step NANO \
        --conditions $CONDITIONS \
        --era $ERA \
        --eventcontent NANOAODSIM \
        -n ${NUM_EVENTS} \
        --datatier NANOAODSIM \
        --filein file://miniAOD_PAT.root \
        --fileout file://nanoAOD.root
