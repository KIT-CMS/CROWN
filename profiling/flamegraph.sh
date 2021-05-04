#!/bin/bash

EXECUTABLE=$1
INPUTFILE=$2
OUTPUTFILE=$3

# Record samples with perf
perf record -g $EXECUTABLE $INPUTFILE $OUTPUTFILE

# Convert the data to be readable for flamegraphs
perf script > out.perf

# Perform the stack collapse
curl -Os https://raw.githubusercontent.com/brendangregg/FlameGraph/master/stackcollapse-perf.pl > stackcollapse-perf.pl
perl stackcollapse-perf.pl out.perf > out.folded

# Generate the flamegraph
curl -Os https://raw.githubusercontent.com/brendangregg/FlameGraph/master/flamegraph.pl > flamegraph.pl
perl flamegraph.pl out.folded > flamegraph.svg
