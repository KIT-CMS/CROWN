#!/bin/bash

EXECUTABLE=$1
OUTPUTFILE=$2
INPUTFILE=$3
GRAPHNAME=$4

# Record samples with perf
# perf record -g $EXECUTABLE
perf record --call-graph dwarf $EXECUTABLE $OUTPUTFILE $INPUTFILE
# Convert the data to be readable for flamegraphs
perf script > out.perf
BASE_URL=https://raw.githubusercontent.com/eguiraud/FlameGraph/160b531f4c5ef0fec37e2b719ec609842a02aa99/
# Perform the stack collapse
curl -OLs ${BASE_URL}/stackcollapse-perf.pl > stackcollapse-perf.pl
perl stackcollapse-perf.pl out.perf > out.folded

# Generate the flamegraph
curl -OLs ${BASE_URL}/flamegraph.pl > flamegraph.pl
perl flamegraph.pl out.folded > flamegraph_${GRAPHNAME}.svg

#cleanup
rm out.perf out.folded perf.data perf.data.old
