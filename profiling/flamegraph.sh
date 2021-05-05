#!/bin/bash

EXECUTABLE=$1
INPUTFILE=$2
OUTPUTFILE=$3

# Record samples with perf
# perf record -g $EXECUTABLE
perf record --call-graph dwarf $EXECUTABLE $INPUTFILE $OUTPUTFILE
# Convert the data to be readable for flamegraphs
perf script > out.perf

# Perform the stack collapse
curl -Os https://raw.githubusercontent.com/eguiraud/FlameGraph/160b531f4c5ef0fec37e2b719ec609842a02aa99/stackcollapse-perf.pl > stackcollapse-perf.pl
perl stackcollapse-perf.pl out.perf > out.folded

# Generate the flamegraph
curl -Os https://raw.githubusercontent.com/brendangregg/FlameGraph/a258e78f17abdf2ce21c2515cfe8306a44774e2a/flamegraph.pl > flamegraph.pl
perl flamegraph.pl out.folded > flamegraph.svg
