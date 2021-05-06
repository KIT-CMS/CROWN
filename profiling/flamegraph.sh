#!/bin/bash

EXECUTABLE=$1
INPUTFILE=$2
OUTPUTFILE=$3

# Record samples with perf
perf record --call-graph dwarf $EXECUTABLE $INPUTFILE $OUTPUTFILE

# Convert the data to be readable for flamegraphs
perf script > out.perf

# Perform the stack collapse
# The repo of Enrico (eguiraud) holds a small fix printing long templates in a nicer way.
BASE_URL=https://raw.githubusercontent.com/eguiraud/FlameGraph/dc043011dd74c9bdda9a9ec4254a17e45b76f560
curl -Os ${BASE_URL}/stackcollapse-perf.pl > stackcollapse-perf.pl
perl stackcollapse-perf.pl out.perf > out.folded

# Generate the flamegraph
curl -Os ${BASE_URL}/flamegraph.pl > flamegraph.pl
perl flamegraph.pl out.folded > flamegraph.svg
