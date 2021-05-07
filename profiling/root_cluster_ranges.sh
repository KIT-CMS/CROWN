#!/bin/bash

FILENAME=$1
TREENAME=$2

# This script uses ROOT's rootls tool to get the number of clusters in a file.
# The number of clusters is the basic limitation to multi-thread scaling.
# As a rule of thumb, you want about N x10 the number of clusters running on
# parallel on N threads.

rootls -t ${FILENAME}:${TREENAME} | grep "Cluster INCLUSIVE ranges:" -A 100000
