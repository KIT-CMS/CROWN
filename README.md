# ToBeNamedPrototype

## Environment setup

Current requirement is a working mount to CERN's `/cvmfs/sft.cern.ch`

```bash
source init.sh
```

## Compilation with cmake

```bash
mkdir build/ && cd build/
cmake .. # or cmake3 on centos7
make -j8 # or 'VERBOSE=1 make' for the verbose output
cd ..
```
Executable will then be `build/Analysis` instead of `a.out`

## Compilation for main executable

```bash
g++ analysis.cxx $(root-config --cflags --libs) -fconcepts
```

## Profiling with perf & flamegraph for CPU

See the script [`profiling/flamegraph.sh`](profiling/flamegraph.sh).

Running profiling on executable

```bash
perf record ./a.out -g
```

Print out the report

```bash
perf report
```

Get flamegraph repo

```bash
git clone https://github.com/brendangregg/FlameGraph
```

Create flamegraph

```bash
perf script > out.perf
FlameGraph/stackcollapse-perf.pl out.perf > out.folded
FlameGraph/flamegraph.pl out.folded > out.svg
```

## Profiling with valgrind massif for Memory

See the script [`profiling/massif.sh`](profiling/massif.sh).

```bash
valgrind --tool=massif ./a.out
ms_print massif.out.4103388 > massif.log
```
