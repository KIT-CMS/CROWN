# ToBeNamedPrototype

## Environment setup

Current requirement is a working mount to CERN's `/cvmfs/sft.cern.ch`

```bash
source init.sh
```

## Compilation for main executable

```bash
g++ analysis.cxx $(root-config --cflags --libs) -fconcepts
```
