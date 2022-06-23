# KingMaker



KingMaker is the workflow management for producing ntuples with the [CROWN](github.com/KIT-CMS/CROWN) framework. The workflow management is based on [law](github.com/riga/law), which is using [luigi](https://github.com/spotify/luigi) as backend.

---
## Setup

Setting up KingMaker should be straight forward:

```sh
git clone --recursive git@github.com:KIT-CMS/KingMaker.git
cd KingMaker
source setup.sh <Analysis Name>
```

this should setup the environment specified in the luigi.cfg file (located at `lawluigi_configs/<Analysis Name>_luigi.cfg`), which includes all needed packages.
The environment is sourced from the conda instance located at `/cvmfs/etp.kit.edu/LAW_envs/conda_envs/miniconda/` if possible. 
If the relevant environment is not available this way, the environment will be set up in a local conda instance.
The environment files are located at `conda_environments/<Analysis Name>_env.cfg`.
In addition other files are installed dependding on the analysis.

A list of available analyses can be found in the `setup.sh` skript or by running 
```sh
source setup.sh -l
```

In addition a `luigid` scheduler is also started if there isn't one running already. 

When setting up an already cloned version, a
```sh
source setup.sh <Analysis Name>
```
is enough.
---
# The KingMaker analysis
---
## Workflow

Currently, the workflow of the KingMaker analysis consists of four distinct tasks:

1. [ProduceSamples](processor/tasks/ProduceSamples.py)
    The main task, which is used to steer the Production of multiple samples at once
2. [CROWNRun](processor/tasks/CROWNRun.py)
    The task used to run CROWN with a specific file
3. [CROWNBuild](processor/tasks/CROWNBuild.py)
    This task is used to compile CROWN from source, and create a tarball, which is used by CROWNRun
4. [ConfigureDatasets](processor/tasks/ConfigureDatasets.py)
    This task is used to create NanoAOD filelists (if not existent) and readout the needed configuration parameters for each sample. This determines the CROWN tarball that is used for that job

---
## Run KingMaker

Normally, `KingMaker` can be run by running the `ProduceSamples` task. This is done using e.g.

```bash

law run ProduceSamples --local-scheduler False --analysis config --sample-list samples.txt --workers 1 --production-tag TestingCROWN

```
The required paramters for the task are:

1. `--local-scheduler False` - With this setting, the luigid scheduler is used  
2. `--analysis config` - The CROWN config to be used
3. `--sample-list samples.txt` - path to a txt file, which contains a list of nicks to be processed
4. `--production-tag TestingCROWN`

Additionally, some optional paramters are beneficial:

1. `--workers 1` - number of workers, currently, this number should not be larger than the number of tarballs to be built
2. `--print-status 2` - print the current status of the task
3. `--remove-output 2` - remove all output files
4. `--CROWNRun-workflow local` - run everything local instead of using HTCondor

---
## Tracking of Samples

The Samples are tracked and handled via nicks. For each sample, a unique nick has to be used. A collection of all samples and their settings is stored in the `datasets.yaml` file found in the `sample_database` folder. Additionally, a `nick.yaml` is generated for each individual sample, which contains all sample settings and a filelist of all `.root` files belonging to this sample.

```yaml
campaign: RunIISummer20UL18NanoAODv2
datasetname: DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8
datatier: NANOAODSIM
dbs: /DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM
energy: 13
era: 2018
extension: ''
filelist:
- root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18NanoAODv2/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/270000/D1972EE1-2627-2D4E-A809-32127A576CF2.root
- root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18NanoAODv2/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/50000/394775F4-CEDE-C34D-B56E-6C4839D7A027.root
- root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18NanoAODv2/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/50000/359B11BC-AE08-9B45-80A2-CC5EED138AB7.root
[....]
generators: amcatnloFXFX-pythia8
nevents: 85259315
nfiles: 85
nick: temp_nick
prepid: SMP-RunIISummer20UL18NanoAODv2-00030
sample_type: mc
status: VALID
version: 1
```

If a sample specific config is not available yet, `ConfigureDatasets` will perform a DAS query to get a filelist for this sample.

---
## Configuration

The default configuration provided in the repository should work out of the box. However some parameters might be changed. The configuration is spread across two files `lawluigi_configs/KingMaker_luigi.cfg` and `lawluigi_configs/KingMaker_law.cfg`. The HTCondor setting can also be found there.

---
## Other analyses
Analyses apart from KingMaker itself are still beeing worked on.
The `ML_LAW` analysis is an example for an analysis that aims to utilize remote GPU resources (e.g. machine learning applications).

