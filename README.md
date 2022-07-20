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

# The ML_train workflow
---
## Workflow

The ML_train workflow currently contains a number of the tasks necessary for the `htt-ml NMSSM` analysis. Non-`NMSSM` analyses are currently not yet supported. The workflow uses the `PuppetMaster` task to speed up status queries. It should also be noted, that all created files are stored in remote storage and might be subject to file caching under certain circumstances.
The tasks are:

1. [CreateTrainingDataShardConfig](processor/tasks/CreateTrainingDatasets.py)
    Task that creates configuration files from which the process-datasets for the machine learning tasks are created. Uses the [write_datashard_config](sm-htt-analysis/ml_datasets/write_datashard_config.py) script. \
    Some aspects of the process-datasets depend on the decay channel. For this reason the datasets for different channels are handled in seperate tasks. Although not specifically necessary, the different run eras are also handled in seperate tasks.
2. [CreateTrainingDataShard](processor/tasks/CreateTrainingDatasets.py)
    Remote workflow task that creates the process-datasets for the machine learning tasks from the config files created by the `CreateTrainingDataShardConfig` Task. The task uses the `ntuples` and `friend trees` described in the [Sample setup](sm-htt-analysis/utils/setup_samples.sh). These dependencies are currently not checked by LAW. Also uses the [create_training_datashard](sm-htt-analysis/ml_datasets/create_training_datashard.py) script. \
    The task branches each return a root file that consists of only one fold of one process of one run era and one decay channel. These files can then be used for the machine learning tasks. Some aspects of the process-datasets depend on the decay channel. For this reason the datasets for different channels are handled in seperate tasks. Although not specifically necessary, the different run eras are also handled in seperate tasks.
3. [CreateTrainingConfig](processor/tasks/RunTraining.py)
    Task that creates configuration files used for the machine learning tasks. Uses the [write_datashard_config](sm-htt-analysis/ml_trainings/write_training_config.py) and the [create_combined_config](sm-htt-analysis/ml_trainings/create_combined_config.py) scripts. \
    The resulting config files contain information about which processes are used in the training, the combined event weights of the classes, and multiple hyperparameters that influence the training process. The training task is able to merge the datasets of the different run eras for a single training. For this reason the config file for such a training is a combination of the config files of the different run eras. At this point the seperate run era tasks are merged and only the different decay channels are still treated as seperate tasks.
4. [RunTraining](processor/tasks/RunTraining.py)
    Remote workflow task that performs the neural network training using GPU resources if possible. Uses the [keras_training](sm-htt-analysis/ml_trainings/keras_training.py) script. The hyperparameters of this training are provided by the config files of the `CreateTrainingConfig` task. The [get_processes](sm-htt-analysis/utils/get_processes.py) script is used to select the signal and background processes for the individual training groups. The script can currntly only handle the `NMSSM` groups. \
    Each branch task returns a set of files for one fold of one training group of one decay channel. Each set includes the trained `.h5` model, the preprocessing object as a `.pickle` file and a graph of the loss as a `.pdf` and `.png`. This task can run all trainings in parallel without splitting them into seperate tasks, so the tasks of the different decay channels are merged here.
5. [RunAllTrainings](processor/tasks/RunTraining.py)
    Task to run all possible trainings.

## Run ML_train

Normally, the `ML_train` workflow can be run by running the `RunTraining` task. This is done using e.g.

```bash

law run RunTraining --workers 2

```
There are a number of parameters to be set in the [luigi](lawluigi_configs/ML_train_luigi.cfg) and [law](lawluigi_configs/ML_train_law.cfg) config files:

1. The run era. Can be either `2016`, `2017`, `2018` or `all_eras`.
2. The decay channels. A list of channels that can include `tt`, `et` and `mt`.
3. The masses of the heavy additional higgs boson. A list that can include `240`, `280`, `320`, `360`, `400`, `450`, `500`, `550`, `600`, `700`, `800`, `900`, `1000` and `heavier`.
4. The groups of masses of the light additional higgs boson. A list that can include `1`, `2`, `3`, `4`, `5`, `6` and `7`.\
**Note**: Only valid combinations of heavy and light masses are used, as determined by the `valid_batches` function of `RunTraining`.\
By default these parameters are already set to show the proper syntax.\
If the `RunAllTrainings` task is used, all of the parameters are set, ignoring the config file.
5. Optional: The production_tag. Can be any string. Used to differentiate the runs. Default is a unique timestamp.

Command line arguments:
1. `--workers`; The number of tasks that are handled simultaneously. Due to the usage of the `PuppetMaster` Task in this workflow, this parameter should be set to 2 times the maximum number of tasks that will be run in parallel. For `RunAllTrainings` this means `3(eras) x 3(channels) x 2 = 18`. Default is 1.
2. `--print-status -1`; Return the current status of all tasks involved in the workflow. Shows only the status of the `PuppetMaster` Tasks for all requirements, which is a summary of the puppeteered task statuses.
3. `--remove-output -1`; Remove all `PuppetMaster` output files as well as the output files of the explicitly called task (like `RunTraining`). Output files of puppeteered tasks are currently not removed and have to be cleaned up by hand.

## Not yet implemented
A number of necessary scripts are not yet implemented. This includes the ML testing scripts as well as the conversion of the `.h5` network files to a different format.

---
# Other analyses
Analyses apart from KingMaker itself are still beeing worked on.

---
# Additional Features
A collection of additional function and tasks have added to the central [framework](processor/framework.py).

## PuppetMaster
A task to speed up the gathering of task statuses when many files are involved. More relevant for remote file storage.
Acts as the tasks it is puppeteering during the status query. Dynamically adds its task to the shedduler at runtime. Checks if given task is still the same as during previous (succesfull) executions by comparing their output targets.
Only prints full representation of puppet task during status queries, if the `fulltask` parameter is set.
Used by giving the required task to the puppet before returning it in a tasks requires function. If multiple tasks of the same kind are used in a workflow the `identifier` parameter has to be used to distinguish them.
Example:
```python
return PuppetMaster(puppet_task=Task(**requirements), identifier=[channel])
```
For workflow tasks only the `workflow_requires` should require a `PuppetMaster` tasks as the individual branches require the actual task data of the puppet. In addition the `workflow_requires` function should be skipped in remote branch tasks as the `PuppetMaster` output files are only stored locally. This can be done by adding 
```python
if self.is_branch():
    return None
```
to the top of the `workflow_requires` function.
Some syntax has to be altered to recieve the correct input data in non-workflow tasks.
Instead of using `self.input()["Name"]`, `self.requires()["Name"].give_puppet_outputs()` should be used.
An example how all of this should be used can be found in the `ML_train` workflow.
