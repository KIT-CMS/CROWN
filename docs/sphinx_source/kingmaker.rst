KingMaker
===========

KingMaker is a workflow management for producing ntuples with the CROWN framework. The workflow management is based on law (https://github.com/riga/law), which uses luigi (https://github.com/spotify/luigi) as the backend.

Setup
-----

.. code-block:: bash

    git clone --recursive git@github.com:KIT-CMS/KingMaker.git
    cd KingMaker
    source setup.sh KingMaker

This should install all required packages and set up the environment. In addition, a ``luigid`` scheduler is started, if not already running. The required port is set to the ```LUIGIPORT``` environment variable.

Management of samples
---------------------

Samples can be managed manually or using the ``sample_manager``, which can be started with

.. code-block:: bash

    sample_manager

This starts a CLI, which can be used to add more samples to the database, update samples or quickly generate a sample list for producing ntuples.

Addition of new Samples
~~~~~~~~~~~~~~~~~~~~~~~

When adding a new sample, follow the instructions of the ``sample_manager``. In the background, the DAS database of CMS is queried, to get samples, matching the provided dataset name:

.. code-block::

    Starting up Datasetmanager
    A working version of the database exists
    ? Load working version of database ? No
    Database loaded
    The database contains 581 samples, split over 4 era(s) and 22 sampletype(s)
    ? What do you want to do? Add a new sample
    ? Enter a DAS nick to add /DYJetsToLL_M-50_*/RunIISummer20UL16NanoAOD*v9-106X*/NANOAODSIM
    Multiple results found
    ? Which dataset do you want to add ? (Use arrow keys to move, <space> to select, <a> to toggle, <i> to invert)
    » ○ Nick: /DYJetsToLL_M-50_TuneCH3_13TeV-madgraphMLM-herwig7/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM - last changed: 03 Jan 2023 11:05 - created: 30 Nov 2022 14:26
    ○ Nick: /DYJetsToLL_M-50_TuneCH3_13TeV-madgraphMLM-herwig7/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM - last changed: 08 Nov 2022 13:17 - created: 08 Nov 2022 05:15
    ○ Nick: /DYJetsToLL_M-500to700_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM - last changed: 05 Nov 2022 22:12 - created: 04 Nov 2022 00:52
    ○ Nick: /DYJetsToLL_M-500to700_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM - last changed: 15 Sep 2022 00:09 - created: 14 Sep 2022 22:32
    ○ Nick: /DYJetsToLL_M-50_Zpt-100to200_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM - last changed: 05 May 2022 07:44 - created: 26 Apr 2022 06:
    ○ Nick: /DYJetsToLL_M-50_Zpt-100to200_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM - last changed: 28 Mar 2022 21:51 - created: 28 Mar 2022 19:42
    ○ Nick: /DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM - last changed: 20 Feb 2022 06:54 - created: 17 Feb 2022 22:29
    ○ Nick: /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM - last changed: 29 Nov 2021 11:10 - created: 28 Nov 2021 07:54

The results will be sorted by time, starting the the newest samples on top. The query name has to match the CMS conventions ``/*/*/*``. Select all samples that you want to add. Afterwards, set the correct sample type. Optionally, the ``sample_manager`` can calculate the GeneratorWeight associated to the sample. Since this process can take some time, the task can also be triggered afterwards.


Generation of sample lists
~~~~~~~~~~~~~~~~~~~~~~~~~~

To generate a sample list select ``Create a production file``

.. code-block::

    The database contains 581 samples, split over 4 era(s) and 22 sampletype(s)
    ? What do you want to do? (Use arrow keys)
        ○ Add a new sample
        ○ Edit a sample (not implemented yet)
        ○ Delete a sample
        ○ Find samples (by nick)
        ○ Find samples (by DAS name)
        ○ Print details of a sample
        » ○ Create a production file
        ○ Update genweight
        ○ Save and Exit
        ○ Exit without Save
        Answer: Create a production file


In the next step, select the eras you want to use using the arrow keys and space bar

.. code-block::

    Select eras to be added  (Use arrow keys to move, <space> to select, <a> to toggle, <i> to invert)
    ○ 2016postVFP
    ○ 2016preVFP
    ● 2018
    » ● 2017


and then select the sample types you want to process. The output file will be a ``.txt`` file, which can be used to produce ntuples.


Submission of ntuples
---------------------

In Kingmaker, two tasks are separated:

1. The production of ntuples
2. The production of friend trees

The first task is handled by the ``ProduceFriends`` task, the latter by the ``ProduceFriends`` task. In the case of friend trees, missing Ntuples are generated automatically.

.. warning::
    By default, KingMaker will write all outputs to the GridKA NRG storage. As a result, the user has to provide a valid X509 proxy, and the environment variable ``X509_USER_PROXY`` has to be set. The proxy can be created using ``voms-proxy-init``. The proxy has to be valid for at least 24 hours. The proxy can be checked using ``voms-proxy-info``.
    To use a different output storage, the KingMaker configuration has to be adapted, more details can be found in the :ref:`KingMaker Configuration` section

Production of NTuples
~~~~~~~~~~~~~~~~~~~~~

To trigger a production of ntuples run

.. code-block:: bash

    law run ProduceSamples --analysis tau --config config --sample-list samples.txt --production-tag debug_2 --workers 10 --scopes mt --print-status -1


The different options are:

- ``--analysis``: The analysis to be used. The name corresponds to the analysis folder in the ``CROWN/analysis_configurations`` folder.
- ``--config``: The config file to be used. The config file contains the information about the samples, the input files, the output files, the friend trees, the branches to be read, etc. The config file is located in the ``CROWN/analysis_configurations/<analysis>/config`` folder.
- ``--sample-list``: The sample list to be used. The sample list can be generated by the ``sample_manager`` and contains the information about the samples to be processed. The sample nicks can also be provided as a comma-separated list.
- ``--production-tag``: The production tag is used to identify the production. It is used to create the output folder and the output files. The output files are stored in the ``/<base>/<production-tag>/CROWNRun/`` folder. The ``base`` variable is set using the Configuration. By default, it is set to ``root://cmsxrootd-kit.gridka.de//store/user/${USER}/CROWN/ntuples/``. Within the folder, the different samples are stored, matching the ``<era>/<samplenick>/<channel>/<samplenick>_<counter>.root`` pattern.
- ``--workers``: The number of workers to be used. Each worker is responsible for the submission and handling of one sample. The number of workers should be at least the number of samples.
- ``--scopes``: The scopes to be used, provided as a comma-separated list.
- ``--shifts``: The shifts to be used, provided as a comma-separated list. If no shifts are provided, no shifts are applied. If ``All`` is provided, all shifts are applied, if ``None`` is provided, no shifts are applied.

.. warning::
    The law processes can get stuck after building the tarball when trying to upload it to the dCache when using more than 1 worker. The task will be stuck indefinitely. To avoid this, the user must cancel the running law command using ``Ctrl+C``. Afterwards, the task can be restarted using the same command. The task will then continue with the upload of the tarball. The reason for this behaviour is unknown.


Additionally, the following options can be useful:

- ``--print-status -1``: Print the status of the tasks. If ``-1`` is provided, the status of every task is printed.
- ``--remove-output -1``: Remove the output files. This option is useful if the production fails and the output files should be removed. This will trigger an interactive CLI, where only parts of the production can be removed as well.
- ``--CROWNRun-workflow local``: This option can be used to run the production locally. This is useful for debugging purposes if the batch system is currently not available. However, be aware, that this option should only run with a limited amount of workers and samples since it is very easy to overload the local machine.

.. warning::
    When using the dCache as Ntuple storage, the remove option should be used with care. Since the dCache caches files without checking, if the file content changes, overwriting files can lead to errors, where the old file is still cached. The saver option is to remove the old files and store the new files using a separate ``production-tag``.


Production of friend trees
~~~~~~~~~~~~~~~~~~~~~~~~~~

For the production of friend trees, the same options as for the production of ntuples are available. An example command is given below:

.. code-block:: bash

    law run ProduceFriends --analysis tau --config config --friend-config tau_friends --sample-list samples.txt --shifts None --friend-name test --production-tag debugging_v81 --workers 2

Some additional options are required:

- ``--friend-config``: The friend config file to be used. The friend config file contains the information about the friend trees to be produced. The friend config file is located in the ``CROWN/analysis_configurations/<analysis>/config`` folder.
- ``--friend-name``: The name of the friend tree to be produced. The name has to match the name in the friend config file. The resulting friend trees will be stored in the ``/<base>/<production-tag>/CROWNFriends/<friend-name>/`` folder.

The resulting folder structure for the command listed above will be

.. code-block::

    /<base>/<production-tag>/
        |- CROWNRun/
                        |- <era>/<samplenick>/<channel>/<samplenick>_<counter>.root
        |- CROWNFriends/
                        |- test/<era>/<samplenick>/<channel>/<samplenick>_<counter>.root

Production of friend trees with additional friends as input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the requested friend tree depends on additional friend trees the ``ProduceMultiFriends`` workflow is used for the final friend tree. The command

.. code-block:: bash

    law run ProduceFriends --analysis tau --config config --friend-config tau_classifier --friend-dependencies tau_friends,tau_friends_2 --sample-list samples.txt --shifts None --friend-name special_tau_classifier --production-tag debugging_v81 --workers 2

contains an additional option

- ``--friend-dependencies``: A list of additional configurations to be run. The list has to be provided as a comma-separated list. The resulting friend trees will be stored in the ``/<base>/<production-tag>/CROWNFriends/<friend-config>/`` folder. To set the name for the intermediate friend trees, two options are available. By default, the name of the configuration will be used as the name of the friend tree. Alternatively, the parameter ``friend_mapping`` can be used, to define a dictionary, where a mapping between the friend config name and the friend tree name can be defined. The dictionary has to be provided as a JSON string. An example is given below:

.. code-block:: python

    friend_mapping = {"tau_friends": "tau_friends_leptonscalefactors", "tau_friends_2": "tau_friends_svfit"}

As an example, the command listed above will produce not only ntuples for all samples specified in ``samples.txt`` using the config but also the friend trees ``tau_friends`` and ``tau_friends_2``. All those three inputs will then be used, to produce the final friend tree ``special_tau_classifier``. The resulting folder structure will be

.. code-block::

    /<base>/<production-tag>/
        |- CROWNRun/
                        |- <era>/<samplenick>/<channel>/<samplenick>_<counter>.root
        |- CROWNFriends/
                        |- tau_friends/<era>/<samplenick>/<channel>/<samplenick>_<counter>.root      (name automatically generated)
                        |- tau_friends_2/<era>/<samplenick>/<channel>/<samplenick>_<counter>.root    (name automatically generated)
        |- CROWNMultiFriends/
                        |- tau_classifier/<era>/<samplenick>/<channel>/<samplenick>_<counter>.root

if no ``friend_mapping`` is used, or

.. code-block::

    /<base>/<production-tag>/
        |- CROWNRun/
                        |- <era>/<samplenick>/<channel>/<samplenick>_<counter>.root
        |- CROWNFriends/
                        |- tau_friends_leptonscalefactors/<era>/<samplenick>/<channel>/<samplenick>_<counter>.root      (name automatically generated)
                        |- tau_friends_svfit/<era>/<samplenick>/<channel>/<samplenick>_<counter>.root    (name automatically generated)
        |- CROWNMultiFriends/
                        |- tau_classifier/<era>/<samplenick>/<channel>/<samplenick>_<counter>.root

with the exmaple ``friend_mapping`` mentioned above.

To perform the generation of friend trees locally, use

- ``--CROWNFriends-workflow local --CROWNRun-workflow local``: This option can be used to run the production locally. This is useful for debugging purposes if the batch system is currently not available. However, be aware, that this option should only run with a limited amount of workers and samples since it is very easy to overload the local machine.



KingMaker Configuration
-----------------------

The two relevant configuration files can be found in the ``lawluigi_configs`` folder. They are called ``KingMaker_law.cfg`` and ``KingMaker_luigi.cfg``.

.. warning::
    Most default parameters in the Configuration are chosen such that only minimal changes are required. Nevertheless, the user should check the configuration files before running KingMaker.

In the ``KingMaker_law.cfg`` file, the different tasks are defined. Also, the remote filesystem is defined here:

.. code-block::

    [wlcg_fs]
    base: root://cmsxrootd-kit.gridka.de//store/user/${USER}/CROWN/ntuples/
    use_cache: True
    cache_root: /tmp/${USER}/
    cache_max_size: 20000

In general, it is good practice to use the ``use_cache`` option. This will cache the files locally, which can speed up the processing. The ``cache_max_size`` option defines the maximum size of the cache in MB. If the cache is full, the oldest files are removed from the cache.

The ``base`` option defines the base path of the remote filesystem. The ``${USER}`` variable is replaced by the username. The ``base`` path is used to define the output path of the ntuples and friend trees. The ``base`` path is also used to define the input path of the friend trees. The ``base`` path should be set to the path of the dCache storage.

The ``KingMaker_luigi.cfg`` file contains the configuration of the different tasks. The most important options are defined in the ``[DEFAULT]`` section and include setting for the HTCondor job submission. Parameters defined in the ``[DEFAULT]`` section can be overwritten in the task-specific sections.

.. code-block::

    name = KingMaker
    ENV_NAME = KingMaker
    wlcg_path = root://cmsxrootd-kit.gridka.de//store/user/${USER}/CROWN/ntuples/
    htcondor_accounting_group = cms.higgs
    htcondor_remote_job = True
    htcondor_universe = docker
    htcondor_docker_image = mschnepf/slc7-condocker:latest
    transfer_logs = True
    local_scheduler = True
    tolerance = 0.00
    acceptance = 1.00
    ; submit only missing htcondor workflow branches (should always be true)
    only_missing = True

    ; bootstrap file to be sourced at beginning of htcondor jobs (relative PATH to framework.py)
    bootstrap_file = setup_law_remote.sh
    files_per_task = 10
    ; scopes and shifts are to be provided in the config, or as command line arguments via --scope and --shift
    ; in both cases, the values are expected to be comma-separated lists without spaces or quotes
    scopes = mt,et
    shifts = None

Here, the ``wlcg_path`` option should be set to the same path, as the ``base`` path in the ``KingMaker_law.cfg``. The different ``htconddor_`` parameters have to be adopted according to the requirements of the batch system. For the two tasks, that are run remotely, different job requirements can be set. The ``files_per_task`` option defines the number of files to be processed per task. The ``scopes`` and ``shifts`` options define the scopes and shifts to be used. These two parameters can also be provided as command line arguments, which is the recommended way.

.. code-block::

    [CROWNRun]
    ; HTCondor
    htcondor_walltime = 10800
    htcondor_request_memory = 16000
    htcondor_requirements = TARGET.ProvidesCPU && TARGET.ProvidesIO
    htcondor_request_disk = 20000000
    htcondor_request_cpus = 4
    # for these eras, only one file per task is processed
    problematic_eras = ["2018B", "2017C", "2016B-ver2"]

    [CROWNFriends]
    ; HTCondor
    htcondor_walltime = 10800
    htcondor_request_memory = 16000
    htcondor_requirements = TARGET.ProvidesCPU && TARGET.ProvidesIO
    htcondor_request_disk = 20000000
    # friends have to be run in single core mode to ensure a correct order of the tree entries
    htcondor_request_cpus = 1

The ``problematic_eras`` option is used to define eras, where only one file per task is processed. This can be required, if the NanoAOD input files have a change in their structure, e.g. if trigger paths are modified. To avoid problems in these cases, jobs can be processed with only one input file. This will slow down the processing but ensures that the processing is not stopped by a single file. Disk, wall time and other requirements can be set in the task-specific sections.

.. warning::
    For friend trees, multiprocessing is not possible, since the resulting friend tree must have the same order as the input tree. Therefore, the ``htcondor_request_cpus`` option has to be set to 1, which will disable multiprocessing.

For a more complete description of the different options, please refer to the overcomplete configuration in the law repository (https://github.com/riga/law/blob/master/law.cfg.example).