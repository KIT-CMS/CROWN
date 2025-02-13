Welcome to The CROWN documentation!
####################################

The **C** ++-based **RO** OT **W** orkflow for **N** -tuples (CROWN) is a fast new way to convert NanoAOD samples into flat :code:`TTrees` to be used in further analysis. The main focus of the framework is to provide a fast and clean way of selecting events and calculating quantities and weights. The framework has minimal dependencies and only uses ROOT and it's Dataframe as a backend.

.. note::
   To get started, go here: :ref:`Getting started`.

.. note::
   To read about recent changes and new features, go here: :ref:`changelog`.


Available Analyses
*******************
The following analysis configurations are currently available in CROWN. If you want to add your analysis configuration, contact the developers.

.. list-table:: Available Analyses Configurations for CROWN
   :widths: 25 150
   :header-rows: 1

   * - Analysis name
     - Repository
   * - H to TauTau
     - https://github.com/KIT-CMS/TauAnalysis-CROWN
   * - BSM di-Higgs to TauTauBB
     - https://github.com/KIT-CMS/TauAnalysis-CROWN/tree/nmssm_devs
   * - W/Z early Run3
     - https://github.com/KIT-CMS/earlyRun3Analysis-CROWN
   * - W + H to TauTau
     - https://github.com/KIT-CMS/WHTauTauAnalysis-CROWN
   * - Boosted H to TauTau
     - https://github.com/KIT-CMS/BoostedHiggsTauTauAnalysis-CROWN
   * - Single top
     - https://github.com/nfaltermann/CROWNs


Documentation Content
######################

.. toctree::
   :maxdepth: 2

   introduction.rst
   changelog.rst
   kingmaker.rst
   postprocessing.rst
   friend_trees.rst

.. toctree::
   :maxdepth: 2
   :caption: Setup your own Configuration

   contrib.rst
   py_configuration.rst
   correction_manager.rst
   cpp_addons.rst
   nanoAODversions.rst

.. toctree::
   :maxdepth: 2
   :caption: Documentation

   py_classes.rst
   namespaces.rst

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   build_root.rst
   create_nanoaod.rst


Index
******

* :ref:`genindex`
* :ref:`search`
