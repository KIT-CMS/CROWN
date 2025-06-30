Welcome to The CROWN documentation!
####################################

The **C** ++-based **RO** OT **W** orkflow for **N** -tuples (CROWN) is a fast new way to convert NanoAOD samples into flat :code:`TTrees` to be used in further analysis. The main focus of the framework is to provide a fast and clean way of selecting events and calculating quantities and weights. The framework has minimal dependencies and only uses ROOT and it's Dataframe as a backend.

.. note::
   To get started with CROWN, go here: :ref:`Getting started`.

.. note::
  To get started with an ntuple production workflow that uses CROWN, go here: :ref:`Workflow Management`.

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
   * - BSM di-Higgs search X &rightarrow; YH &rightarrow; bb&tau;&tau;
     - https://github.com/KIT-CMS/XYHBBTauTauAnalysis-CROWN
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
   kingmaker.rst
   postprocessing.rst
   friend_trees.rst
   cms_analysis.rst
   changelog.rst

.. toctree::
   :maxdepth: 2
   :caption: Setup your own Configuration

   py_configuration.rst
   contrib.rst
   cpp_addons.rst
   correction_manager.rst
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
