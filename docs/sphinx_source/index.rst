Welcome to The CROWN documentation!
####################################

The The **C** ++-based **RO** OT **W** orkflow for **N** -tuples (CROWN) is a fast new way to convert NanoAOD samples into flat :code:`TTrees` to be used in further analysis. The main focus of the framework is to provide a fast and clean way of selecting events, calculating quantities and weights. The framework has minimal dependencies and only uses ROOT and it's Dataframe as a backend.

.. note::
   In order to get started, go here: :ref:`Getting started`.

.. note::
   To read about recent changes and new features, go here: :ref:`changelog`.


Available Analyses
*******************
   The following analyses configurations are currently available in CROWN. If you want to add your own analysis configuration, contact the developers.

.. list-table:: Available Analyses Configurations for CROWN
   :widths: 25 150
   :header-rows: 1

   * - Analysis name
     - Repository
   * - ``HTauTau``
     - https://github.com/KIT-CMS/TauAnalysis-CROWN
   * - ``earlyrun3``
     - https://github.com/khaosmos93/CROWN-config-earlyRun3
   * - ``WHTauTau``
     - https://github.com/KIT-CMS/WHTauTauAnalysis-CROWN


Documentation Content
######################

.. toctree::
   :maxdepth: 2

   introduction.rst
   changelog.rst
   friend_trees.rst
   kingmaker.rst
   postprocessing.rst

Documentation
**************

.. toctree::
   :maxdepth: 2
   :caption: Python Configuration

   py_configuration.rst
   py_classes.rst

.. toctree::
   :maxdepth: 2
   :caption: C++ Configuration

   c_functions.rst
   namespaces.rst


.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   build_root.rst
   create_nanoaod.rst

.. toctree::
   :maxdepth: 2
   :caption: How to contribute

   contrib.rst


Index
******

* :ref:`genindex`
* :ref:`search`
