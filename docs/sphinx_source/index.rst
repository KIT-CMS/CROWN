Welcome to The CROWN documentation!
========================================

The The **C** ++-based **RO** OT **W** orkflow for **N** -tuples (CROWN) is a fast new way to convert NanoAOD samples into flat :code:`TTrees` to be used in further analysis. The main focus of the framework is to provide a fast and clean way of selecting events, calculating quantities and weights. The framework has minimal dependencies and only uses ROOT and it's Dataframe as a backend.

.. note::
   In order to get started, go here: :ref:`Getting started`.


.. image:: ../images/framework_workflow.svg
  :width: 900
  :align: center
  :alt: CROWN Workflow sketch

.. list-table:: Available Analyses Configurations for CROWN
   :widths: 25 150
   :header-rows: 1

   * - Analysis name
     - Repository
   * - ``tau``
     - https://github.com/KIT-CMS/TauAnalysis-CROWN
   * - ``earlyrun3``
     - https://github.com/khaosmos93/CROWN-config-earlyRun3



.. toctree::
   :maxdepth: 2

   introduction.rst

Documentation
==============

.. toctree::
   :maxdepth: 2
   :caption: How to contribute

   contrib.rst

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


Index
======

* :ref:`genindex`
* :ref:`search`
