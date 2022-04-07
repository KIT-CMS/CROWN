Introduction
=============

The The **C** ++-based **RO** OT **W** orkflow for **N** -tuples (CROWN) is a fast new way to convert NanoAOD samples into flat :code:`TTrees` to be used in further analysis. The main focus of the framework is to provide a fast and clean way of selecting events, calculating quantities and weights. The framework has minimal dependencies and only uses ROOT and it's Dataframe as a backend.


Design Idea
************

The framework consists of two main parts, a python configuration and a set of C++ functions. The python configuration is used to automatically generate a C++ script, which is then compiled to an executable using :code:`cmake` and all available compiler optimizations. This has the main advantage that the compiled executable is very fast and efficient in calculating the output `TTree`. In the following sketch, the overall workflow of CROWN is illustrated.

.. image:: ../images/framework_workflow.svg
  :width: 900
  :align: center
  :alt: CROWN Workflow sketch


Getting started
****************

Setting up the framework is very straight forward.

First clone the Repository

.. code-block:: console

   git clone --recurse-submodules git@github.com:KIT-CMS/CROWN.git

and source the current LCG stack (at the moment we use a nightly build)

.. code-block:: console

   source init.sh

after this, the framework should be installed

Running the framework
**********************

In order to create a new executable, first create a build directory

.. code-block:: console

   mkdir build && cd build

and then run `cmake` to setup the Makefiles. For the cmake command a minimal set of options has to be provided:

.. code-block:: console

   cmake .. -DANALYSIS=config -DSAMPLES=emb -DERAS=2018

The options that are currently available are:

   * :code:`-DANALYSIS=config`: The analysis configuration to be used. This is the name of the python configuration file. The file has to be located in the :code:`config` directory and the path is provided in the python import syntax so e.g. :code:`subfolder.myspecialconfig`
   * :code:`-DSAMPLES=emb`: The samples to be used. This is a single sample or a comma separated list of sample names.
   * :code:`-DERAS=2018`: The era to be used. This is a single era or a comma separated list of era names.
   * :code:`-DTHREADS=20`: The number of threads to be used. Defaults to single threading.
   * :code:`-DSHIFTS=all`: The shifts to be used. Defaults to all shifts. If set to :code:`all`, all shifts are used, if set to :code:`none`, no shifts are used, so only nominal is produced. If set to a comma separated list of shifts, only those shifts are used. If set to only a substring matching multiple shifts, all shifts matching that string will be produced e.g. :code:`-DSHIFTS=tauES` will produce all shifts containing :code:`tauES` in the name.
   * :code:`-DDEBUG=true`: If set to true, the code generation will run with debug information and the executable will be compiled with debug flags
   * :code:`-DOPTIMIZED=true`: If set to true, the compiler will run with :code:`-O3`, resulting in slower build times but faster runtimes. Should be used for developements, but not in production.
   * :code:`-DGENERATOR=Ninja`: The generator to be used. Defaults to Ninja. to set the generator to regular make files use :code:`-DGENERATOR="Unix Makefiles"`

and compile the executable using

.. code-block:: console

   ninja install

By default, the ninja_ build system is used for CROWN. However, the usage of other build systems is also possible and can be specified using the :code:`-G=` option, e.g. for regular makefiles use :code:`-DGENERATOR="Unix Makefiles"`, and then use the :code:`make install` command to compile the executable.
.. _ninja: https://ninja-build.org/

Creating Documentation
***********************

The Web documentation at readthedocs is updated automatically. However, if you want to create the documentation locally you have to first create a new build directory like :code:`build_docs`

.. code-block:: console

   mkdir build_docs && cd build_docs


then run :code:`cmake` to setup the documentation building process

.. code-block:: console

   cmake ../docs

and build the documentation using

.. code-block:: console

   make

The resulting documentation can than be found in

.. code-block:: console

   build_docs/docs/index.html

