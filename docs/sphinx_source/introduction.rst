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

   git clone git@github.com:KIT-CMS/CROWN.git

and source the current LCG stack (at the moment we use a nightly build)

.. code-block:: console

   source init.sh

after this, the framework should be installed

Running the framework
**********************

In order to create a new executable, first create a build directory

.. code-block:: console

   mkdir build && cd build

and then run `cmake` to setup the Makefiles

.. code-block:: console

   cmake ..

and compile the executable using

.. code-block:: console

   make install


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

It may require installing the following python packages beforehand

.. code-block:: console

   pip3 install --user breathe
   pip3 install --user sphinx_rtd_theme

The resulting documentation can than be found in

.. code-block:: console

   build_docs/docs/index.html

