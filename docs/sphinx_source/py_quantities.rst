Quantities in the python part
=============================

Quantities_ are objects in the python part that are used to trace the dependency between physical quantities and for bookkeeping, which systematic variations of a quantity exist.
Each physical quantity that is subject to systematic variations needs to be represented by such a python object. For others it can also be comfortable.
The output collection is defined as a list of such quantities and an individual branch is created in the ROOT tree for each systematic variation.

.. _Quantities: https://github.com/KIT-CMS/CROWN/blob/main/code_generation/quantities.py