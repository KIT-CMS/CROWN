C++ Addons
==========

In some cases, the core codebase of CROWN may not include all the features required for an analysis. To address this, users can add custom C++ code in their analysis configurations. These addons are automatically included in the generated C++ code during the code generation process.

Location and directory structure
--------------------------------

The expected structure within the analysis configuration:

.. code-block:: console

    analysis_configurations
    └── <analysis>
        └── cpp_addons
            ├── include
            │   └── <file1>.hxx
            |   └── <file2>.hxx
            |   └── ...
            └── src
                └── <file1>.cxx
                └── <file2>.cxx
                └── ...


If an analysis does not require any additional C++ code and can rely solely on the core codebase, the ``cpp_addons`` folder can be omitted entirely from the analysis configuration.

``cxx`` and ``hxx`` structure
-----------------------------

To avoid redefinition conflicts, each file should use a unique guard, especially if the addon file names are not unique w.r.t. the core codebase of CROWN. Functionalities from the core codebase can be imported using the corresponding relative path to the header files as shown in the example below.


.. code-block:: cpp

    #ifndef UNIQUE_GUARD_NAME_H  // should be unique for each file
    #define UNIQUE_GUARD_NAME_H 
    
    #include "../../../../include/utility/CorrectionManager.hxx"  // from the core codebase
    #include "../../../../include/utility/Logger.hxx"  // from the core codebase
    
    #include "ROOT/RDataFrame.hxx"  // from the ROOT framework
    #include "correction.h"

    /* Your code here */

    // End of the file
    #endif // UNIQUE_GUARD_NAME_H

For the ``hxx`` files, the same unique guard should be used as in the ``cxx`` file.

Function definitions and namespaces
-----------------------------------

All function and namespace definitions follow the standard C++ practices. If cpp addons uses one or multiple namespaces, that are also present in the core codebase of CROWN the function definitions must differ in at least their signature or name to avoid naming conflicts during compilation that will abort the compilation.
