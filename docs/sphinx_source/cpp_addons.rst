C++ Addons
==========

In some cases, the core codebase may not include all the features required for an analysis. To address this, users can add custom C++ code in their analysis configurations. These addons are automatically included in the generated C++ code during the code generation process of CROWN.

Location and directory structure
--------------------------------

The expected structure within the analysis configuration is as follows:

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

To avoid redefinition conflicts, each file should use a unique guard, especially if the addon file names are not unique w.r.t. the core code base of CROWN. Functionalities from the core codebase can be imported using the corresponding relative path to the header files as shown in the example below.


.. code-block:: cpp

    // UNIQUE_GUARD_NAME_H should be unique for each file
    #ifndef UNIQUE_GUARD_NAME_H
    #define UNIQUE_GUARD_NAME_H 
    
    // Include the header files from the core codebase
    #include "../../../../include/utility/CorrectionManager.hxx"
    #include "../../../../include/utility/Logger.hxx"
    
    // Include the header files from the ROOT framework
    #include "ROOT/RDataFrame.hxx"
    #include "correction.h"

    /* Your code here */

    // End of the file
    #endif // UNIQUE_GUARD_NAME_H

For the ``hxx`` files, the same unique guard should be used as in the ``cxx`` file.

Namespaces within the addon ``cxx`` files do not need to be unique w.r.t. the core codebase if the functions are not conflicting in their name and signature with the CROWN definitions.
Otherwise a unique namespaces or names should be used for the addon C++ code.
