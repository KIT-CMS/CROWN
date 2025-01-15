C++ Addons
==========

In some cases, the core codebase of CROWN (CROWNLIB) may not include all the features required for an analysis. To address this, users can add custom C++ code within their analysis configurations. These addons are automatically integrated to the C++ code during the code generation process.

Location and directory structure
--------------------------------

The expected structure within the analysis configuration is as follows:

.. code-block:: console

    analysis_configurations
    └── <analysis>
        └── cpp_addons
            ├── include
            │   ├── <file1>.hxx
            │   ├── <file2>.hxx
            │   └── ...
            └── src
                ├── <file1>.cxx
                ├── <file2>.cxx
                └── ...


If an analysis does not require any additional C++ code and can rely solely on CROWNLIB, the ``cpp_addons`` folder can be omitted entirely from the analysis configuration.

``.cxx`` and ``.hxx`` File structure
------------------------------------

This functionality considers files in ``analysis_configuration/<analysis>/cpp_addons/src`` and ``analysis_configuration/<analysis>/cpp_addons/include`` during compilation. The following points should be followed when adding and using custom C++ code:

* Use unique guards for each ``.cxx`` file you introduce, especially concerning CROWNLIB. For the correpsonding ``.hxx`` file(s), the same unique guard(s) should be applied.
* Use a unique function name or function signature if the custom function needs to reside in a namespace that allready exists in CROWNLIB
* Use ``../../../../include/<filename>.hxx`` if you explicetly want to import functionality from CROWNLIB. Importing CROWNLIB files using different relative paths can lead to unexpected behaviour. 

A example ``.cxx`` file could have the following structure:


.. code-block:: cpp

    #ifndef UNIQUE_GUARD_NAME_H  // unique w.r.t. CROWNLIB and other files in cpp_addons
    #define UNIQUE_GUARD_NAME_H 
    
    // Include CROWNLIB funtionalities
    #include "../../../../include/utility/CorrectionManager.hxx"
    #include "../../../../include/utility/Logger.hxx"
    
    // Feature.hxx file defined in cpp_addons
    #include "../Feature.hxx"
    
    // Globally present i.e. from the ROOT framework
    #include "ROOT/RDataFrame.hxx"
    #include "correction.h"

    /* Your code here */

    // End of the file
    #endif // UNIQUE_GUARD_NAME_H

