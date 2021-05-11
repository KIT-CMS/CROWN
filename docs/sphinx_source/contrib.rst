Writing a new producer
=======================

Writing a new producer requires two main parts, adding the :ref:`C++ function<New C++ functions>` and the required :ref:`python part<Implementing new producers in the python part>`.



Best Practices for Contributions
=================================

C++ functions
**************

The main purpose of the framework is to be efficient and fast. Therefore it is essential to write clear and fast C++ functions, with as little overhead as possible. We try to enforce the following minimal requirements for new functions:

* If possible, no jitting should be used. Although RDataFrames support jitted functions, this should be avoided if possible, since a jitted function can not be optimized at compile time and will slow down the execution time of the framework.

* Use `const` references wherever possible

* Documentation via docstrings directly in the code. These docstrings are then used to automatically generate the documentation.

* Check the performance using the methods described below. Try to avoid adding functions that will be "fixed later down the line". This will be the beginning of the end of the frameworks performance.


Python producers
*****************

To be added


Debugging
**********

A more verbose version of the framework can be activated by setting a higher debug level. This can be done for the framework directly using

.. code-block:: cpp

    Logger::setLevel(Logger::LogLevel::DEBUG);

and for the RDataFrame using

.. code-block:: cpp

    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(
            ROOT::Detail::RDF::RDFLogChannel(),
            ROOT::Experimental::ELogLevel::kDebug + 10);

in the beginning of the :code:`main()` function in :code:`analysis_template.cxx`

Profiling
**********

Profiling with perf & flamegraph for CPU
-----------------------------------------

See the script https://github.com/KIT-CMS/CROWN/blob/main/profiling/flamegraph.sh.


Profiling with valgrind massif for Memory
------------------------------------------

See the script https://github.com/KIT-CMS/CROWN/blob/main/profiling/massif.sh.
