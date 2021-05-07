Best Practices for Contributions
=================================

C++ functions
**************

The main purpose of the framework is to be effieicnt and fast. Therefore it is essential to write clear and fast C++ functions, with as little overhead as possible. We try to enforce the following minimal requirements for new functions:

* If possible, no jitting should be used. Although RDataFrames support jitted functions, this should be avoided if possible, since a jitted function can not be optimized at compile time and will slow down the execution time of the framework.

* Use `const` references wherever possible

* Documentation via docstrings directly in the code. These docstrings are then used to automatically generate the documentation.

* Check the performance using the methods described below. Try to avoid adding functions that will be "fixed later down the line". This will be the beginning of the end of the frameworks performance.


Debugging
==========

A more verbose version of the framework can be activated by setting a higher debug level. This can be done for the framework directly using

.. code-block:: cpp

    Logger::setLevel(Logger::LogLevel::DEBUG);

and for the RDataFrame using

.. code-block:: cpp

    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(
            ROOT::Detail::RDF::RDFLogChannel(),
            ROOT::Experimental::ELogLevel::kDebug + 10);

in the beginning of the `main()` function in `analysis_template.cxx`

Profiling
==========

Profiling with perf & flamegraph for CPU
*****************************************

See the script [`profiling/flamegraph.sh`](profiling/flamegraph.sh).

Running profiling on executable

.. code-block:: console

    perf record --call-graph dwarf $EXECUTABLE $INPUTFILE $OUTPUTFILE


If you want to print out the report

.. code-block:: console

    perf report

Get flamegraph repo

.. code-block:: console

    BASE_URL=https://raw.githubusercontent.com/eguiraud/FlameGraph/160b531f4c5ef0fec37e2b719ec609842a02aa99/
    # Perform the stack collapse
    curl -Os ${BASE_URL}/stackcollapse-perf.pl > stackcollapse-perf.pl


and create the flamegraph

.. code-block:: console

    perf script > out.perf
    perl stackcollapse-perf.pl out.perf > out.folded
    curl -Os ${BASE_URL}/flamegraph.pl > flamegraph.pl
    perl flamegraph.pl out.folded > flamegraph.svg


Profiling with valgrind massif for Memory
*******************************************

See the script [`profiling/massif.sh`](profiling/massif.sh).


.. code-block:: console

    valgrind --tool=massif ./a.out
    ms_print massif.out.4103388 > massif.log