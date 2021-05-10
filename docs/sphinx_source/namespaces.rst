New C++ functions
==================


Writing a new producer consists of two parts. First, the C++ function to act on the ROOT Dataframe has to be created. In principle, the interface of the function is completely arbitrary, but a few things should be remembered, to make the new producer align with the design of the existing ones.

* The producer should live in a well defined namespace. If a producer is meant for electrons only, if should be contained in an electron namespace, rather than putting electron in the function name.

* The call of the dataframe object  :code:`(auto df,)` should always the the first argument and be the only return object of the function

* If possible, use references only for the function arguments (--> faster).

A generic new skeleton for a function looks like this:

.. code-block:: cpp

    auto GenericFunction(auto df, const std::string& parameter_1, const std::string& parameter_2,
           const float value_1) {
        auto df1 = df.Define(....);
        return df1;
    }



Documentation of Namespaces from C++
=====================================

Quantities
***********
.. doxygennamespace:: quantities
   :members:

Metfilter
***********
.. doxygennamespace:: metfilter
   :members:

Basefunctions
*************
.. doxygennamespace:: basefunctions
   :members:

Pairselection
*************
.. doxygennamespace:: pairselection
   :members:

Physicsobjects
***************
.. doxygennamespace:: physicsobject
   :members:


Lorentzvectors
***************
.. doxygenfile:: lorentzvectors.hxx
.. literalinclude:: ../../lorentzvectors.hxx
   :language: cpp
   :linenos: