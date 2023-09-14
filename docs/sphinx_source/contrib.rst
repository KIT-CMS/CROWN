Writing a new producer
=======================

Writing a new producer requires two main parts, adding the :ref:`C++ function<Writing a new C++ function>` and the required :ref:`python part<Defining a new python Producer>`.

If the C++ function is written generally enough, it can be used in multiple producers and multiple purposes in the end.
For example, the producer generating the pt_1 quantity can be used regardless of what particle is being considered.

In the following, an introduction on how to add a new producer is given. As an example, we will add a new producer, which can be used to calculate the Lorentz vectors of particles, in our case electrons. For simplicity, we only want to calculate one single Lorentz vector for a given index. First, we will do the C++ implementation of the function followed by the Python definition. Keep in mind, that those two parts are connected.

Writing a new C++ function
============================

For a new C++ function, a definition in the header file, and the implementation in the source file are required. As good practice, we will add the function to a namespace called ``lorentzvector``, and call the function ``build``.
The return type of any function in CROWN should always be ``ROOT::RDF::RNode`` and the first argument of the function should always be the RDataframe, where we want to Define our new quantity. This means the basic definition of the function should look like this:

.. code-block:: cpp

    ROOT::RDF::RNode buildMet(ROOT::RDF::RNode df, ...);

This will go into the corresponding header file located within the ``include`` folder, and the implementation will go into the source file in the ``src`` folder.

The C++ functions used in CROWN always act as wrapper functions surrounding the actual RDataFrame commands that are used. Most of the time, this will be either a ``Define`` or a ``Filter``. The next step is to add all the inputs, that are needed to build the Lorentz vector. In our case, we will need the ``pt``, ``eta``, ``phi``, and ``mass`` of the electron collection from the nanoAOD and the index of which of those electrons we want to build the vector. Since this function is very generic, we could use it for all kinds of particles, not just electrons. The names of these inputs have to be added to the function. Since we just need the names of the columns, these are ``string`` arguments, and the index is an ``int`` argument.

.. code-block:: cpp

    ROOT::RDF::RNode build(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &pts, const std::string &etas,
                        const std::string &phis, const std::string &masses,
                        const int &index);

Note, that we can add these arguments as ``const string &``, since we do not need to change the value of the argument, also a reference is enough. Additionally, we added the ``outputname``, which is the name of the output quantity. This is also a string and is used to name the output quantity. The next step is to add the actual RDataFrame command. In our case, we will use the ``Define`` command, which is used to define a new quantity. The following is the definition of the command:

.. code-block:: cpp

    ROOT::RDF::RNode build(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &pts, const std::string &etas,
                        const std::string &phis, const std::string &masses,
                        const int &index) {
        return df.Define(outputname, build_vector, {pts, etas, phis, masses});
    };

Here, we add the ``Define`` command, which tells the RDataFrame that we want to define a new quantity. The first argument is the name of the output quantity, and the second argument is the name of the lambda function that is used to build the quantity. The third argument is a vector of the column names of the inputs, that we want to use within our lambda function. The parameter ``index`` is passed as a capture to the lambda, using the ``[]`` brackets. Now, the main exercise of writing a new producer is to add the implementation of the lambda function.

.. warning::
    Within the lambda is the code, that will be run on an event-by-event basis.

In our example, the lambda function is straightforward. We use our input values and calculate the Lorentz vector. First, we start with the skeleton of the lambda function:

.. code-block:: cpp

    ROOT::RDF::RNode build(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &pts, const std::string &etas,
                        const std::string &phis, const std::string &masses,
                        const int &index) {
        auto build_vector = [index](ROOT::RVec<float> pts, ROOT::RVec<float> etas,
                                    ROOT::RVec<float> phis,
                                    ROOT::RVec<float> masses) {
            // lambda code here
        };
        return df.Define(outputname, build_vector, {pts, etas, phis, masses});
    };

We define the lambda function with the type ``auto`` and pass it all our input columns. Since this is now the event-by-event part, we have to specify the correct type for all columns. In our case, since we take the quantities directly from the NanoAOD file, they all have the type ``ROOT::RVec<float>``. For the actual implementation, we can simply define a new object with the type ``ROOT::Math::PtEtaPhiMVector`` and the ``ROOT`` backend will do the rest. Since we only want to get the Lorentz vector for a single electron, we can use ``pts.at(index, -10.)`` to select the correct ``float`` value. If this value cannot be found, e.g. when there is no electron within the event, the default value ``-10`` is used. After adding this implementation, our function is complete and should look like this:

.. code-block:: cpp

    ROOT::RDF::RNode build(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &pts, const std::string &etas,
                        const std::string &phis, const std::string &masses,
                        const int &index) {
        auto build_vector = [index](ROOT::RVec<float> pts, ROOT::RVec<float> etas,
                                    ROOT::RVec<float> phis,
                                    ROOT::RVec<float> masses) {
            // Create the Lorentz vector for each particle
            Logger::get("build")->debug("size of pt {}, eta {}, phi {}, mass {}",
                                    pts.size(), etas.size(), phis.size(),
                                    masses.size());
            auto fourVec = ROOT::Math::PtEtaPhiMVector(
                pts.at(index, -10.), etas.at(index, -10.), phis.at(index, -10.),
                masses.at(index, -10.));
            return fourVec;
        };
        return df.Define(outputname, build_vector, {pts, etas, phis, masses});
    };

We also added a simple debug statement here, to print the size of the ``RVec`` objects. This concludes the implementation of the new producer.

.. warning::
    Remember to add both the definition within the header file and the implementation within the source file. Also, add docstrings to the source file as documentation of what the function does.


Defining a New Python Producer
================================

Now, after we have finished our new C++ function, we want to add it to our configuration. Therefore, we must define a new Python producer. There are several types of Producers available, more information can be found in the documentation of the Producer classes :ref:`here <Producers>`. In our example, a regular :py:class:`~code_generation.producer.Producer` is sufficient. For the complete definition of the producer, we have to define

1. The name of the producer
2. The function call, representing the mapping between the C++ function and the Python function
3. The input quantities
4. The output quantities
5. The scopes in with the producer can run

In our case, the Producer will look like this:

.. code-block:: python

    ElectronLV = Producer(
        name="ElectronLV",
        call="lorentzvector::build({df}, {output}, {input}, {electron_index_to_use})",
        input=[
            nanoAOD.Electron_pt,
            nanoAOD.Electron_eta,
            nanoAOD.Electron_phi,
            nanoAOD.Electron_mass,
        ],
        output=[q.Electron_p4],
        scopes=["ee", "em", "et"],
    )

We set the name of the Producer to be ``ElectronLV``. The call corresponds to the C++ function that is used to build the Lorentz vector. The two keywords ``input`` and ``output`` are used to specify the input and output columns. During the code generation, this will be filled with the quantities defined as the input and output of the producer.  In this example, we also use a configuration parameter called ``electron_index_to_use``. This parameter has to be defined in the configuration file and could look something like this

.. code-block:: python

    configuration.add_config_parameters(
        ["et"],
        {
            "electron_index_to_use": 0,
        },
    )


The last keyword ``scopes`` is used to specify the scopes in which the producer is available. This makes sense to prevent errors, where the producer is used in a scope that is not specified, e.g. in a final state without any electrons, we would not need to run this producer. Note that the output of this producer is of type ``ROOT::Math::PtEtaPhiMVector``, so it always makes sense to represent that in the name of the quantity in some way for easier understanding.

.. warning::
    The definition of the producer should be put into a corresponding file in the ``code_generation/producers`` directory.

The quantities themselves that are used also have to be defined. Within CROWN, systematic shifts are tracked within these quantity objects, so if a systematic shift is defined, the quantity object will register the shift. During the code generation, this allows to automatic create the necessary code to calculate all needed systematic shifts. Quantities are defined as


.. code-block:: python

    Electron_pt = NanoAODQuantity("Electron_pt")
    Electron_eta = NanoAODQuantity("Electron_eta")
    Electron_phi = NanoAODQuantity("Electron_phi")
    Electron_mass = NanoAODQuantity("Electron_mass")

.. code-block:: python

    Electron_p4s = Quantity("Electron_p4")

The only argument here is the column name of the quantity. The same goes for our new output quantity, however, since it is a new quantity it should be of type :py:class:`~code_generation.quantity.Quantity`, not :py:class:`~code_generation.quantity.NanoAODQuantity`. The quantities are defined in the files found in the ``code_generation/quantities`` directory.

After this, our new producer is now ready to be added to the configuration. To get the producer running, we have to add it to the set of producers, and we have to add the output quantity to the set of required outputs. To learn more about writing a configuration check out :ref:`Writing a CROWN Configuration<Writing a CROWN Configuration>`.

.. code-block:: python

    configuration.add_producers(
        "et",
        [
            electrons.ElectronLV,
        ],
    )

        configuration.add_outputs(
        "et",
        [
            q.Electron_p4,
        ],
    )





Best Practices for Contributions
=================================

C++ functions
**************

The main purpose of the framework is to be efficient and fast. Therefore, it is essential to write clear and fast C++ functions, with as little overhead as possible. We try to enforce the following minimal requirements for new functions:

* The producer should live in a well-defined namespace. If a producer is meant for electrons only, it should be contained in an electron namespace, rather than putting an electron in the function name.

* If possible, no jitting should be used. Although RDataFrames support jitted functions, this should be avoided if possible, since a jitted function can not be optimized at compile time and will slow down the execution time of the framework.

* Use `const` references wherever possible

* Documentation via docstrings directly in the code. These docstrings are then used to automatically generate the documentation.

* Check the performance using the methods described below. Try to avoid adding functions that will be "fixed later down the line". This will be the beginning of the end of the frameworks' performance.

* The return ``type`` of a new CROWN function should always be ``ROOT::RDF::RNode``

* The first argument of the function should always be the dataframe, again with the type ``ROOT::RDF::RNode``.

* Add meaningful debug messages to the code, using the provided logging functions.



Python Producers
*******************

There are different types of producers available

Producer: This is the standard producer class and takes the following arguments:

  * ``<string> name``: Name of the producer showing up in error messages of the Python workflow
  * ``<string> call``: Function call to be embedded into the C++ template. Use curly brackets like ``{parameter_name}`` to mark places where parameters
    of the configuration shall be written. The following keys fulfil special roles and are reserved, therefore:

    * ``{output}``: to be filled with names of output quantities (see :ref:`Python Quantities`) as strings separated by commas
    * ``{output_vec}``: like output but with curly brackets around it representing a C++ vector
    * ``{input}``: to be filled with names of input quantities as strings separated by commas
    * ``{input_vec}``: like input but with curly brackets around it representing a C++ vector
    * ``{df}``: to be filled with the input dataframe

  * ``<list of quantities> inputs``: input quantities, which are used to fill ``{input}`` and/or ``{input_vec}``. The list can be empty if no inputs are required.
  * ``<list of quantities> outputs / None``: is used to fill ``{output}`` (not usable if None). Use None (not an empty list) if no output is generated.
  * ``<list of strings> scopes``: Scopes define certain sections of the production chain. ``global`` is the initial scope, and it can be split into multiple custom scopes working on individual dataframe branches and writing out separate ROOT trees. This list of scopes defines, which scopes the producer can be used in. Dependencies between quantities will be traced separately for each scope. For example, properties of the tau candidates may be generated with the same producer but in different decay channels, which are represented by separate scopes.

VectorProducer: This is an extension of the standard producer class which can be used for C++ producers that need to be called several times with various parameter values.
  It takes the same arguments as the standard producer plus the following additional one:

  + ``<list of strings> vec_configs``: names of config parameters which contain a list of values and of which one value is supposed to substitute the corresponding placeholder in the call for each instance of the VectorProducer. Note that for VectorProducers the output argument can only be None or a list of quantities where the list must have the same length as vec_configs such that each instance will produce one of the outputs.

ProducerGroup: This object can be used to collect several producers for simplifying the configuration.
  It takes the same arguments as the standard producer plus the following additional one:

  * ``<list of producers> subproducers``: Producers can be any of the three types listed here.
    The producer group executes the subproducers first. Optionally, a closing call can be added by filling the ``call``, ``inputs``, and ``output`` arguments accordingly. If set to None, no closing call is added and only the subproducers are executed. A closing call is used to process the outputs of the subproducers forming a new output. In this case, the outputs of the subproducers can be regarded as internal quantities and be set automatically. Initialize the output of subproducers as an empty list if this automated generation of the output quantity is intended. All output quantities of the subproducers (generated automatically or by hand) are appended to the inputs of the closing call.

Python Quantities
******************

Quantities_ are objects in the Python part that are used to trace the dependency between physical quantities and for bookkeeping, in which systematic variations of a quantity exist.
Each physical quantity needs to be represented by such a Python object.
The output collection is defined as a list of such quantities and an individual branch is created in the ROOT tree for each systematic variation.

.. _Quantities: https://github.com/KIT-CMS/CROWN/blob/main/code_generation/quantity.py

Debugging
**********

A more verbose version of the framework can be activated by setting a higher debug level. This can be done by setting the argument ``-DDebug=True`` during the cmake build. This will make the code generation, as well as the executable much more verbose.

Profiling
**********

Profiling with perf & flamegraph for CPU
-----------------------------------------

See the script https://github.com/KIT-CMS/CROWN/blob/main/profiling/flamegraph.sh.


Profiling with valgrind massif for Memory
------------------------------------------

See the script https://github.com/KIT-CMS/CROWN/blob/main/profiling/massif.sh.
