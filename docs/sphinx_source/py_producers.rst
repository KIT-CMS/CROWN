Implementing new producers in the python part
=============================================

C++ functions that are supposed to be used by the code generation, referred to as producers, need a corresponding python object defined producer.py_.
There are different types of producer classes available.

.. _producer.py: https://github.com/KIT-CMS/CROWN/blob/main/code_generation/producer.py

- Producer: This is the standard producer class and takes the following arguments:

  - ``<string> name``: Name of the producer showing up in error messages of the python workflow
  - ``<string> call``: Function call to be embedded into the C++ template. Use curly brackets like ``{parameter_name}`` in order to mark places where parameters of the configuration shall be written. The following keys fulfill special roles and are reserved therefore:

    - ``{output}``: to be filled with name of output quantity as string or vector of strings, depending on number of outputs.
    - ``{input}``: to be filled with names of input quantities as strings separated by commas
    - ``{input_vec}``: like input but with curly brackets around it representing a C++ vector
    - ``{df}``: to be filled with the input dataframe

  - ``<list of quantities> inputs``: input quantities, which are used to fill ``{input}`` and/or ``{input_vec}``, list can be empty.
  - ``output``: can be None, quantity or list of quantities. Is used to fill ``{output}`` (not usable if None). A single quantity will appear as single string. A list of quantities will be filled as c++ vector of strings.
  - ``<list of strings> scopes``: Scopes define certain sections of the production chain. ``global`` is the initial scope and it can be split into multiple custom scopes working on individual dataframe branches and writing out separate ROOT trees. This list of scopes defines, in which scopes the producer can be used. Dependencies between quantities will be traced separately for each scope. For example properties of the tau candidates may be generated with the same producer but in different decay channels, which are represented by separate scopes.

- VectorProducer: This is an extension of the standard producer class which can be used for C++ producers that need to be called several times with various parameter values. It takes the same arguments as the standard producer plus the following additional one:

  - ``<list of strings> vec_configs``: names of config parameters which contain a list of values and of which one value is supposed to substitute the corresponding placeholder in the call for each instance of the VectorProducer.

  Note that for VectorProducers the output argument can only be None or a list of quantities where the list must have the same length as vec_configs such that each instance will produce one of the outputs.

- ProducerGroup: This object can be used to collect several producers for simplifying the configuration. It takes the same the same arguments as the standard producer plus the following additional one:

  - ``<list of producers> subproducers``: Producers can be any of the three types listed here. The producer group executes the subproducers first. Optionally, a closing call can be added by filling the ``call``, ``inputs``, and ``output`` arguments accordingly. If set to None, no closing call is added and only the subproducers executed. A closing call is used to process the outputs of the subproducers forming a new output. In this case, the outputs of the subproducers are regarded as internal quantities and will be set automatically. The output argument of the subproducers must be None, the internal quantities are filled in by the producer group and appended to the inputs of the closing call.
