==================================
Visualization of the Configuration
==================================

CROWN offers a visualization for configuration files. For this, the ``create_graph`` function can be used, which creates a set of JSON files that can then be visualized with the ``CROWN_visualization.html`` file.

Setup
-----
The ``create_graph`` function can be used as follows:

.. code-block:: python

    from code_generation.utility.generate_DAG import create_graph
    
    # Config definition
    # ...

    configuration.optimize()
    configuration.validate()
    configuration.report()
    if DAG_dir:
        NanoAOD_inputs = [n for n in dir(nanoAOD) if not n.startswith("__")]
        create_graph(configuration, NanoAOD_inputs, DAG_dir, "CROWNelements")
    configuration = configuration.expanded_configuration()
    return configuration

The ``create_graph`` function takes the following parameters:

* ``configuration``: The configuration object (Either a ``Configuration`` or ``FriendTreeConfiguration``).
* ``NanoAOD_inputs``: A list of strings representing inputs available from the nanoAOD files.
* ``DAG_dir``: The directory where the JSON files will be saved. It is recommended to follow the example of the template analysis and point this directory to the build directory.
* ``json_name``: The name prepended to the generated DAG JSON files.

Each call of the ``create_graph`` function performs three actions:

1. It creates or extends the set of available JSON files in the ``<DAG_dir>/DAG_files.json`` file. This includes all combinations of ``era``, ``sample type``, and ``scope`` that are present in the configuration.
2. It creates a JSON file with the DAG information for each ``scope`` in the configuration. The name of the JSON files is formatted as ``<DAG_dir>/<json_name>_<era>_<sample_type>_<scope>.json``.
3. It copies the ``CROWN_visualization.html`` file from the CROWN repository to the ``DAG_dir`` and renames it to ``index.html``.

After running the ``create_graph`` function, the DAG can be visualized by opening the directory in a web browser. The easiest ways to do this is to rely on preexisting webservers like ``https://web.etp.kit.edu/`` or to use the ``python -m http.server`` command in the terminal and navigate to the ``DAG_dir`` directory. An example for the visualization of the ``tau`` analysis can be found at ``https://web.etp.kit.edu/~tvoigtlaender/CROWN_visualization/tau/main/``.

Use with KingMaker is supported with no additional changes necessary.


Node Legend
-----------
The visualization is intended to provide a complete overview of the configuration and the dependencies between the different producers and quantities. Each node in the DAG represents a producer, with the labeled edges representing the quantity requirements between them. The width of the edges represents for how many producers a quantity is required.

Bounding boxes represent the groupings of the configuration.
* **Scopes:** Indicated by grey outlines.
* **Producer Groups:** Indicated by blue dashed boxes.
* **Vector Groups:** Indicated by purple dashed boxes.

Nodes utilize specific shapes and colors to denote their function and data sources:
* **Producers (Rectangles):** Core C++ function calls. Colored sides indicate file I/O (reads/writes). Yellow left edges represent the need for NanoAOD inputs, while blue left edges represent the need for Ntuples inputs derived from other CROWN configurations. Similarly, right blue edges illustrate that the producer contributes a quantity to the output file.
* **Filters (Orange Octagons):** Nodes that remove events based on input quantities. Therefore, they do not have output quantities and only have incoming edges.
* **Unused Producers (Pink Nodes):** Producers with no clear purpose, as they do not provide an output quantity or the produced quantity is unused. As vector producers cannot be clearly identified as producers or filters by object type alone, they have to be designated as filters during the producer definition with the ``is_filter`` parameter or will otherwise be marked as pink nodes.
* **Quantities (Yellow Labels and Edges):** Physics quantities produced by one producer and utilized by one or more subsequent producers.

Object Inspection
-----------------
By clicking on a node or edge label, a sidebar opens that provides detailed information about the corresponding producer or quantity. 

For producers this includes the producer name, the producer type, the name of the overlying group, the C++ function call it evokes, the config parameters used in the call, as well as the required input quantities and the provided output quantities. For the input quantities there is a distinction between producer and file inputs. Furthermore, in cases where quantities are read from CROWN Ntuples, the name of the original configuration is used to distinguish them. If an input is labeled as a "Filtered Source," it indicates that the data comes from a producer currently hidden by your active UI filters.

For quantities, the sidebar shows the quantity name, the producer that provides it, and the producers that require it. With both the producers and the quantities, the providing/requiring quantities/producers are interactable to jump to the respective element in the graph.

UI Controls
-----------
The UI provides several controls to customize the visualization. 

* **Dataset Selection:** The dropdown menu in the top row allow for the selection of the ``era``, ``sample type``, and ``scope`` to be visualized, based on the data in the ``DAG_files.json`` file.
* **Search & Navigation:** The regex-based search bar allows for searching for specific producers or quantities by name. For producers this also includes the producer call string and quantities read from/written to a file. The left and right arrow buttons can be used to jump between matching search results, and the cancel button can be used to clear the search.
* **Shift Selection:** The shift drop down menu allows for the selection of a systematic shift. When a shift is selected, the graph will automatically center on the affected head node. It highlights all directly affected producers, as well as all downstream producers that are affected by the shift in red. In addition the specific config parameter changed by the shift is shown in the sidebar when clicking on a producer. Downstream shifts of Friend configs take shifted input quantities into account.
* **I/O Filtering:** The input and output quantities menus show the set of inputs and outputs relevant for the current visualization. Any input/output affected by an active shift is highlighted in red. In addition, both input and output visibility can be toggled to reduce the visual complexity. Any node that does not directly interact with the toggled inputs/outputs is removed from the visualization upon clicking "Apply". 

The control UI also includes a question mark button at the top for a quick reminder of the UI controls, and a hide button to hide the control UI. Finally, the reset button in the bottom left corner can be used to reset all UI settings and return to the default visualization state.