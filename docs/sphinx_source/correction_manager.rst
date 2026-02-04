The Correction Manager
=======================

With CROWN 0.4, we introduce the Correction Manager, a new feature that signifiantly improves the initial setup time of the dataframe when running a CROWN executable. 
The ``correctionManager::CorrectionManager`` is a global object, that is responsible for loading correction files. 
This object is always called ``correctionManager``. If a correction file is already loaded, the ``correctionManager`` will not load it again but instead return a pointer to the already loaded correction file. 
This way, the same correction file is not loaded multiple times, which can save a lot of time when running a CROWN executable.

Supported Correction Files
****************************

For now, the CorrectionManager supports the following correction files:

- correctionlib files of type ``correction::Correction`` using the :cpp:func:`correctionManager::CorrectionManager::loadCorrection` function
- correctionlib files of type ``correction::CompoundCorrection`` using the :cpp:func:`correctionManager::CorrectionManager::loadCompoundCorrection` function
- json files using the :cpp:func:`correctionManager::CorrectionManager::loadjson` function

A Documentation of all Correction Manager functions can be found in :ref:`Namespace: Correctionmanager`

Centrally provided Correction Files for CMS Analyses
******************************************************

In CMS analyses, many correction files are centrally provided by the CMS collaboration.
A list of all available correction files can be found on the CMS Analysis Corrections Wiki:
https://cms-analysis-corrections.docs.cern.ch/

The easiest way to access these correction files in CROWN is to use ``cvmfs``. 
All POGs provide their correction files on cvmfs under the path: ``/cvmfs/cms-griddata.cern.ch/cat/metadata/*``.
The exact paths for each correction can be found in the CMS Analysis Corrections Wiki linked above.

The old way of accessing the correction files via the ``jsonPOG-integration`` repository is depricated and will be removed in future versions of CROWN.

Required Changes
******************

Using the Correction Manager slightly changes the signature of CROWN functions. In the following, one simple example is shown, how to use the Correction Manager in a CROWN executable.

Python Producer Old

.. code-block:: python

   PUweights = Producer(
    name="PUweights",
    call='event::reweighting::Pileup({df}, {output}, {input}, "{PU_reweighting_file}", "{PU_reweighting_era}", "{PU_reweighting_variation}")',
    input=[nanoAOD.Pileup_nTrueInt],
    output=[q.puweight],
    scopes=["global"],
    )

Python Producer New - note the additional argument `correctionManager`

.. code-block:: python

    PUweights = Producer(
        name="PUweights",
        call='event::reweighting::Pileup({df}, correctionManager, {output}, {input}, "{PU_reweighting_file}", "{PU_reweighting_era}", "{PU_reweighting_variation}")',
        input=[nanoAOD.Pileup_nTrueInt],
        output=[q.puweight],
        scopes=["global"],
    )


In C++ calls the user has to provide the CorrectionManager object as an additional argument. In addition, the correction files are no longer loaded directly, but using one of the CorrectionManager functions. Afterwards, the correction can be used as before.

Old Function

.. code-block:: cpp

    ROOT::RDF::RNode Pileup(ROOT::RDF::RNode df, 
                            const std::string &outputname,
                            const std::string &true_pileup_number,
                            const std::string &corr_file,
                            const std::string &corr_name,
                            const std::string &variation) {

        auto evaluator =
            correction::CorrectionSet::from_file(corr_file)->at(corr_name); // old way of loading the correction file
        auto df1 =
            df.Define(outputname,
                    [evaluator, variation](const float &pu) {
                        double weight = evaluator->evaluate({pu, variation});
                        return weight;
                    },
                    {true_pileup_number});
        return df1;
    }

New Function

.. code-block:: cpp

    ROOT::RDF::RNode Pileup(ROOT::RDF::RNode df,
                            correctionManager::CorrectionManager &correction_manager,
                            const std::string &outputname,
                            const std::string &true_pileup_number,
                            const std::string &corr_file,
                            const std::string &corr_name,
                            const std::string &variation) {
        auto evaluator = correction_manager.loadCorrection(corr_file, corr_name); // new loading function
        auto df1 =
            df.Define(outputname,
                    [evaluator, variation](const float &pu) {
                        double weight = evaluator->evaluate({pu, variation});
                        return weight;
                    },
                    {true_pileup_number});
        return df1;
    }
