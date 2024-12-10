The Correction Manager
=======================

With CROWN 0.4, we introduce the Correction Manager, a new feature that signifiantly improves the initial setup time of the dataframe when running a CROWN executable. The ``correctionManager::CorrectionManager`` is a global object, that is responsible for loading correction files. This object is always called ``correctionManager``. If a correction file is already loaded, the ``correctionManager`` will not load it again but instead return a pointer to the already loaded correction file. This way, the same correction file is not loaded multiple times, which can save a lot of time when running a CROWN executable.


Supported Correction Files
****************************

For now, the CorrectionManager supports the following correction files:

- correctionlib files of type ``correction::Correction`` using the :cpp:func:`correctionManager::CorrectionManager::loadCorrection` function
- correctionlib files of type ``correction::CompoundCorrection`` using the :cpp:func:`correctionManager::CorrectionManager::loadCompoundCorrection` function
- json files using the :cpp:func:`correctionManager::CorrectionManager::loadjson` function

A Documentation of all Correction Manager functions can be found in :ref:`Namespace: Correctionmanager`

Required Changes
******************

Using the Correction Manager slightly changes the signature of CROWN functions. In the following, one simple example is shown, how to use the Correction Manager in a CROWN executable.

Python Producer Old

.. code-block:: python

   PUweights = Producer(
    name="PUweights",
    call='reweighting::puweights({df}, {output}, {input}, "{PU_reweighting_file}", "{PU_reweighting_era}", "{PU_reweighting_variation}")',
    input=[nanoAOD.Pileup_nTrueInt],
    output=[q.puweight],
    scopes=["global"],
    )

Python Producer New - note the additional argument `correctionManager`

.. code-block:: python

    PUweights = Producer(
        name="PUweights",
        call='reweighting::puweights({df}, correctionManager, {output}, {input}, "{PU_reweighting_file}", "{PU_reweighting_era}", "{PU_reweighting_variation}")',
        input=[nanoAOD.Pileup_nTrueInt],
        output=[q.puweight],
        scopes=["global"],
    )


In C++ calls the user has to provide the CorrectionManager object as an additional argument. In addition, the correction files are no longer loaded directly, but using one of the CorrectionManager functions. Afterwards, the correction can be used as before.

Old Function

.. code-block:: cpp

    ROOT::RDF::RNode puweights(ROOT::RDF::RNode df, const std::string &weightname,
                            const std::string &truePU,
                            const std::string &filename,
                            const std::string &eraname,
                            const std::string &variation) {

        auto evaluator =
            correction::CorrectionSet::from_file(filename)->at(eraname); // old way of loading the correction file
        auto df1 =
            df.Define(weightname,
                    [evaluator, variation](const float &pu) {
                        double weight = evaluator->evaluate({pu, variation});
                        return weight;
                    },
                    {truePU});
        return df1;
    }

New Function

.. code-block:: cpp

    ROOT::RDF::RNode
    puweights(ROOT::RDF::RNode df, correctionManager::CorrectionManager &correctionManager,
            const std::string &weightname, const std::string &truePU,
            const std::string &filename, const std::string &eraname,
            const std::string &variation) {
        auto evaluator = correctionManager.loadCorrection(filename, eraname); // new loading function
        auto df1 =
            df.Define(weightname,
                    [evaluator, variation](const float &pu) {
                        double weight = evaluator->evaluate({pu, variation});
                        return weight;
                    },
                    {truePU});
        return df1;
    }
