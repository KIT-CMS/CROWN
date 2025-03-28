#ifndef GUARD_MUONS_H
#define GUARD_MUONS_H

#include "../include/RoccoR.hxx"
#include "../include/basefilters.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
#include "correction.h"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

namespace physicsobject {
namespace muon {
/// Function to cut on muons based on the muon ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe
/// \param[in] nameID name of the ID column in the NanoAOD
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID) {
    auto df1 = df.Define(
        maskname,
        [](const ROOT::RVec<Bool_t> &id) { return (ROOT::RVec<int>)id; },
        {nameID});
    return df1;
}
/// Function to cut on muons based on the muon isolation using
/// basefilters::Max
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe
/// \param[in] Threshold maximal isolation threshold
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIsolation(ROOT::RDF::RNode df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold) {
    auto df1 = df.Define(maskname, basefilters::Max<float>(Threshold),
                         {isolationName});
    return df1;
}
/// Function to cut on muons based on the muon isolation using
/// basefilters::Min
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe
/// \param[in] Threshold minimal isolation threshold
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode AntiCutIsolation(ROOT::RDF::RNode df,
                                  const std::string &maskname,
                                  const std::string &isolationName,
                                  const float &Threshold) {
    auto df1 = df.Define(maskname, basefilters::Min<float>(Threshold),
                         {isolationName});
    return df1;
}
/// Function to cut on muons based on the muon signature: isTracker or isGlobal
///
/// \param[in] df the input dataframe
/// \param[in] isTracker name of the signature column in the NanoAOD
/// \param[in] isGlobal name of the signature column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIsTrackerOrIsGlobal(ROOT::RDF::RNode df,
                                        const std::string &isTracker,
                                        const std::string &isGlobal,
                                        const std::string &maskname) {
    auto lambda = [](const ROOT::RVec<Bool_t> &tracker,
                     const ROOT::RVec<Bool_t> &global) {
        ROOT::RVec<int> mask = (tracker == 1 || global == 1);
        Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
            ->debug("istracker {}, isglobal {}, mask {}", tracker, global,
                    mask);
        return mask;
    };
    auto df1 = df.Define(maskname, lambda, {isTracker, isGlobal});
    return df1;
}
/// Function to create a column of vector of random numbers between 0 and 1
/// with size of the input object collection
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] objCollection the name of the input object collection
/// \param[in] seed the seed of the random number generator
///
/// \return a dataframe with the new column
ROOT::RDF::RNode GenerateRndmRVec(ROOT::RDF::RNode df,
                                  const std::string &outputname,
                                  const std::string &objCollection, int seed) {
    gRandom->SetSeed(seed);
    auto lambda = [](const ROOT::RVec<int> &objects) {
        const int len = objects.size();
        float rndm[len];
        gRandom->RndmArray(len, rndm);
        ROOT::RVec<float> out = {};
        for (auto &x : rndm) {
            out.push_back(x);
        }
        return out;
    };
    return df.Define(outputname, lambda, {objCollection});
}

/// Function to create a column of Rochester correction applied transverse
/// momentum for data, see https://gitlab.cern.ch/akhukhun/roccor
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] filename the name of Rochester correction file
/// \param[in] position index of the position in the input vector
/// \param[in] objCollection the name of the column containing the input vector
/// \param[in] chargColumn the name of the column containing muon charges
/// \param[in] ptColumn the name of the column containing muon charge values
/// \param[in] etaColumn the name of the column containing muon eta values
/// \param[in] phiColumn the name of the column containing muon phi values
/// \param[in] error_set the error set number
/// \param[in] error_member the error member number
///
/// \return a dataframe with the new column
ROOT::RDF::RNode
applyRoccoRData(ROOT::RDF::RNode df, const std::string &outputname,
                const std::string &filename, const int &position,
                const std::string &objCollection,
                const std::string &chargColumn, const std::string &ptColumn,
                const std::string &etaColumn, const std::string &phiColumn,
                int error_set, int error_member) {
    RoccoR rc(filename);
    auto lambda = [rc, position, error_set,
                   error_member](const ROOT::RVec<int> &objects,
                                 const ROOT::RVec<int> &chargCol,
                                 const ROOT::RVec<float> &ptCol,
                                 const ROOT::RVec<float> &etaCol,
                                 const ROOT::RVec<float> &phiCol) {
        const int index = objects.at(position);
        double pt_rc =
            ptCol.at(index) * rc.kScaleDT(chargCol.at(index), ptCol.at(index),
                                          etaCol.at(index), phiCol.at(index),
                                          error_set, error_member);
        return pt_rc;
    };

    return df.Define(
        outputname, lambda,
        {objCollection, chargColumn, ptColumn, etaColumn, phiColumn});
}

/// Function to create a column of Rochester correction applied transverse
/// momentum for MC, see https://gitlab.cern.ch/akhukhun/roccor
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] filename the name of Rochester correction file
/// \param[in] position index of the position in the input vector
/// \param[in] objCollection the name of the column containing the input vector
/// \param[in] chargColumn the name of the column containing muon charges
/// \param[in] ptColumn the name of the column containing muon charge values
/// \param[in] etaColumn the name of the column containing muon eta values
/// \param[in] phiColumn the name of the column containing muon phi values
/// \param[in] genPtColumn the name of the column containing gen-level
/// transverse momentum value of the target muon \param[in] nTrackerLayersColumn
/// the name of the column containing number of tracker layers values \param[in]
/// rndmColumn the name of the column containing random number generated for
/// each muon \param[in] error_set the error set number \param[in] error_member
/// the error member number
///
/// \return a dataframe with the new column
ROOT::RDF::RNode
applyRoccoRMC(ROOT::RDF::RNode df, const std::string &outputname,
              const std::string &filename, const int &position,
              const std::string &objCollection, const std::string &chargColumn,
              const std::string &ptColumn, const std::string &etaColumn,
              const std::string &phiColumn, const std::string &genPtColumn,
              const std::string &nTrackerLayersColumn,
              const std::string &rndmColumn, int error_set, int error_member) {
    RoccoR rc(filename);
    auto lambda = [rc, position, error_set, error_member](
                      const ROOT::RVec<int> &objects,
                      const ROOT::RVec<int> &chargCol,
                      const ROOT::RVec<float> &ptCol,
                      const ROOT::RVec<float> &etaCol,
                      const ROOT::RVec<float> &phiCol, const float &genPt,
                      const ROOT::RVec<int> &nTrackerLayersCol,
                      const ROOT::RVec<float> &rndmCol) {
        double pt_rc = default_float;
        const int index = objects.at(position);
        if (genPt > 0.) {
            pt_rc = ptCol.at(index) *
                    rc.kSpreadMC(chargCol.at(index), ptCol.at(index),
                                 etaCol.at(index), phiCol.at(index), genPt,
                                 error_set, error_member);
        } else {
            pt_rc = ptCol.at(index) *
                    rc.kSmearMC(chargCol.at(index), ptCol.at(index),
                                etaCol.at(index), phiCol.at(index),
                                nTrackerLayersCol.at(index),
                                rndmCol.at(position), error_set, error_member);
        }

        return pt_rc;
    };

    return df.Define(outputname, lambda,
                     {objCollection, chargColumn, ptColumn, etaColumn,
                      phiColumn, genPtColumn, nTrackerLayersColumn,
                      rndmColumn});
}
} // end namespace muon
} // namespace physicsobject
#endif /* GUARD_MUONS_H */