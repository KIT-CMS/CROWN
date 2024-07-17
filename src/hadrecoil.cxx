#ifndef GUARDHADRECOIL_H
#define GUARDHADRECOIL_H

#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>

/// Namespace used for calculating quantities related to the hadronic recoil

namespace hadrecoil {

ROOT::RDF::RNode scalar_ht(ROOT::RDF::RNode df, const std::string &outputname,
                           const std::string &jet_pt, const std::string &jet_collection) {
    auto scalar_ht_func = [](const ROOT::RVec<int> &jet_pt, const ROOT::RVec<int> &jet_collection) {
        // add up the magnitudes of the pt of all selected jets
        float ht = 0;
        for (const int &jet_index : jet_collection) {
            ht += jet_pt.at(jet_index);
        }
        return ht;
    };

    return df.Define(outputname, scalar_ht_func, {jet_pt, jet_collection});
}

ROOT::RDF::RNode vectorial_mht(ROOT::RDF::RNode df, const std::string &outputname,
                               const std::string &jet_pt, const std::string &jet_eta, const std::string &jet_collection) {
    auto vectorial_mht_func = [](const ROOT::RVec<int> &jet_pt, const ROOT::RVec<int> &jet_eta, const ROOT::RVec<int> &jet_collection) {
        // add up the negative pt vectors of all selected jets 
        auto mht = ROOT::Math::PtEtaPhiEVector(0., 0., 0., 0.);
        for (const int &jet_index : jet_collection) {
            auto jet_pt_vec = ROOT::Math::PtEtaPhiEVector(jet_pt.at(jet_index), jet_eta.at(jet_index), 0., jet_pt.at(jet_index));
            mht -= jet_pt_vec;
        }
        return mht;
    };

    return df.Define(outputname, vectorial_mht_func, {jet_pt, jet_eta, jet_collection});
}

} // namespace hadrecoil
#endif /* GUARDHADRECOIL_H */
