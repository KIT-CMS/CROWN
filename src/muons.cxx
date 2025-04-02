#ifndef GUARD_MUONS_H
#define GUARD_MUONS_H

#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/RoccoR.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

namespace physicsobject {
namespace muon {

/**
 * @brief This function defines a new column in the dataframe that contains the 
 * corrected transverse momentum (pt) of a muon, calculated using the Rochester 
 * correction for data. 
 * 
 * The correction is taken from https://gitlab.cern.ch/akhukhun/roccor (Run2) and
 * is evaluated using the `RoccoR` correction function which is also taken from
 * this repository.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the corrected pt values
 * @param charge name of the column containing the muon charges
 * @param pt name of the column containing the muon transverse momenta
 * @param eta name of the column containing the muon eta values
 * @param phi name of the column containing the muon phi values
 * @param index_vector name of the column containing index values
 * @param position position within the index vector used to retrieve the index of
 * the wanted muon
 * @param filename file path to the Rochester correction
 * @param error_set error set number that should be used
 * @param error_member error member number that should be used
 * 
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode
PtCorrectionData(ROOT::RDF::RNode df, const std::string &outputname,
                 const std::string &charge, const std::string &pt, 
                 const std::string &eta, const std::string &phi,
                 const std::string &index_vector, const int &position, 
                 const std::string &filename, int error_set, int error_member) {
    RoccoR rc(filename);
    auto lambda = [rc, position, error_set,
                   error_member](const ROOT::RVec<int> &charges,
                                 const ROOT::RVec<float> &pts,
                                 const ROOT::RVec<float> &etas,
                                 const ROOT::RVec<float> &phis,
                                 const ROOT::RVec<int> &index_vector) {
        const int index = index_vector.at(position);
        double corrected_pt =
            pts.at(index) * rc.kScaleDT(charges.at(index), pts.at(index),
                                          etas.at(index), phis.at(index),
                                          error_set, error_member);
        Logger::get("MuonPtCorrectionData")
            ->debug("muon pt before {}, muon pt after {}",
                    pts.at(index), corrected_pt);
        return corrected_pt;
    };

    return df.Define(
        outputname, lambda,
        {charge, pt, eta, phi, index_vector});
}

/**
 * @brief This function defines a new column in the dataframe that contains the 
 * corrected transverse momentum (pt) of a muon, calculated using the Rochester 
 * correction for Monte Carlo (MC) simulation. 
 *
 * The correction is taken from https://gitlab.cern.ch/akhukhun/roccor (Run2) and
 * is evaluated using the `RoccoR` correction function which is also taken from
 * this repository.
 *
 * The correction is applied depending on whether the muon could be matched to 
 * a generator level muon or not. If the muon is matched, the correction is done 
 * using `kSpreadMC`, otherwise `kSmearMC` is used. 
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the corrected pt values
 * @param charge name of the column containing the muon charges
 * @param pt name of the column containing the muon transverse momenta
 * @param eta name of the column containing the muon eta values
 * @param phi name of the column containing the muon phi values
 * @param gen_pt name of the column containing the generator level transverse 
 * momentum of the matched muon, if no muon can be match a negative default value
 * has to be provide
 * @param n_tracker_layers name of the column containing the number of tracker layers
 * @param rndm_vector name of the column containing random values for the smearing
 * @param index_vector name of the column containing index values
 * @param position position within the index vector used to retrieve the index of
 * the wanted muon
 * @param filename file path to the Rochester correction
 * @param error_set error set number that should be used
 * @param error_member error member number that should be used
 * 
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df, const std::string &outputname,
               const std::string &charge, const std::string &pt, 
               const std::string &eta, const std::string &phi, 
               const std::string &gen_pt, const std::string &n_tracker_layers,
               const std::string &rndm_vector, const std::string &index_vector, 
               const int &position, const std::string &filename,
               const int error_set, const int error_member) {
    RoccoR rc(filename);
    auto lambda = [rc, position, error_set, error_member](
                      const ROOT::RVec<int> &charges,
                      const ROOT::RVec<float> &pts,
                      const ROOT::RVec<float> &etas,
                      const ROOT::RVec<float> &phis, 
                      const float &gen_pt,
                      const ROOT::RVec<int> &n_tracker_layers,
                      const ROOT::RVec<float> &rndm_values,
                      const ROOT::RVec<int> &index_vector) {
        double corrected_pt = default_float;
        const int index = index_vector.at(position);
        if (gen_pt > 0.) {
            corrected_pt = pts.at(index) *
                    rc.kSpreadMC(charges.at(index), pts.at(index),
                                 etas.at(index), phis.at(index), gen_pt,
                                 error_set, error_member);
        } else {
            corrected_pt = pts.at(index) *
                    rc.kSmearMC(charges.at(index), pts.at(index),
                                etas.at(index), phis.at(index),
                                n_tracker_layers.at(index),
                                rndm_values.at(position), error_set, error_member);
        }
        Logger::get("MuonPtCorrectionMC")
            ->debug("muon pt before {}, muon pt after {}",
                    pts.at(index), corrected_pt);
        return corrected_pt;
    };

    return df.Define(outputname, lambda,
                     {charge, pt, eta, phi, gen_pt, n_tracker_layers,
                      rndm_vector, index_vector});
}
} // end namespace muon
} // namespace physicsobject
#endif /* GUARD_MUONS_H */