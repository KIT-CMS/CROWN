#ifndef GUARD_LORENTZVECTOR_H
#define GUARD_LORENTZVECTOR_H

#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>

namespace lorentzvector {

/**
 * @brief This function builds a Lorentz vector for a single object with given 
 * kinematic values for \f$p_T\f$, \f$\eta\f$, \f$\phi\f$ and mass. The object 
 * is specified by an index from an `index_vector`. This function utilizes the 
 * [PtEtaPhiMVector from ROOT](https://root.cern/doc/master/namespaceROOT_1_1Math.html#a6cea5921731c7ac99dea921fb188df31).
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the Lorentz vector
 * @param pt name of the column containing the \f$p_T\f$ values of an object 
 * for the Lorentz vector
 * @param eta name of the column containing the \f$\eta\f$ values of an object 
 * for the Lorentz vector
 * @param phi name of the column containing the \f$\phi\f$ values of an object 
 * for the Lorentz vector
 * @param mass name of the column containing the mass values of an object 
 * for the Lorentz vector
 * @param index_vector name of the column containing indices of objects
 * @param position position in the index vector that specifies for which object 
 * in the object vector the Lorentz vector should be build
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Build(ROOT::RDF::RNode df,
                       const std::string &outputname,
                       const std::string &pt,
                       const std::string &eta,
                       const std::string &phi,
                       const std::string &mass,
                       const std::string &index_vector,
                       const int position) {
    Logger::get("lorentzvector::Build")->debug("Building {}", outputname);
    auto df1 = df.Define(
        outputname,
        [position](const ROOT::RVec<float> &pts, const ROOT::RVec<float> &etas, const ROOT::RVec<float> &phis,
            const ROOT::RVec<float> &masses, const ROOT::RVec<int> &index_vector) {
            ROOT::Math::PtEtaPhiMVector p4;
            try {
                const int index = index_vector.at(position);
                Logger::get("lorentzvector::Build")->debug("Used object {} with quantities pt {}, eta {}, phi {}, mass {}", 
                                                index, pts.at(index), etas.at(index), phis.at(index), masses.at(index));
                p4 = ROOT::Math::PtEtaPhiMVector(pts.at(index), etas.at(index), phis.at(index), masses.at(index));
            }
            catch (const std::out_of_range &e) {
                p4 = ROOT::Math::PtEtaPhiMVector(default_float, default_float,
                                                 default_float, default_float);
            }
            Logger::get("lorentzvector::Build")
                ->debug("P4 : {}", p4);
            return p4;
        },
        {pt, eta, phi, mass, index_vector});
    return df1;
}

/**
 * @brief This function builds a Lorentz vector for a single object with given 
 * kinematic values for \f$p_T\f$, \f$\eta\f$, \f$\phi\f$ and mass. The object is specified
 * by an `index`. This function utilizes the [PtEtaPhiMVector from ROOT](https://root.cern/doc/master/namespaceROOT_1_1Math.html#a6cea5921731c7ac99dea921fb188df31).
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the Lorentz vector
 * @param pt name of the column containing the \f$p_T\f$ values of an object 
 * for the Lorentz vector
 * @param eta name of the column containing the \f$\eta\f$ values of an object 
 * for the Lorentz vector
 * @param phi name of the column containing the \f$\phi\f$ values of an object 
 * for the Lorentz vector
 * @param mass name of the column containing the mass values of an object 
 * for the Lorentz vector
 * @param index index of an object that specifies for which object in the object vector 
 * the Lorentz vector should be build
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Build(ROOT::RDF::RNode df,
                       const std::string &outputname,
                       const std::string &pt,
                       const std::string &eta,
                       const std::string &phi,
                       const std::string &mass,
                       const int index) {
    Logger::get("lorentzvector::Build")->debug("Building {}", outputname);
    auto df1 = df.Define(
        outputname,
        [index](const ROOT::RVec<float> &pts, const ROOT::RVec<float> &etas, const ROOT::RVec<float> &phis,
            const ROOT::RVec<float> &masses) {
            ROOT::Math::PtEtaPhiMVector p4;
            try {
                Logger::get("lorentzvector::Build")->debug("Used object {} with quantities pt {}, eta {}, phi {}, mass {}", 
                                                index, pts.at(index), etas.at(index), phis.at(index), masses.at(index));
                p4 = ROOT::Math::PtEtaPhiMVector(pts.at(index), etas.at(index), phis.at(index), masses.at(index));
            }
            catch (const std::out_of_range &e) {
                p4 = ROOT::Math::PtEtaPhiMVector(default_float, default_float,
                                                 default_float, default_float);
            }
            Logger::get("lorentzvector::Build")
                ->debug("P4 : {}", p4);
            return p4;
        },
        {pt, eta, phi, mass});
    return df1;
}

/**
 * @brief This function builds a Lorentz vector for a single particle with given 
 * kinematic values for \f$p_T\f$, \f$\eta\f$, \f$\phi\f$ and mass. This function 
 * utilizes the [PtEtaPhiMVector from ROOT](https://root.cern/doc/master/namespaceROOT_1_1Math.html#a6cea5921731c7ac99dea921fb188df31).
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the Lorentz vector
 * @param pt name of the column containing the \f$p_T\f$ value for the Lorentz vector
 * @param eta name of the column containing the \f$\eta\f$ value for the Lorentz vector
 * @param phi name of the column containing the \f$\phi\f$ value for the Lorentz vector
 * @param mass name of the column containing the mass value for the Lorentz vector
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Build(ROOT::RDF::RNode df,
                       const std::string &outputname,
                       const std::string &pt,
                       const std::string &eta,
                       const std::string &phi,
                       const std::string &mass) {
    Logger::get("lorentzvector::Build")->debug("Building {}", outputname);
    auto df1 = df.Define(
        outputname,
        [](const float &pt, const float &eta, const float &phi,
            const float &mass) {
            Logger::get("lorentzvector::Build")->debug("Used object quantities pt {}, eta {}, phi {}, mass {}", 
                                                pt, eta, phi, mass);
            ROOT::Math::PtEtaPhiMVector p4;
            if (pt >= 0.) {
                p4 = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
            }
            else {
                p4 = ROOT::Math::PtEtaPhiMVector(default_float, default_float,
                                                 default_float, default_float);
            }
            Logger::get("lorentzvector::Build")
                ->debug("P4 : {}", p4);
            return p4;
        },
        {pt, eta, phi, mass});
    return df1;
}

/**
 * @brief This function builds a Lorentz vector for the missing transverse 
 * momentum/energy (MET) with given kinematic values for \f$p_T\f$ and \f$\phi\f$. 
 * For MET, \f$\eta\f$ and the mass are assumed to be zero. This function 
 * utilizes the [PtEtaPhiMVector from ROOT](https://root.cern/doc/master/namespaceROOT_1_1Math.html#a6cea5921731c7ac99dea921fb188df31).
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the missing transverse
 * energy Lorentz vector
 * @param met_pt name of the column containing the \f$p_T\f$ value for the Lorentz vector
 * @param met_phi name of the column containing the \f$\phi\f$ value for the Lorentz vector
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode BuildMET(ROOT::RDF::RNode df, const std::string &outputname,
                            const std::string &met_pt, const std::string &met_phi) {
    auto construct_metvector = [](const float &pt, const float &phi) {
        // for Met, eta and mass is zero
        auto met = ROOT::Math::PtEtaPhiMVector(pt, 0., phi, 0.);
        return met;
    };
    return df.Define(outputname, construct_metvector, {met_pt, met_phi});
}

/**
 * @brief Build a new column with a collection of Lorentz vectors per dataframe entry, which are
 * created from the \f$p_T\f$, \f$\eta\f$, \f$\phi\f$ and mass columns of a collection (e.g. jets or muons).
 *
 * For instance, this can be used to create four-vector objects from the four-vector component
 * columns of a NanoAOD collection.
 *
 * The function expects \f$p_T\f$, \f$\eta\f$, \f$\phi\f$ and mass columns, which contain `ROOT::VecOps::RVec<float>`
 * objects. In addition, the argument `object_mask` must point to a column which contains
 * index lists of type `ROOT::VecOps::RVec<int>`. These index lists contain the indices of
 * elements, for which four-vectors should be build. The output column contains a vector of
 * four-vectors `ROOT::VecOps::RVec<PtEtaPhiMVector>` and only contains
 * the elements, which have been selected in by applying the `object_mask`.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the collection vector 
 * of Lorentz vectors
 * @param pt name of the column containing the \f$p_T\f$ values of the object collection
 * @param eta name of the column containing the \f$\eta\f$ values of the object collection
 * @param phi name of the column containing the \f$\phi\f$ values of the object collection
 * @param mass name of the column containing the mass values of the object collection
 * @param object_mask name of the mask column indicating selected objects to only calculate
 * the four-vectors for selected objects
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode BuildCollection(ROOT::RDF::RNode df,
                                   const std::string &outputname,
                                   const std::string &pt,
                                   const std::string &eta,
                                   const std::string &phi,
                                   const std::string &mass,
                                   const std::string &object_mask) {
    auto build_p4_collection = [](const ROOT::RVec<float> &pts,
                                  const ROOT::RVec<float> &etas,
                                  const ROOT::RVec<float> &phis,
                                  const ROOT::RVec<float> &masses,
                                  const ROOT::RVec<int> &object_mask) {
        auto pts_ = ROOT::VecOps::Take(pts, object_mask);
        auto etas_ = ROOT::VecOps::Take(etas, object_mask);
        auto phis_ = ROOT::VecOps::Take(phis, object_mask);
        auto masses_ = ROOT::VecOps::Take(masses, object_mask);
        return ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(pts_, etas_, phis_, masses_);
    };
    return df.Define(outputname, build_p4_collection, {pt, eta, phi, mass, object_mask});
}

/** 
 * @brief Build a new column with a collection of Lorentz vectors per dataframe entry, which are
 * created from the \f$p_T\f$, \f$\eta\f$, \f$\phi\f$ and mass columns of a collection (e.g. jets or muons).
 *
 * For instance, this can be used to create four-vector objects from the four-vector component
 * columns of a NanoAOD collection.
 *
 * The function expects \f$p_T\f$, \f$\eta\f$, \f$\phi\f$ and mass columns, which contain `ROOT::VecOps::RVec<float>`
 * objects. The output column contains four-vectors `ROOT::VecOps::RVec<PtEtaPhiMVector>`.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the collection vector 
 * of Lorentz vectors
 * @param pt name of the column containing the \f$p_T\f$ values of the object collection
 * @param eta name of the column containing the \f$\eta\f$ values of the object collection
 * @param phi name of the column containing the \f$\phi\f$ values of the object collection
 * @param mass name of the column containing the mass values of the object collection
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode BuildCollection(ROOT::RDF::RNode df,
                                const std::string &outputname,
                                const std::string &pt,
                                const std::string &eta,
                                const std::string &phi,
                                const std::string &mass) {
    auto build_p4_collection = [](const ROOT::RVec<float> &pts,
                                  const ROOT::RVec<float> &etas,
                                  const ROOT::RVec<float> &phis,
                                  const ROOT::RVec<float> &masses) {
        return ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(pts, etas, phis, masses);
    };
    return df.Define(outputname, build_p4_collection, {pt, eta, phi, mass});
}

/**
 * @brief This function is scaling the \f$p_T\f$ and mass (therefore also energy)
 * of a Lorentz vector by a given scale factor.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the scaled Lorentz vector 
 * @param vector name of the column containing the Lorentz vector to be scaled
 * @param scalefactor scale factor value that should be applied
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Scale(ROOT::RDF::RNode df, const std::string &outputname,
                         const std::string &vector,
                         const float &scalefactor) {
    return df.Define(
        outputname,
        [scalefactor](const ROOT::Math::PtEtaPhiMVector &p4) {
            if (p4.pt() < 0.0)
                return ROOT::Math::PtEtaPhiMVector(
                    default_float, default_float, default_float, default_float);
            auto p4_scaled = p4 * scalefactor;
            return p4_scaled;
        },
        {vector});
}

/**
 * @brief This function get the \f$p_T\f$ value of a Lorentz vector and saves
 * it into a new column.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the \f$p_T\f$ value 
 * of a Lorentz vector 
 * @param vector name of the column containing the Lorentz vector
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode GetPt(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &vector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.pt(); },
        {vector});
}

/**
 * @brief This function get the \f$\eta\f$ value of a Lorentz vector and saves
 * it into a new column.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the \f$\eta\f$ value 
 * of a Lorentz vector 
 * @param vector name of the column containing the Lorentz vector
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode GetEta(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &vector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.eta(); },
        {vector});
}

/**
 * @brief This function get the \f$\phi\f$ value of a Lorentz vector and saves
 * it into a new column.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the \f$\phi\f$ value 
 * of a Lorentz vector 
 * @param vector name of the column containing the Lorentz vector
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode GetPhi(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &vector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.phi(); },
        {vector});
}

/**
 * @brief This function get the mass value of a Lorentz vector and saves
 * it into a new column.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the mass value 
 * of a Lorentz vector 
 * @param vector name of the column containing the Lorentz vector
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode GetMass(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::string &vector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.mass(); },
        {vector});
}

/**
 * @brief This function get the energy value of a Lorentz vector and saves
 * it into a new column.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the energy value 
 * of a Lorentz vector 
 * @param vector name of the column containing the Lorentz vector
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode GetEnergy(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::string &vector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.energy(); },
        {vector});
}
} // end namespace lorentzvector
#endif /* GUARD_LORENTZVECTOR_H */