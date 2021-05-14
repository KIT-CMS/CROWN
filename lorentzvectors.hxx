#include "utility/Logger.hxx"
#include <Math/Vector4D.h>

/// Namespace used for lorentzvector operations

namespace lorentzvectors {

/// Function to build the lorentzvector from the pt, eta, phi and mass of a
/// particle. This utilizes the [PtEtaPhiMVector from
/// ROOT](https://root.cern/doc/master/namespaceROOT_1_1Math.html#a6cea5921731c7ac99dea921fb188df31)
///
/// \param df The input dataframe
/// \param quantities This vector contains the names of the columns:
///    -# pair - In this vector the particle index is stored. This index is used
///    to select the correct particle from the four particle quantity vectors.
///    -# pts - In this vector, the pts of the particle are stored.
///    -# etas - In this vector, the etas of the particle are stored.
///    -# phis - In this vector, the phis of the particle are stored.
///    -# masses - In this vector, the masses of the particle are stored.
///      This order has to be kept!
/// \param outputname The name of the output column in the new dataframe
/// \param position The position in the pair vector, which is used to store the
/// index of the particle in the particle quantity vectors.
///
/// \returns a new dataframe, which contains the new lorentz vector
auto buildparticle(auto df, const std::vector<std::string> quantities,
                   const std::string outputname, const int &position) {
    auto df1 = df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<float> &pts,
                   const ROOT::RVec<float> &etas, const ROOT::RVec<float> &phis,
                   const ROOT::RVec<float> &masses) {
            // the index of the particle is stored in the pair vector
            ROOT::Math::PtEtaPhiMVector p4;
            try {
                const int index = pair.at(position);
                Logger::get("lorentzvectors")->debug("pair {}", pair);
                Logger::get("lorentzvectors")->debug("pts {}", pts);
                Logger::get("lorentzvectors")->debug("etas {}", etas);
                Logger::get("lorentzvectors")->debug("phis {}", phis);
                Logger::get("lorentzvectors")->debug("masses {}", masses);
                Logger::get("lorentzvectors")->debug("Index {}", index);

                p4 = ROOT::Math::PtEtaPhiMVector(pts.at(index), etas.at(index),
                                                 phis.at(index),
                                                 masses.at(index));
            } catch (const std::out_of_range &e) {
                p4 = ROOT::Math::PtEtaPhiMVector(-999, -999, -999, -999);
                Logger::get("lorentzvectors")
                    ->debug("Index not found, retuning dummy vector !");
            }
            Logger::get("lorentzvectors")
                ->debug("P4 - Particle {} : {}", position, p4);
            return p4;
        },
        quantities);
    return df1;
}

auto build(auto df, const std::vector<std::string> &obj_quantities,
           const int pairindex, const std::string &obj_p4_name) {
    for (auto i : obj_quantities)
        Logger::get("lorentzvectors")->debug("Used object quantities {}", i);
    return lorentzvectors::buildparticle(df, obj_quantities, obj_p4_name,
                                         pairindex);
}

/// namespace used for mutau lorentzvectors
namespace mutau {
auto build(auto df, const std::string &pairname,
           const std::vector<std::string> &muon_quantities,
           const std::vector<std::string> &tau_quantities,
           const std::string &muon_p4_name, const std::string &tau_p4_name) {
    // first the muon
    // The quantities list should contain the paircolumn, and pt, eta, phi, mass
    std::vector<std::string> muons;
    muons.push_back(pairname);
    for (auto i : muon_quantities)
        muons.push_back(i);
    for (auto i : muons)
        Logger::get("lorentzvectors")->debug("Used muon quantities {}", i);
    auto df1 = lorentzvectors::buildparticle(df, muons, muon_p4_name, 0);
    // now the tau
    std::vector<std::string> taus;
    taus.push_back(pairname);
    for (auto i : tau_quantities)
        taus.push_back(i);
    for (auto i : taus)
        Logger::get("lorentzvectors")->debug("Used tau quantities {}", i);
    auto df2 = lorentzvectors::buildparticle(df1, taus, tau_p4_name, 1);
    return df2;
}
} // namespace mutau
} // namespace lorentzvectors
