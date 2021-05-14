#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>

namespace jet {
auto VetoOverlappingJets(auto df, const std::string &output_col,
                         const std::string &jet_eta, const std::string &jet_phi,
                         const std::string &p4_1, const std::string &p4_2,
                         const float deltaRmin) {
    auto df1 = df.Define(
        output_col,
        [deltaRmin](const ROOT::RVec<float> &jet_eta,
                    const ROOT::RVec<float> &jet_phi,
                    const ROOT::Math::PtEtaPhiMVector &p4_1,
                    ROOT::Math::PtEtaPhiMVector &p4_2) {
            Logger::get("VetoOverlappingJets")->debug("Checking jets");
            ROOT::RVec<int> mask(jet_eta.size(), 1);
            for (std::size_t idx = 0; idx < mask.size(); ++idx) {
                ROOT::Math::RhoEtaPhiVectorF jet(0, jet_eta[idx], jet_phi[idx]);
                Logger::get("VetoOverlappingJets")
                    ->debug("Jet:  Eta: {} Phi: {} ", jet.Eta(), jet.Phi());
                Logger::get("VetoOverlappingJets")
                    ->debug("Letpon 1 {}:  Eta: {} Phi: {}, Pt{}", p4_1,
                            p4_1.Eta(), p4_1.Phi(), p4_1.Pt());
                Logger::get("VetoOverlappingJets")
                    ->debug("Lepton 2 {}:  Eta: {} Phi: {}, Pt{}", p4_2,
                            p4_2.Eta(), p4_2.Phi(), p4_2.Pt());
                auto deltaR_1 = ROOT::Math::VectorUtil::DeltaR(jet, p4_1);
                auto deltaR_2 = ROOT::Math::VectorUtil::DeltaR(jet, p4_2);
                Logger::get("VetoOverlappingJets")
                    ->debug("DeltaR 1 {}", deltaR_1);
                Logger::get("VetoOverlappingJets")
                    ->debug("DeltaR 2 {}", deltaR_2);
                mask[idx] = (deltaR_1 > deltaRmin && deltaR_2 > deltaRmin);
            }
            return mask;
        },
        {jet_eta, jet_phi, p4_1, p4_2});
    return df1;
}

auto OrderJetsByPt(auto df, const std::string &output_col,
                   const std::string &jetmask, const std::string &jet_pt) {
    auto df1 = df.Define(
        output_col,
        [](const ROOT::RVec<int> &jetmask, const ROOT::RVec<float> &jet_pt) {
            Logger::get("OrderJetsByPt")->debug("Ordering good jets by pt");
            Logger::get("OrderJetsByPt")->debug("Jetpt before {}", jet_pt);
            Logger::get("OrderJetsByPt")->debug("Mask {}", jetmask);
            auto good_jets_pt =
                ROOT::VecOps::Where(jetmask > 0, jet_pt, (float)0.);
            Logger::get("OrderJetsByPt")->debug("Jetpt after {}", good_jets_pt);
            // we have to convert the result into an RVec of ints since argsort
            // gives back an unsigned long vector
            auto temp =
                ROOT::VecOps::Argsort(ROOT::VecOps::Nonzero(good_jets_pt));
            ROOT::RVec<int> result;
            std::transform(temp.begin(), temp.end(), result.begin(),
                           [](unsigned long x) { return (int)x; });
            return result;
        },
        {jetmask, jet_pt});
    return df1;
}
} // end namespace jet

namespace physicsobject {
namespace jet {

/// Function to filter jets based on the tau ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// idxID bitvalue of the WP the has to be passed
///
/// \return a dataframe containing the new mask
auto FilterID(auto df, const std::string maskname, const std::string nameID,
              const int idxID) {
    auto df1 = df.Define(maskname, basefunctions::FilterJetID(idxID), {nameID});
    return df1;
}
} // end namespace jet
} // end namespace physicsobject

namespace quantities {
namespace jet {
auto NumberOfJets(auto df, const std::string &outputname,
                  const std::string &jetcollection) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &jetcollection) {
                         return jetcollection.size();
                     },
                     {jetcollection});
}
} // end namespace jet
} // end namespace quantities