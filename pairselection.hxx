#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include "utility/Logger.hxx"
#include "utility/utility.hxx"

// TODO general: enable proper logging at COMPILE (!!!) time.
//               use namespaces appropriately in functions, and use "using" to
//               make types shorter.
//
// Examples: void foo() {
//    using namespace ROOT::VecOps;
//    }
//
//    using VecF = const ROOT::RVec<float>&;

/// Namespace used for lorentzvector operations
namespace pairselection {
/**
 * @brief Function used to build a pair of GenParticles from the selected
 * DiTauPair. This uses the references of the reco particles to the gen
 * particles.
 *
 * @param df the Dataframe
 * @param recopair the column containing the DiTauPair vector
 * @param genindex_particle1 the column containing the index of the GenParticle
 * reference for the first pair particle
 * @param genindex_particle2 the column containing the index of the GenParticle
 * reference for the second pair particle
 * @param genpair name of the new column containing the GenDiTauPair
 * @return a new Dataframe with the GenDiTauPair column
 */
auto buildgenpair(auto &df, const std::string &recopair,
                  const std::string &genindex_particle1,
                  const std::string &genindex_particle2,
                  const std::string &genpair) {
    auto getGenPair = [](const ROOT::RVec<int> &recopair,
                         const ROOT::RVec<int> &genindex_particle1,
                         const ROOT::RVec<int> &genindex_particle2) {
        ROOT::RVec<int> genpair = {-1, -1};
        Logger::get("buildgenpair")->debug("existing DiTauPair: {}", recopair);
        genpair[0] = genindex_particle1[recopair.at(0)];
        genpair[1] = genindex_particle2[recopair.at(1)];
        Logger::get("buildgenpair")
            ->debug("matching GenDiTauPair: {}", genpair);
        return genpair;
    };
    return df.Define(genpair, getGenPair,
                     {recopair, genindex_particle1, genindex_particle2});
}
/// This function flags events, where a suitable particle pair is found. A
/// pair is considered suitable, if a PairSelectionAlgo (like
/// pairselection::mutau::PairSelectionAlgo) returns indices, that are not -1.
/// Events, where any of the particle indices is -1 are vetoed by this filter.
///
/// \param df The input dataframe
/// \param flagname The name of the generated flag column
/// \param pairname The name of the column, containing the particle indices
/// index of the particle in the particle quantity vectors.
///
/// \returns a dataframe with the new flag
auto flagGoodPairs(auto &df, const std::string &flagname,
                   const std::string &pairname) {
    using namespace ROOT::VecOps;
    return df.Define(
        flagname,
        [](const ROOT::RVec<int> &pair) { return bool(Min(pair) >= 0); },
        {pairname});
}

/// Function used to sort two particles based on the isolation and the pt of the
/// two particles. The function is used as the ordering function for the
/// [ROOT::VecOps::Sort()](https://root.cern.ch/doc/master/group__vecops.html#ga882439c2ff958157d2990b52dd76f599)
/// algorithm. If two quantities are the same within a given epsilon of 1e-5,
/// the next criterion is applied. The sorting is done using the following
/// criterion odering:
/// -# Isolation of the first particle
/// -# pt of the first particle
/// -# Isolation of the second particle
/// -# pt of the second particle
///
/// \param lep1pt `ROOT::RVec<float>` containing pts of the first particle
/// \param lep1iso `ROOT::RVec<float>` containing isolations of the first
/// particle \param lep2pt `ROOT::RVec<float>` containing pts of the second
/// particle \param lep2iso `ROOT::RVec<float>` containing isolations of the
/// second particle
///
/// \returns true or false based on the particle ordering.
auto compareForPairs(const ROOT::RVec<float> &lep1pt,
                     const ROOT::RVec<float> &lep1iso,
                     const ROOT::RVec<float> &lep2pt,
                     const ROOT::RVec<float> &lep2iso) {
    return [lep1pt, lep1iso, lep2pt, lep2iso](auto value_next,
                                              auto value_previous) {
        Logger::get("PairSelectionCompare")->debug("lep1 Pt: {}", lep1pt);
        Logger::get("PairSelectionCompare")->debug("lep1 Iso: {}", lep1iso);
        Logger::get("PairSelectionCompare")->debug("lep2 Pt: {}", lep2pt);
        Logger::get("PairSelectionCompare")->debug("lep2 Iso: {}", lep2iso);

        Logger::get("PairSelectionCompare")
            ->debug("Next pair: {}, {}", std::to_string(value_next.first),
                    std::to_string(value_next.second));
        Logger::get("PairSelectionCompare")
            ->debug("Previous pair: {}, {}",
                    std::to_string(value_previous.first),
                    std::to_string(value_previous.second));
        const auto i1_next = value_next.first;
        const auto i1_previous = value_previous.second;

        // start with lep1 isolation
        const auto iso1_next = lep1iso[i1_next];
        const auto iso1_previous = lep1iso[i1_previous];
        Logger::get("PairSelectionCompare")
            ->debug("Isolations: {}, {}", iso1_next, iso1_previous);
        if (not utility::ApproxEqual(iso1_next, iso1_previous)) {
            return iso1_next > iso1_previous;
        } else {
            // if too similar, compare lep1 pt
            Logger::get("PairSelectionCompare")
                ->debug("Isolation lep 1 too similar, taking pt");
            const auto pt1_next = lep1pt[i1_next];
            const auto pt1_previous = lep1pt[i1_previous];
            if (not utility::ApproxEqual(pt1_next, pt1_previous)) {
                return pt1_next > pt1_previous;
            } else {
                // if too similar, compare lep2 iso
                const auto i2_next = value_next.first;
                const auto i2_previous = value_previous.second;
                Logger::get("PairSelectionCompare")
                    ->debug("Pt lep 1 too similar, taking lep2 iso");
                const auto iso2_next = lep2iso[i2_next];
                const auto iso2_previous = lep2iso[i2_previous];
                if (not utility::ApproxEqual(iso2_next, iso2_previous)) {
                    return iso2_next > iso2_previous;
                } else {
                    // if too similar, compare lep2 pt
                    Logger::get("PairSelectionCompare")
                        ->debug("Isolation lep 2 too similar, taking pt");
                    const auto pt2_next = lep2pt[i2_next];
                    const auto pt2_previous = lep2pt[i2_previous];
                    return pt2_next > pt2_previous;
                }
            }
        }
    };
}

/// namespace for pairs in the MuTau channel
namespace mutau {

/// Implementation of the pair selection algorithm. First, only events that
/// contain at least one goodMuon and one goodTau are considered. Events contain
/// at least one good muon and one good tau, if the taumask and the mounmask
/// both have nonzero elements. These masks are constructed using the functions
/// from the physicsobject namespace (e.g. physicsobject::CutPt).
///
/// \returns an `ROOT::RVec<int>` with two values, the first one beeing the muon
/// index and the second one beeing the tau index.
auto PairSelectionAlgo() {
    Logger::get("PairSelection")->debug("Setting up algorithm");
    return [](const ROOT::RVec<float> &taupt, const ROOT::RVec<float> &tauiso,
              const ROOT::RVec<float> &muonpt, const ROOT::RVec<float> &muoniso,
              const ROOT::RVec<int> &taumask, const ROOT::RVec<int> &muonmask) {
        ROOT::RVec<int> selected_pair; // first entry is the muon index, second
                                       // entry is the tau index
        const auto original_tau_indices = ROOT::VecOps::Nonzero(taumask);
        const auto original_muon_indices = ROOT::VecOps::Nonzero(muonmask);

        if (original_tau_indices.size() == 0 or
            original_muon_indices.size() == 0) {
            selected_pair = {-1, -1};
            return selected_pair;
        }
        Logger::get("PairSelection")
            ->debug("Running algorithm on good taus and muons");

        const auto selected_taupt =
            ROOT::VecOps::Take(taupt, original_tau_indices);
        const auto selected_tauiso =
            ROOT::VecOps::Take(tauiso, original_tau_indices);
        const auto selected_muonpt =
            ROOT::VecOps::Take(muonpt, original_muon_indices);
        const auto selected_muoniso =
            ROOT::VecOps::Take(muoniso, original_muon_indices);

        const auto pair_indices = ROOT::VecOps::Combinations(
            selected_muonpt, selected_taupt); // Gives indices of mu-tau pair
        Logger::get("PairSelection")->debug("Pairs: ", pair_indices);

        // TODO, try out std::pair<UInt_t>, or std::tuple<UInt_t>.
        const auto pairs = ROOT::VecOps::Construct<std::pair<UInt_t, UInt_t>>(
            pair_indices[0], pair_indices[1]);
        Logger::get("PairSelection")->debug("Pairs size: {}", pairs.size());
        Logger::get("PairSelection")
            ->debug("Constituents pair 0: {} {}", pairs[0].first,
                    pairs[0].second);

        if (pairs.size() > 1) {
            Logger::get("PairSelection")
                ->debug("Constituents pair 1: {} {}",
                        std::to_string(pairs[1].first),
                        std::to_string(pairs[1].second));
        }

        const auto sorted_pairs = ROOT::VecOps::Sort(
            pairs, compareForPairs(selected_muonpt, -1. * selected_muoniso,
                                   selected_taupt, selected_tauiso));

        Logger::get("PairSelection")->debug("TauPt: {}", selected_taupt);
        Logger::get("PairSelection")->debug("TauIso: {}", selected_tauiso);
        Logger::get("PairSelection")->debug("MuonPt: {}", selected_muonpt);
        Logger::get("PairSelection")->debug("MuonIso: {}", selected_muoniso);

        const auto selected_mu_index = sorted_pairs[0].first;
        const auto selected_tau_index = sorted_pairs[0].second;
        selected_pair = {
            static_cast<int>(original_muon_indices[selected_mu_index]),
            static_cast<int>(original_tau_indices[selected_tau_index])};
        Logger::get("PairSelection")
            ->debug("Selected original pair indices: mu = {} , tau = {}",
                    selected_pair[0], selected_pair[1]);
        Logger::get("PairSelection")
            ->debug("MuonPt = {} , TauPt = {} ",
                    muonpt[static_cast<UInt_t>(selected_pair[0])],
                    taupt[static_cast<UInt_t>(selected_pair[1])]);
        Logger::get("PairSelection")
            ->debug("MuonIso = {} , TauIso = {} ",
                    muoniso[static_cast<UInt_t>(selected_pair[0])],
                    tauiso[static_cast<UInt_t>(selected_pair[1])]);

        return selected_pair;
    };
}

/// Function to add the Pairselection result to a dataframe
///
///  TODO add documentation here
///
/// \returns a dataframe containing the new pairname column
auto PairSelection(auto &df, const std::vector<std::string> &input_vector,
                   const std::string &pairname) {
    Logger::get("PairSelection")->debug("Setting up mutau pair building");
    auto df1 = df.Define(pairname, pairselection::mutau::PairSelectionAlgo(),
                         input_vector);
    return df1;
}

} // end namespace mutau

} // end namespace pairselection
