#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include "bitset"
#include "utility/Logger.hxx"
#include "utility/utility.hxx"

typedef std::bitset<15> StatusBits;
/**
 * @brief A struct to store the information of a generator pair for easier
 * access
 *
 */
struct GenParticle {
    int index;
    int status;
    std::bitset<15> statusflag;
    int pdgid;
    int motherid;
};
/**
 * @brief Function used to propagate back the mother particles of a given
 * particle, and check if a mother particle with a given PDG ID is found.
 *
 * @param genparticles An RVec of GenParticle objects.
 * @param index The index of the particle to be checked.
 * @param mother_pdgid The PDG ID of the mother particle to be checked.
 * @return True if the mother particle with the given PDG ID is found, false
 * otherwise.
 */
bool check_mother(ROOT::RVec<GenParticle> genparticles, const int index,
                  const int mother_pdgid) {
    GenParticle mother = genparticles.at(genparticles.at(index).motherid);
    Logger::get("check_mother")->debug("Testing particle: {}", index);
    Logger::get("check_mother")->debug("-->Mother PDGID: {}", mother.pdgid);
    if (mother.pdgid == mother_pdgid) {
        Logger::get("check_mother")->debug("-->found ");
        return true;
    } else if (mother.motherid == -1) {
        Logger::get("check_mother")->debug("--> no compatible mother found");
        return false;
    } else {
        Logger::get("check_mother")->debug("going deeper.... ");
        return check_mother(genparticles, mother.index, mother_pdgid);
    }
}

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
        genpair[0] = genindex_particle1.at(recopair.at(0), -1);
        genpair[1] = genindex_particle2.at(recopair.at(1), -1);
        Logger::get("buildgenpair")
            ->debug("matching GenDiTauPair: {}", genpair);
        return genpair;
    };
    return df.Define(genpair, getGenPair,
                     {recopair, genindex_particle1, genindex_particle2});
}

/**
 * @brief Function to get the true gen-level di-tau pair from the event. The
 * pair is build by searching for the gen mother particle and the two requested
 * daughter particles. For each stable daughter particle found in the collection
 of gen particles, it is checked if the particle is a daughter of the requested
 mother particle.
 *
 * @param df the Dataframe
 * @param statusflags the column containing the status flags of the gen
 particles
 * @param status the column containing the status of the genparticles (status=1
 means stable)
 * @param pdgids the column containing the PDGID of the gen particles
 * @param motherids the column containing the index of the mother particle of
 the gen particles
 * @param pts the column containing the pt of the gen particles (used for
 sorting the particles by pt)
 * @param genpair the output column containing the index of the two selected gen
 particles
 * @param mother_pdgid the PDGID of the mother particle
 * @param daughter_1_pdgid the PDGID of the first daughter particle
 * @param daughter_2_pdgid the PDGID of the second daughter particle
 * @return auto the new Dataframe with the genpair column
 */
auto buildtruegenpair(auto &df, const std::string &statusflags,
                      const std::string &status, const std::string &pdgids,
                      const std::string &motherids, const std::string &pts,
                      const std::string &genpair, const int mother_pdgid,
                      const int daughter_1_pdgid, const int daughter_2_pdgid) {

    auto getTrueGenPair = [mother_pdgid, daughter_1_pdgid,
                           daughter_2_pdgid](const ROOT::RVec<int> &statusflags,
                                             const ROOT::RVec<int> &status,
                                             const ROOT::RVec<int> &pdgids,
                                             const ROOT::RVec<int> &motherids,
                                             const ROOT::RVec<float> &pts) {
        ROOT::RVec<int> genpair = {-1, -1};

        // first we build structs, one for each genparticle
        Logger::get("buildtruegenpair")
            ->debug("Starting to build True Genpair for event");
        ROOT::RVec<GenParticle> genparticles;
        for (int i = 0; i < statusflags.size(); ++i) {
            GenParticle genparticle;
            genparticle.index = i;
            genparticle.status = status.at(i);
            genparticle.statusflag = StatusBits(statusflags.at(i));
            genparticle.pdgid = abs(pdgids.at(i));
            genparticle.motherid = motherids.at(i);
            genparticles.push_back(genparticle);
        }
        Logger::get("buildtruegenpair")->debug("genparticles: ");
        for (auto &genparticle : genparticles) {
            Logger::get("buildtruegenpair")
                ->debug("|--------------------------------------------");
            Logger::get("buildtruegenpair")
                ->debug("|    Index: {}", genparticle.index);
            Logger::get("buildtruegenpair")
                ->debug("|    Status: {}", genparticle.status);
            Logger::get("buildtruegenpair")
                ->debug("|     Statusflag: {}", genparticle.statusflag);
            Logger::get("buildtruegenpair")
                ->debug("|    Pdgid: {}", genparticle.pdgid);
            Logger::get("buildtruegenpair")
                ->debug("|    motherid: {}", genparticle.motherid);
            Logger::get("buildtruegenpair")
                ->debug("|--------------------------------------------");
        }
        auto gen_candidates_1 = ROOT::VecOps::Filter(
            genparticles, [daughter_1_pdgid](const GenParticle &genparticle) {
                return (daughter_1_pdgid == genparticle.pdgid) &&
                       (genparticle.status == 1);
            });
        auto gen_candidates_2 = ROOT::VecOps::Filter(
            genparticles, [daughter_2_pdgid](const GenParticle &genparticle) {
                return (daughter_2_pdgid == genparticle.pdgid) &&
                       (genparticle.status == 1);
            });
        if (daughter_1_pdgid == daughter_2_pdgid) {
            // if the two daughter particles are the same, we need to find
            // two valid ones there
            for (const auto &gen_candidate_1 : gen_candidates_1) {
                bool found = check_mother(genparticles, gen_candidate_1.index,
                                          mother_pdgid);
                Logger::get("buildtruegenpair")
                    ->debug("Checking Daughter Candidates");
                Logger::get("buildtruegenpair")
                    ->debug("|--------------------------------------------");
                Logger::get("buildtruegenpair")
                    ->debug("|    Index: {}", gen_candidate_1.index);
                Logger::get("buildtruegenpair")
                    ->debug("|    Status: {}", gen_candidate_1.status);
                Logger::get("buildtruegenpair")
                    ->debug("|     Statusflag: {}", gen_candidate_1.statusflag);
                Logger::get("buildtruegenpair")
                    ->debug("|    Pdgid: {}", gen_candidate_1.pdgid);
                Logger::get("buildtruegenpair")
                    ->debug("|    motherid: {}", gen_candidate_1.motherid);
                Logger::get("buildtruegenpair")
                    ->debug("|    found_correct_mother: {}", found);
                Logger::get("buildtruegenpair")
                    ->debug("|    motherPdgid: {}", mother_pdgid);
                Logger::get("buildtruegenpair")
                    ->debug("|--------------------------------------------");
                if (found && genpair[0] == -1) {
                    genpair[0] = gen_candidate_1.index;
                } else if (found && genpair[1] == -1) {
                    genpair[1] = gen_candidate_1.index;
                }
            }
        } else {
            // in the case of two different daughter particles, we need find
            // each one seperately
            for (const auto &gen_candidate_1 : gen_candidates_1) {
                if (check_mother(genparticles, gen_candidate_1.index,
                                 mother_pdgid)) {
                    genpair[0] = gen_candidate_1.index;
                    break;
                }
            }
            for (const auto &gen_candidate_2 : gen_candidates_2) {
                if (check_mother(genparticles, gen_candidate_2.index,
                                 mother_pdgid)) {
                    genpair[1] = gen_candidate_2.index;
                    break;
                }
            }
        }
        Logger::get("buildtruegenpair")
            ->debug("Selected Particles: {} {}", genpair[0], genpair[1]);
        if (genpair[0] == -1 || genpair[1] == -1) {
            Logger::get("buildtruegenpair")
                ->debug("no viable daughter particles found");
            return genpair;
        }
        if (daughter_1_pdgid == daughter_2_pdgid) {
            if (pts.at(genpair[0]) < pts.at(genpair[1])) {
                std::swap(genpair[0], genpair[1]);
            }
        }

        return genpair;
    };
    return df.Define(genpair, getTrueGenPair,
                     {statusflags, status, pdgids, motherids, pts});
}
/// This function flags events, where a suitable particle pair is found.
/// A pair is considered suitable, if a PairSelectionAlgo (like
/// pairselection::mutau::PairSelectionAlgo) returns indices, that are
/// not -1. Events, where any of the particle indices is -1 are vetoed
/// by this filter.
///
/// \param df The input dataframe
/// \param flagname The name of the generated flag column
/// \param pairname The name of the column, containing the indices of the
/// particles in the particle quantity vectors.
/// \returns a dataframe with the
/// new flag
auto flagGoodPairs(auto &df, const std::string &flagname,
                   const std::string &pairname) {
    using namespace ROOT::VecOps;
    return df.Define(
        flagname,
        [](const ROOT::RVec<int> &pair) { return bool(Min(pair) >= 0); },
        {pairname});
}

/// Function used to sort two particles based on the isolation and the
/// pt of the two particles. The function is used as the ordering
/// function for the
/// [ROOT::VecOps::Sort()](https://root.cern.ch/doc/master/group__vecops.html#ga882439c2ff958157d2990b52dd76f599)
/// algorithm. If two quantities are the same within a given epsilon of
/// 1e-5, the next criterion is applied. The sorting is done using the
/// following criterion odering:
/// -# Isolation of the first particle
/// -# pt of the first particle
/// -# Isolation of the second particle
/// -# pt of the second particle
///
/// \param lep1pt `ROOT::RVec<float>` containing pts of the first
/// particle
/// \param lep1iso `ROOT::RVec<float>` containing isolations of
/// the first particle
/// \param lep2pt `ROOT::RVec<float>` containing pts
/// of the second particle
/// \param lep2iso `ROOT::RVec<float>` containing
/// isolations of the second particle
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
        bool result = false;
        Logger::get("PairSelectionCompare")
            ->debug("Next pair: {}, {}", std::to_string(value_next.first),
                    std::to_string(value_next.second));
        Logger::get("PairSelectionCompare")
            ->debug("Previous pair: {}, {}",
                    std::to_string(value_previous.first),
                    std::to_string(value_previous.second));
        const auto i1_next = value_next.first;
        const auto i1_previous = value_previous.first;
        Logger::get("PairSelectionCompare")
            ->debug("i1_next: {}, i1_previous : {}", i1_next, i1_previous);
        // start with lep1 isolation
        const auto iso1_next = lep1iso.at(i1_next);
        const auto iso1_previous = lep1iso.at(i1_previous);
        Logger::get("PairSelectionCompare")
            ->debug("Isolations: {}, {}", iso1_next, iso1_previous);
        if (not utility::ApproxEqual(iso1_next, iso1_previous)) {
            result = iso1_next > iso1_previous;
        } else {
            // if too similar, compare lep1 pt
            Logger::get("PairSelectionCompare")
                ->debug("Isolation lep 1 too similar, taking pt 1");
            const auto pt1_next = lep1pt.at(i1_next);
            const auto pt1_previous = lep1pt.at(i1_previous);
            if (not utility::ApproxEqual(pt1_next, pt1_previous)) {
                result = pt1_next > pt1_previous;
            } else {
                // if too similar, compare lep2 iso
                const auto i2_next = value_next.second;
                const auto i2_previous = value_previous.second;
                Logger::get("PairSelectionCompare")
                    ->debug("Pt lep 1 too similar, taking lep2 iso");
                const auto iso2_next = lep2iso.at(i2_next);
                const auto iso2_previous = lep2iso.at(i2_previous);
                if (not utility::ApproxEqual(iso2_next, iso2_previous)) {
                    result = iso2_next > iso2_previous;
                } else {
                    // if too similar, compare lep2 pt
                    Logger::get("PairSelectionCompare")
                        ->debug("Isolation lep 2 too similar, taking pt 2");
                    const auto pt2_next = lep2pt.at(i2_next);
                    const auto pt2_previous = lep2pt.at(i2_previous);
                    result = pt2_next > pt2_previous;
                }
            }
        }
        Logger::get("PairSelectionCompare")
            ->debug("Returning result {}", result);
        return result;
    };
}

/// namespace for pairs in the MuTau channel
namespace semileptonic {

/// Implementation of the pair selection algorithm. First, only events
/// that contain at least one goodlepton and one goodTau are considered.
/// Events contain at least one good lepton and one good tau, if the
/// tau_mask and the mounmask both have nonzero elements. These masks are
/// constructed using the functions from the physicsobject namespace
/// (e.g. physicsobject::CutPt).
///
/// \returns an `ROOT::RVec<int>` with two values, the first one beeing
/// the lepton index and the second one beeing the tau index.
auto PairSelectionAlgo(const float &mindeltaR) {
    Logger::get("semileptonic::PairSelectionAlgo")
        ->debug("Setting up algorithm");
    return [mindeltaR](const ROOT::RVec<float> &tau_pt,
                       const ROOT::RVec<float> &tau_eta,
                       const ROOT::RVec<float> &tau_phi,
                       const ROOT::RVec<float> &tau_mass,
                       const ROOT::RVec<float> &tau_iso,
                       const ROOT::RVec<float> &lepton_pt,
                       const ROOT::RVec<float> &lepton_eta,
                       const ROOT::RVec<float> &lepton_phi,
                       const ROOT::RVec<float> &lepton_mass,
                       const ROOT::RVec<float> &lepton_iso,
                       const ROOT::RVec<int> &tau_mask,
                       const ROOT::RVec<int> &lepton_mask) {
        // first entry is the lepton index,
        // second entry is the tau index
        ROOT::RVec<int> selected_pair = {-1, -1};
        const auto original_tau_indices = ROOT::VecOps::Nonzero(tau_mask);
        const auto original_lepton_indices = ROOT::VecOps::Nonzero(lepton_mask);

        if (original_tau_indices.size() == 0 or
            original_lepton_indices.size() == 0) {
            return selected_pair;
        }
        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Running algorithm on good taus and leptons");

        const auto selected_tau_pt =
            ROOT::VecOps::Take(tau_pt, original_tau_indices);
        const auto selected_tau_iso =
            ROOT::VecOps::Take(tau_iso, original_tau_indices);
        const auto selected_lepton_pt =
            ROOT::VecOps::Take(lepton_pt, original_lepton_indices);
        const auto selected_lepton_iso =
            ROOT::VecOps::Take(lepton_iso, original_lepton_indices);

        const auto pair_indices = ROOT::VecOps::Combinations(
            selected_lepton_pt,
            selected_tau_pt); // Gives indices of mu-tau pair
        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Pairs: {} {}", pair_indices[0], pair_indices[1]);

        const auto pairs = ROOT::VecOps::Construct<std::pair<UInt_t, UInt_t>>(
            pair_indices[0], pair_indices[1]);
        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Pairs size: {}", pairs.size());
        int counter = 0;
        for (auto &pair : pairs) {
            counter++;
            Logger::get("semileptonic::PairSelectionAlgo")
                ->debug("Constituents pair {}. : {} {}", counter, pair.first,
                        pair.second);
        }

        const auto sorted_pairs = ROOT::VecOps::Sort(
            pairs,
            compareForPairs(selected_lepton_pt, -1. * selected_lepton_iso,
                            selected_tau_pt, selected_tau_iso));

        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Original TauPt: {}", tau_pt);
        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Original TauIso: {}", tau_iso);
        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Original leptonPt: {}", lepton_pt);
        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Original leptonIso: {}", lepton_iso);

        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Selected TauPt: {}", selected_tau_pt);
        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Selected TauIso: {}", selected_tau_iso);
        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Selected leptonPt: {}", selected_lepton_pt);
        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Selected leptonIso: {}", selected_lepton_iso);

        // construct the four vectors of the selected leptons and taus to check
        // deltaR and reject a pair if the candidates are too close

        for (auto &candidate : sorted_pairs) {
            auto leptonindex = original_lepton_indices[candidate.first];
            ROOT::Math::PtEtaPhiMVector lepton = ROOT::Math::PtEtaPhiMVector(
                lepton_pt.at(leptonindex), lepton_eta.at(leptonindex),
                lepton_phi.at(leptonindex), lepton_mass.at(leptonindex));
            Logger::get("semileptonic::PairSelectionAlgo")
                ->debug("{} lepton vector: {}", leptonindex, lepton);
            auto tauindex = original_tau_indices[candidate.second];
            ROOT::Math::PtEtaPhiMVector tau = ROOT::Math::PtEtaPhiMVector(
                tau_pt.at(tauindex), tau_eta.at(tauindex), tau_phi.at(tauindex),
                tau_mass.at(tauindex));
            Logger::get("semileptonic::PairSelectionAlgo")
                ->debug("{} tau vector: {}", tauindex, tau);
            Logger::get("semileptonic::PairSelectionAlgo")
                ->debug("DeltaR: {}",
                        ROOT::Math::VectorUtil::DeltaR(lepton, tau));
            if (ROOT::Math::VectorUtil::DeltaR(lepton, tau) > mindeltaR) {
                Logger::get("semileptonic::PairSelectionAlgo")
                    ->debug(
                        "Selected original pair indices: mu = {} , tau = {}",
                        leptonindex, tauindex);
                Logger::get("semileptonic::PairSelectionAlgo")
                    ->debug("leptonPt = {} , TauPt = {} ", lepton.Pt(),
                            tau.Pt());
                Logger::get("semileptonic::PairSelectionAlgo")
                    ->debug("leptonPt = {} , TauPt = {} ",
                            lepton_pt[leptonindex], tau_pt[tauindex]);
                Logger::get("semileptonic::PairSelectionAlgo")
                    ->debug("leptonIso = {} , TauIso = {} ",
                            lepton_iso[leptonindex], tau_iso[tauindex]);
                selected_pair = {static_cast<int>(leptonindex),
                                 static_cast<int>(tauindex)};
                break;
            }
        }
        Logger::get("semileptonic::PairSelectionAlgo")
            ->debug("Final pair {} {}", selected_pair[0], selected_pair[1]);

        return selected_pair;
    };
}

} // end namespace semileptonic

namespace fullhadronic {

/// Implementation of the ditau pair selection algorithm for the fullhadronic
/// channel. First, only events that contain two goodTaus are considered. Events
/// contain at two good tau, if the tau_mask has at least two nonzero elemts.
/// These mask is contructed constructed using the functions from the
/// physicsobject namespace (e.g. physicsobject::CutPt).
///
/// \returns an `ROOT::RVec<int>` with two values, the first one beeing
/// the leading tau index and the second one beeing trailing tau index.
auto PairSelectionAlgo(const float &mindeltaR) {
    Logger::get("fullhadronic::PairSelectionAlgo")
        ->debug("Setting up algorithm");
    return [mindeltaR](const ROOT::RVec<float> &tau_pt,
                       const ROOT::RVec<float> &tau_eta,
                       const ROOT::RVec<float> &tau_phi,
                       const ROOT::RVec<float> &tau_mass,
                       const ROOT::RVec<float> &tau_iso,
                       const ROOT::RVec<int> &tau_mask) {
        // first entry is the leading tau index,
        // second entry is the trailing tau index
        ROOT::RVec<int> selected_pair = {-1, -1};
        const auto original_tau_indices = ROOT::VecOps::Nonzero(tau_mask);

        if (original_tau_indices.size() < 2) {
            return selected_pair;
        }
        Logger::get("fullhadronic::PairSelectionAlgo")
            ->debug("Running algorithm on good taus");

        const auto selected_tau_pt =
            ROOT::VecOps::Take(tau_pt, original_tau_indices);
        const auto selected_tau_iso =
            ROOT::VecOps::Take(tau_iso, original_tau_indices);

        Logger::get("fullhadronic::PairSelectionAlgo")
            ->debug("Original TauPt: {}", tau_pt);
        Logger::get("fullhadronic::PairSelectionAlgo")
            ->debug("Original TauIso: {}", tau_iso);

        Logger::get("fullhadronic::PairSelectionAlgo")
            ->debug("Selected TauPt: {}", selected_tau_pt);
        Logger::get("fullhadronic::PairSelectionAlgo")
            ->debug("Selected TauIso: {}", selected_tau_iso);

        const auto pair_indices = ROOT::VecOps::Combinations(
            selected_tau_pt, 2); // Gives indices of tau-tau pairs
        Logger::get("fullhadronic::PairSelectionAlgo")
            ->debug("Pairs: {} {}", pair_indices[0], pair_indices[1]);

        const auto pairs = ROOT::VecOps::Construct<std::pair<UInt_t, UInt_t>>(
            pair_indices[0], pair_indices[1]);
        Logger::get("fullhadronic::PairSelectionAlgo")
            ->debug("Pairs size: {}", pairs.size());
        int counter = 0;
        for (auto &pair : pairs) {
            counter++;
            Logger::get("fullhadronic::PairSelectionAlgo")
                ->debug("Constituents pair {}. : {} {}", counter, pair.first,
                        pair.second);
        }

        const auto sorted_pairs = ROOT::VecOps::Sort(
            pairs, compareForPairs(selected_tau_pt, selected_tau_iso,
                                   selected_tau_pt, selected_tau_iso));

        // construct the four vectors of the selected taus to check
        // deltaR and reject a pair if the candidates are too close
        bool found = false;
        for (auto &candidate : sorted_pairs) {
            auto tau_index_1 = original_tau_indices[candidate.first];
            ROOT::Math::PtEtaPhiMVector tau_1 = ROOT::Math::PtEtaPhiMVector(
                tau_pt.at(tau_index_1), tau_eta.at(tau_index_1),
                tau_phi.at(tau_index_1), tau_mass.at(tau_index_1));
            Logger::get("fullhadronic::PairSelectionAlgo")
                ->debug("{} leadint tau vector: {}", tau_index_1, tau_1);
            auto tau_index_2 = original_tau_indices[candidate.second];
            ROOT::Math::PtEtaPhiMVector tau_2 = ROOT::Math::PtEtaPhiMVector(
                tau_pt.at(tau_index_2), tau_eta.at(tau_index_2),
                tau_phi.at(tau_index_2), tau_mass.at(tau_index_2));
            Logger::get("fullhadronic::PairSelectionAlgo")
                ->debug("{} tau vector: {}", tau_index_2, tau_2);
            Logger::get("fullhadronic::PairSelectionAlgo")
                ->debug("DeltaR: {}",
                        ROOT::Math::VectorUtil::DeltaR(tau_1, tau_2));
            if (ROOT::Math::VectorUtil::DeltaR(tau_1, tau_2) > mindeltaR) {
                Logger::get("fullhadronic::PairSelectionAlgo")
                    ->debug("Selected original pair indices: tau_1 = {} , "
                            "tau_2 = {}",
                            tau_index_1, tau_index_2);
                Logger::get("fullhadronic::PairSelectionAlgo")
                    ->debug("Tau_1 Pt = {} , Tau_2 Pt = {} ", tau_1.Pt(),
                            tau_2.Pt());
                Logger::get("fullhadronic::PairSelectionAlgo")
                    ->debug("Tau_1 Iso = {} , Tau_2 Iso = {} ",
                            tau_iso[tau_index_1], tau_iso[tau_index_2]);
                selected_pair = {static_cast<int>(tau_index_1),
                                 static_cast<int>(tau_index_2)};
                found = true;
                break;
            }
        }
        // sort it that the leading tau in pt is first
        if (found) {
            if (tau_pt.at(selected_pair[0]) < tau_pt.at(selected_pair[1])) {
                std::swap(selected_pair[0], selected_pair[1]);
            }
        }
        Logger::get("fullhadronic::PairSelectionAlgo")
            ->debug("Final pair {} {}", selected_pair[0], selected_pair[1]);

        return selected_pair;
    };
}

} // end namespace fullhadronic

namespace mutau {

/**
 * @brief Function used to select the pair of tau leptons based on the standard
 * pair selection algorithm
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the muTau pair selection these values are:
    - tau_pt
    - tau_eta
    - tau_phi
    - tau_mass
    - tau_iso
    - muon_pt
    - muon_eta
    - muon_phi
    - muon_mass
    - muon_iso
    - tau_mask containing the flags whether the tau is a good tau or not
    - muon_mask containing the flags whether the muon is a good muon or not
 * @param pairname name of the new column containing the pair index
 * @param mindeltaR the seperation between the muon at the tau has to be larger
 than
 * this value
 * @return a new dataframe with the pair index column added
 */
auto PairSelection(auto &df, const std::vector<std::string> &input_vector,
                   const std::string &pairname, const float &mindeltaR) {
    Logger::get("mutau::PairSelection")
        ->debug("Setting up MuTau pair building");
    auto df1 = df.Define(
        pairname, pairselection::semileptonic::PairSelectionAlgo(mindeltaR),
        input_vector);
    return df1;
}

} // end namespace mutau

namespace eltau {

/**
 * @brief Function used to select the pair of tau leptons based on the standard
 * pair selection algorithm
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the ElTau pair selection these values are:
    - tau_pt
    - tau_eta
    - tau_phi
    - tau_mass
    - tau_iso
    - electron_pt
    - electron_eta
    - electron_phi
    - electron_mass
    - electron_iso
    - tau_mask containing the flags whether the tau is a good tau or not
    - electron_mask containing the flags whether the electron is a good electron
 or not
 * @param pairname name of the new column containing the pair index
 * @param mindeltaR the seperation between the electron at the tau has to be
 larger than
 * this value
 * @return a new dataframe with the pair index column added
 */
auto PairSelection(auto &df, const std::vector<std::string> &input_vector,
                   const std::string &pairname, const float &mindeltaR) {
    Logger::get("eltau::PairSelection")
        ->debug("Setting up ElTau pair building");
    auto df1 = df.Define(
        pairname, pairselection::semileptonic::PairSelectionAlgo(mindeltaR),
        input_vector);
    return df1;
}

} // end namespace eltau

namespace tautau {

/**
 * @brief Function used to select the pair of tau leptons based on the standard
 * pair selection algorithm
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the TauTau pair selection these values are:
    - tau_pt
    - tau_eta
    - tau_phi
    - tau_mass
    - tau_iso
    - tau_mask containing the flags whether the tau is a good tau or not
 * @param pairname name of the new column containing the pair index
 * @param mindeltaR the seperation between the two tau candidates has to be
 larger than
 * this value
 * @return a new dataframe with the pair index column added
 */
auto PairSelection(auto &df, const std::vector<std::string> &input_vector,
                   const std::string &pairname, const float &mindeltaR) {
    Logger::get("tautau::PairSelection")
        ->debug("Setting up TauTau pair building");
    auto df1 = df.Define(
        pairname, pairselection::fullhadronic::PairSelectionAlgo(mindeltaR),
        input_vector);
    return df1;
}

} // namespace tautau

namespace elmu {

/**
 * @brief Function used to select the pair of tau leptons based on the standard
 * pair selection algorithm
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the ElTau pair selection these values are:
    - electron_pt
    - electron_eta
    - electron_phi
    - electron_mass
    - electron_iso
    - muon_pt
    - muon_eta
    - muon_phi
    - muon_mass
    - muon_iso
    - electron_mask containing the flags whether the electron is a good electron
 or not
    - muon_mask containing the flags whether the muon is a good muon or not
 * @param pairname name of the new column containing the pair index
 * @param mindeltaR the seperation between the electron and the muon has to be
 larger than
 * this value
 * @return a new dataframe with the pair index column added
 */
auto PairSelection(auto &df, const std::vector<std::string> &input_vector,
                   const std::string &pairname, const float &mindeltaR) {
    Logger::get("elmu::PairSelection")->debug("Setting up elmu pair building");
    // TODO Implement the fully leptonic tau pair selection algorithm
    // auto df1 =
    //     df.Define(pairname,
    //     pairselection::semileptonic::PairSelectionAlgo(mindeltaR),
    //               input_vector);
    return df;
}

} // namespace elmu

namespace mumu {
/**
 * @brief Lambda function containg the algorithm to select the pair of muons
 * with the highest pt
 *
 * @param mindeltaR the seperation between the two muons has to be larger than
 this value
 *
 * @return vector with two entries, the first entry is the leading muon index,
    the second entry is the trailing muon index
 */
auto PairSelectionAlgo(const float &mindeltaR) {
    Logger::get("PairSelection")->debug("Setting up algorithm");
    return [mindeltaR](const ROOT::RVec<float> &muon_pt,
                       const ROOT::RVec<float> &muon_eta,
                       const ROOT::RVec<float> &muon_phi,
                       const ROOT::RVec<float> &muon_mass,
                       const ROOT::RVec<int> &muon_mask) {
        // first entry is the leading muon index,
        // second entry is the trailing muon index
        ROOT::RVec<int> selected_pair = {-1, -1};
        const auto original_muon_indices = ROOT::VecOps::Nonzero(muon_mask);
        // we need at least two fitting muons
        if (original_muon_indices.size() < 2) {
            return selected_pair;
        }
        const auto good_pts =
            ROOT::VecOps::Take(muon_pt, original_muon_indices);
        const auto good_etas =
            ROOT::VecOps::Take(muon_eta, original_muon_indices);
        const auto good_phis =
            ROOT::VecOps::Take(muon_phi, original_muon_indices);
        const auto good_masses =
            ROOT::VecOps::Take(muon_mass, original_muon_indices);
        auto fourVecs = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(
            good_pts, good_etas, good_phis, good_masses);
        auto selected_mu_indices = std::vector<int>{-1, -1};
        auto selected_pts = std::vector<float>{-1, -1};
        auto combinations =
            ROOT::VecOps::Combinations(original_muon_indices, 2);
        if (original_muon_indices.size() > 2) {
            Logger::get("MMPairSelectionAlgo")
                ->debug("More than two suitable muons found, printing "
                        "combinations.... ");
            for (auto &comb : combinations) {
                Logger::get("MMPairSelectionAlgo")->debug("index: {}", comb);
            };
            Logger::get("MMPairSelectionAlgo")->debug("---------------------");
        }
        for (int n = 0; n < combinations[0].size(); n++) {
            auto mu1 = fourVecs[combinations[0][n]];
            auto mu2 = fourVecs[combinations[1][n]];
            auto deltaR = ROOT::Math::VectorUtil::DeltaR(mu1, mu2);
            Logger::get("MMPairSelectionAlgo")
                ->debug("deltaR check: {}", deltaR);
            if (deltaR > mindeltaR) {
                if (mu1.Pt() >= selected_pts[0] &&
                    mu2.Pt() >= selected_pts[1]) {
                    selected_pts[0] = mu1.Pt();
                    selected_pts[1] = mu2.Pt();
                    selected_mu_indices[0] = combinations[0][n];
                    selected_mu_indices[1] = combinations[1][n];
                }
            }
        }

        if (good_pts[selected_mu_indices[0]] <
            good_pts[selected_mu_indices[1]]) {
            std::swap(selected_mu_indices[0], selected_mu_indices[1]);
        }
        Logger::get("MMPairSelectionAlgo")->debug("good pts: {}", good_pts);
        Logger::get("MMPairSelectionAlgo")
            ->debug("selected_mu_indices: {}, {}", selected_mu_indices[0],
                    selected_mu_indices[1]);
        selected_pair = {static_cast<int>(selected_mu_indices[0]),
                         static_cast<int>(selected_mu_indices[1])};
        return selected_pair;
    };
}
/**
 * @brief Lambda function containg the algorithm to select the pair
 * of muons closest to the Z mass
 *
 * @param mindeltaR the seperation between the two muons has to be larger than
 * this value
 * @return vector with two entries, the first entry is the leading muon index,
 *  the second entry is the trailing muon index
 */
auto ZBosonPairSelectionAlgo(const float &mindeltaR) {
    Logger::get("PairSelection")->debug("Setting up algorithm");
    return [mindeltaR](const ROOT::RVec<float> &muon_pt,
                       const ROOT::RVec<float> &muon_eta,
                       const ROOT::RVec<float> &muon_phi,
                       const ROOT::RVec<float> &muon_mass,
                       const ROOT::RVec<int> &muon_mask) {
        // first entry is the leading muon index,
        // second entry is the trailing muon index
        ROOT::RVec<int> selected_pair = {-1, -1};
        const auto original_muon_indices = ROOT::VecOps::Nonzero(muon_mask);
        // we need at least two fitting muons
        if (original_muon_indices.size() < 2) {
            return selected_pair;
        }
        Logger::get("ZBosonPairSelectionAlgo")
            ->debug("Running algorithm on good muons");

        const auto good_pts =
            ROOT::VecOps::Take(muon_pt, original_muon_indices);
        const auto good_etas =
            ROOT::VecOps::Take(muon_eta, original_muon_indices);
        const auto good_phis =
            ROOT::VecOps::Take(muon_phi, original_muon_indices);
        const auto good_masses =
            ROOT::VecOps::Take(muon_mass, original_muon_indices);
        auto fourVecs = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(
            good_pts, good_etas, good_phis, good_masses);
        float mass_difference = -1.0;
        float zmass_candidate = -1.0;
        auto selected_mu_indices = std::vector<int>{-1, -1};
        if (original_muon_indices.size() > 2) {
            Logger::get("ZBosonPairSelectionAlgo")
                ->debug("More than two potential muons found. running "
                        "algorithm to find Z Boson muon pairs");
            Logger::get("ZBosonPairSelectionAlgo")
                ->debug("original_muon_indices: {}", original_muon_indices);
            for (auto &fourVec : fourVecs) {
                Logger::get("ZBosonPairSelectionAlgo")
                    ->debug("fourVec: {}", fourVec);
            }
        }
        auto combinations =
            ROOT::VecOps::Combinations(original_muon_indices, 2);
        Logger::get("ZBosonPairSelectionAlgo")
            ->debug("printing combinations.... ");
        for (auto &comb : combinations) {
            Logger::get("ZBosonPairSelectionAlgo")->debug("index: {}", comb);
        };
        Logger::get("ZBosonPairSelectionAlgo")->debug("---------------------");

        for (int n = 0; n < combinations[0].size(); n++) {
            auto mu1 = fourVecs[combinations[0][n]];
            auto mu2 = fourVecs[combinations[1][n]];
            auto deltaR = ROOT::Math::VectorUtil::DeltaR(mu1, mu2);
            zmass_candidate = (mu1 + mu2).M();
            Logger::get("ZBosonPairSelectionAlgo")
                ->debug("deltaR check: {}", deltaR);
            Logger::get("ZBosonPairSelectionAlgo")
                ->debug("mass check: {}", zmass_candidate);
            if (deltaR > mindeltaR) {
                if (std::abs(91.2 - zmass_candidate) < mass_difference ||
                    mass_difference < 0) {
                    mass_difference = std::abs(91.2 - zmass_candidate);
                    selected_mu_indices[0] = combinations[0][n];
                    selected_mu_indices[1] = combinations[1][n];
                }
            }
        }
        if (good_pts[selected_mu_indices[0]] <
            good_pts[selected_mu_indices[1]]) {
            std::swap(selected_mu_indices[0], selected_mu_indices[1]);
        }
        Logger::get("ZBosonPairSelectionAlgo")->debug("good pts: {}", good_pts);
        Logger::get("ZBosonPairSelectionAlgo")
            ->debug("selected_mu_indices: {}, {}", selected_mu_indices[0],
                    selected_mu_indices[1]);

        selected_pair = {static_cast<int>(selected_mu_indices[0]),
                         static_cast<int>(selected_mu_indices[1])};
        return selected_pair;
    };
}
/**
 * @brief Function used to select the pair of muons with the highest
 * pt
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns
 * needed for the alogrithm. For the muon pair selection the required paramters
 are:
    - muon_pt
    - muon_eta
    - muon_phi
    - muon_mass
    - muon_mask containing the flags whether the muon is a good muon or not
 * @param pairname name of the new column containing the pair index
 * @param mindeltaR the seperation between the two muons has to be larger than
 * this value
 * @return a new dataframe with the pair index column added
 */
auto PairSelection(auto &df, const std::vector<std::string> &input_vector,
                   const std::string &pairname, const float &mindeltaR) {
    Logger::get("MuonPairSelection")->debug("Setting up mumu pair building");
    auto df1 =
        df.Define(pairname, pairselection::mumu::PairSelectionAlgo(mindeltaR),
                  input_vector);
    return df1;
}
/**
 * @brief Function used to select the pair of muons closest to the Z
 * mass
 *
 * @param df the input dataframe
 * @param input_vector . For the Z boson muon pair selection the required
 paramters are:
    - muon_pt
    - muon_eta
    - muon_phi
    - muon_mass
    - muon_mask containing the flags whether the muon is a good muon or not
 * @param pairname name of the new column containing the pair index
 * @param mindeltaR the seperation between the two muons has to be larger than
 * this value
 * @return a new dataframe with the pair index column added
 */
auto ZBosonPairSelection(auto &df, const std::vector<std::string> &input_vector,
                         const std::string &pairname, const float &mindeltaR) {
    Logger::get("ZBosonPairSelection")
        ->debug("Setting up Z boson mumu pair building");
    auto df1 = df.Define(
        pairname, pairselection::mumu::ZBosonPairSelectionAlgo(mindeltaR),
        input_vector);
    return df1;
}

} // end namespace mumu
} // end namespace pairselection
