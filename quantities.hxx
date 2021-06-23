#include "defaults.hxx"
#include "utility/Logger.hxx"
#include <Math/Vector4D.h>
/// The namespace that is used to hold the functions for basic quantities that
/// are needed for every event
namespace quantities {

/// Function to add an input quantity under a different name
///
/// \param df the dataframe to add the quantity to
/// \param inputname name of the existing column
/// \param outputname name of the new column
///
/// \returns a dataframe with the new column

template <typename T>
auto rename(auto &df, const std::string &inputname,
            const std::string &outputname) {
    return df.Define(outputname, [](const T &q) { return q; }, {inputname});
}
/// Function to calculate the pt from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

auto pt(auto &df, const std::string &outputname,
        const std::string &inputvector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.pt(); },
        {inputvector});
}
/// Function to calculate the eta from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the eta value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

auto eta(auto &df, const std::string &outputname,
         const std::string &inputvector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.eta(); },
        {inputvector});
}
/// Function to calculate the eta from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the eta value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

auto phi(auto &df, const std::string &outputname,
         const std::string &inputvector) {
    return df.Define(outputname,
                     [](const ROOT::Math::PtEtaPhiMVector &p4) {
                         if (p4.pt() <
                             0.0) // negative pt is used to mark invalid LVs
                             return default_float;
                         return (float)p4.phi();
                     },
                     {inputvector});
}
/// Function to calculate the mass from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the mass value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

auto mass(auto &df, const std::string &outputname,
          const std::string &inputvector) {
    return df.Define(outputname,
                     [](const ROOT::Math::PtEtaPhiMVector &p4) {
                         if (p4.pt() <
                             0.0) // negative pt is used to mark invalid LVs
                             return default_float;
                         return (float)p4.mass();
                     },
                     {inputvector});
}
/// Function to writeout the dxy impact parameter from a particle. The particle
/// is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the dxy value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param dxycolumn name of the column containing the dxy values
///
/// \returns a dataframe with the new column

auto dxy(auto &df, const std::string &outputname, const int &position,
         const std::string &pairname, const std::string &dxycolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<float> &dxy) {
            const int index = pair.at(position);
            return dxy.at(index, default_float);
        },
        {pairname, dxycolumn});
}
/// Function to writeout the dz impact parameter from a particle. The particle
/// is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the dz value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param dzcolumn name of the column containing the dz values
///
/// \returns a dataframe with the new column

auto dz(auto &df, const std::string &outputname, const int &position,
        const std::string &pairname, const std::string &dzcolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<float> &dz) {
            const int index = pair.at(position);
            return dz.at(index, default_float);
        },
        {pairname, dzcolumn});
}
/// Function to writeout the charge of a particle. The particle is identified
/// via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the charge value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param chargecolumn name of the column containing the charge values
///
/// \returns a dataframe with the new column

auto charge(auto &df, const std::string &outputname, const int &position,
            const std::string &pairname, const std::string &chargecolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<int> &charge) {
            const int index = pair.at(position);
            return charge.at(index, default_int);
        },
        {pairname, chargecolumn});
}
/// Function to calculate the visible mass from a pair of lorentz vectors and
/// add it to the dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column

auto m_vis(auto &df, const std::string &outputname,
           const std::vector<std::string> &inputvectors) {
    // build visible mass from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
                return default_float;
            auto const dileptonsystem = p4_1 + p4_2;
            return (float)dileptonsystem.mass();
        },
        inputvectors);
}
/// Function to writeout the isolation of a particle. The particle is identified
/// via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the isolation value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param isolationcolumn name of the column containing the isolation values
///
/// \returns a dataframe with the new column

auto isolation(auto &df, const std::string &outputname, const int &position,
               const std::string &pairname,
               const std::string &isolationcolumn) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<float> &isolation) {
                         const int index = pair.at(position);
                         return isolation.at(index, default_float);
                     },
                     {pairname, isolationcolumn});
}
/// Function to writeout the PDGID from a genparticle. The particle
/// is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the PDGID
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param pdgidcolumn name of the column containing the pdgID
///
/// \returns a dataframe with the new column

auto pdgid(auto &df, const std::string &outputname, const int &position,
           const std::string &pairname, const std::string &pdgidcolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<int> &pdgid) {
            const int index = pair.at(position);
            return pdgid.at(index, default_pdgid);
        },
        {pairname, pdgidcolumn});
}
/// namespace for tau specific quantities
namespace tau {
/// Function to writeout the decaymode of a tau. The particle is identified via
/// the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the decaymode value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param decaymodecolumn name of the column containing the decaymode values
///
/// \returns a dataframe with the new column

auto decaymode(auto &df, const std::string &outputname, const int &position,
               const std::string &pairname,
               const std::string &decaymodecolumn) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &decaymode) {
                         const int index = pair.at(position);
                         return decaymode.at(index, default_int);
                     },
                     {pairname, decaymodecolumn});
}
/// Function to writeout the genmatch of a tau. The particle is identified via
/// the index stored in the pair vector Genmatch values are defined as \code
///   1 = prompt electron,
///   2 = prompt muon,
///   3 = tau->e decay,
///   4 = tau->mu decay,
///   5 = hadronic tau decay,
///   0 = unknown or unmatched
///   \endcode
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the genmatch value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param genmatchcolumn name of the column containing the genmatch values
///
/// \returns a dataframe with the new column

auto genmatch(auto &df, const std::string &outputname, const int &position,
              const std::string &pairname, const std::string &genmatchcolumn) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<UChar_t> &genmatch) {
                         const int index = pair.at(position);
                         return genmatch.at(index, default_uchar);
                     },
                     {pairname, genmatchcolumn});
}
/// Function to writeout the pt of the reco jet associated with the given tau.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the jet pt value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param taujet_index name of the column containing the association between
/// the tau and the reco jet \param jetpt_column name of the column containing
/// the recojet pt values
///
/// \returns a dataframe with the new column
auto matching_jet_pt(auto &df, const std::string &outputname,
                     const int &position, const std::string &pairname,
                     const std::string &taujet_index,
                     const std::string &jetpt_column) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &taujets,
                                const ROOT::RVec<float> &jetpt) {
                         const int tauindex = pair.at(position);
                         const int jetindex = taujets.at(tauindex, -1);
                         return jetpt.at(jetindex, default_float);
                     },
                     {pairname, taujet_index, jetpt_column});
}
/// Function to writeout the pt of the gen jet associated with the reco jet,
/// which is associated with the given tau. \code
///  Tau --> recoJet --> GenJet
///   \endcode
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the jet pt value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param taujet_index name of the column containing the association between
/// the tau and the reco jet \param genjet_index name of the column containing
/// the association between the reco jet and the gen jet \param genjetpt_column
/// name of the column containing the genJet pt values
///
/// \returns a dataframe with the new column
auto matching_genjet_pt(auto &df, const std::string &outputname,
                        const int &position, const std::string &pairname,
                        const std::string &taujet_index,
                        const std::string &genjet_index,
                        const std::string &genjetpt_column) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &taujets,
                                const ROOT::RVec<int> &genjets,
                                const ROOT::RVec<float> &genjetpt) {
                         const int tauindex = pair.at(position);
                         const int jetindex = taujets.at(tauindex, -1);
                         const int genjetindex = genjets.at(jetindex, -1);
                         return genjetpt.at(genjetindex, default_int);
                     },
                     {pairname, taujet_index, genjet_index, genjetpt_column});
}
} // end namespace tau
} // end namespace quantities