#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "utility.hxx"

namespace pairselection {

    auto compareQuantities = [](auto value1, auto value2){
        if(not utility::ApproxEqual(value1, value2))
        {
            return value1 > value2;
        }
        else
        {
            return false;
        }
    };

    namespace mutau {
        // 1. Muon Isolation
        // 2. Muon pt
        // 3. Tau Isolation
        // 4. Tau pt

        auto PairSelectionAlgo(){
            std::cout << "Setting up algorithm \n";
            return [](const ROOT::RVec<float>& taupt, const ROOT::RVec<float>& tauiso, const ROOT::RVec<float>& muonpt,
                      const ROOT::RVec<float>& muoniso, const ROOT::RVec<int>& taumask, const ROOT::RVec<int>& muonmask){
                        std::cout << "Running algorithm \n";
                        ROOT::RVec<int> selected_pair; // first entry is the muon index, second entry is the tau index
                        auto selected_tau_indices = ROOT::VecOps::Nonzero(taumask);
                        auto selected_muon_indices = ROOT::VecOps::Nonzero(muonmask);

                        auto selected_taupt = ROOT::VecOps::Take(taupt, selected_tau_indices);
                        auto selected_tauiso = ROOT::VecOps::Take(tauiso, selected_tau_indices);
                        auto selected_muonpt = ROOT::VecOps::Take(muonpt, selected_muon_indices);
                        auto selected_muoniso = ROOT::VecOps::Take(muoniso, selected_muon_indices);

                        auto sorted_taupt = ROOT::VecOps::Argsort(selected_taupt, pairselection::compareQuantities);
                        auto sorted_tauiso = ROOT::VecOps::Argsort(selected_tauiso, pairselection::compareQuantities); // tau discriminator: best at large values, so opposite to isolation
                        auto sorted_muonpt = ROOT::VecOps::Argsort(selected_muonpt, pairselection::compareQuantities);
                        auto sorted_muoniso = ROOT::VecOps::Argsort(-1. * selected_muoniso, pairselection::compareQuantities);
                        std::cout << "TauPt ranking: " << sorted_taupt << std::endl;
                        std::cout << "TauPt: " << selected_taupt << std::endl;
                        std::cout << "TauIso ranking: " << sorted_tauiso << std::endl;
                        std::cout << "TauIso: " << selected_tauiso << std::endl;
                        std::cout << "MuonPt ranking: " << sorted_muonpt << std::endl;
                        std::cout << "MuonPt: " << selected_muonpt << std::endl;
                        std::cout << "MuonIso ranking: " << sorted_muoniso << std::endl;
                        std::cout << "MuonIso: " << selected_muoniso << std::endl;
                selected_pair.push_back(10);
                return selected_pair;
            };
        }

        auto PairSelection(auto df, const std::string taumask, const std::string muonmask, const std::string pairname, const std::vector<std::string> pairvariables){
            std::cout << "Starting to build pair !\n";
            auto df1 = df.Define(pairname, pairselection::mutau::PairSelectionAlgo(), {"Tau_pt", "Tau_rawDeepTau2017v2p1VSjet", "Muon_pt", "Muon_pfRelIso04_all",taumask, muonmask});
            return df1;
        }




    } // end namespace mutau

} // end namespace pairselection