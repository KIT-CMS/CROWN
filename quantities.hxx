#include <Math/Vector4D.h>
#include "utility/Logger.hxx"

namespace quantities {
    auto pt(auto df){
        // for particle 1
        auto df1 = df.Define("pt_1",[](const ROOT::Math::PtEtaPhiMVector& p4){ return p4.pt();},{"p4_1"});
        // for particle 2
        return df1.Define("pt_2",[](const ROOT::Math::PtEtaPhiMVector& p4){ return p4.pt();},{"p4_2"});
    }
    auto eta(auto df){
        // for particle 1
        auto df1 = df.Define("eta_1",[](const ROOT::Math::PtEtaPhiMVector& p4){ return p4.eta();},{"p4_1"});
        // for particle 2
        return df1.Define("eta_2",[](const ROOT::Math::PtEtaPhiMVector& p4){ return p4.eta();},{"p4_2"});
    }
    auto phi(auto df){
        // for particle 1
        auto df1 = df.Define("phi_1",[](const ROOT::Math::PtEtaPhiMVector& p4){ return p4.phi();},{"p4_1"});
        // for particle 2
        return df1.Define("phi_2",[](const ROOT::Math::PtEtaPhiMVector& p4){ return p4.phi();},{"p4_2"});
    }
    auto mvis(auto df){
        // build visible mass from the two particles
        return df.Define("m_vis",[](const ROOT::Math::PtEtaPhiMVector& p4_1, const ROOT::Math::PtEtaPhiMVector& p4_2){
            auto const dileptonsystem = p4_1 + p4_2;
            return dileptonsystem.mass();
        },{"p4_1", "p4_2"});
    }
} // end namespace quantities