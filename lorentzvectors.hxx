#include <Math/Vector4D.h>
#include "utility/Logger.hxx"

namespace lorentzvectors{

    auto buildparticle(auto df, const std::string& particle, const std::string& pairname, const int& position){
        std::vector<std::string> cols;
        cols.push_back(pairname);
        cols.push_back(particle + "_pt");
        cols.push_back(particle + "_eta");
        cols.push_back(particle + "_phi");
        cols.push_back(particle + "_mass");
        const auto outputcol = "particle_" + std::to_string(position+1) + "_p4";
        auto df1 = df.Define(outputcol,
            [position](const ROOT::RVec<int>& pair, const ROOT::RVec<float>& pts, const ROOT::RVec<float>& etas, const ROOT::RVec<float>& phis, const ROOT::RVec<float>& masses){
                // the index of the particle is stored in the pair vector
                const int index = pair[position];
                // auto p4s = ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(pts, etas, phis, masses);
                //return p4s[index];
                Logger::get("lorentzvectors")->debug("pair {}", pair);
                Logger::get("lorentzvectors")->debug("pts {}", pts);
                Logger::get("lorentzvectors")->debug("etas {}", etas);
                Logger::get("lorentzvectors")->debug("phis {}", phis);
                Logger::get("lorentzvectors")->debug("masses {}", masses);
                Logger::get("lorentzvectors")->debug("Index {}", index);

                const auto p4 = ROOT::Math::PtEtaPhiMVector(pts[index], etas[index], phis[index], masses[index]);
                Logger::get("lorentzvectors")->debug("P4 - Particle {} : {}", position, p4);
                return p4;
                }, cols);
        return df1;
    }

    namespace mutau{

        auto build(auto df, const std::string pairname)
        {
            // first the muon
            auto df1 = lorentzvectors::buildparticle(df, "Muon", pairname, 0);
            // now the tau
            auto df2 = lorentzvectors::buildparticle(df1, "Tau", pairname, 1);
            return df2;
        }
    }
}