#include "utility/Logger.hxx"
#include <Math/Vector4D.h>

namespace lorentzvectors {

auto buildparticle(auto df, const std::vector<std::string> quantities,
                   const std::string outputname, const int &position) {
  // std::vector<std::string> cols;
  // cols.push_back(particle + "_pt");
  // cols.push_back(particle + "_eta");
  // cols.push_back(particle + "_phi");
  // cols.push_back(particle + "_mass");
  // const auto outputcol = "p4_" + std::to_string(position+1);
  auto df1 = df.Define(
      outputname,
      [position](const ROOT::RVec<int> &pair, const ROOT::RVec<float> &pts,
                 const ROOT::RVec<float> &etas, const ROOT::RVec<float> &phis,
                 const ROOT::RVec<float> &masses) {
        // the index of the particle is stored in the pair vector
        const int index = pair[position];
        Logger::get("lorentzvectors")->debug("pair {}", pair);
        Logger::get("lorentzvectors")->debug("pts {}", pts);
        Logger::get("lorentzvectors")->debug("etas {}", etas);
        Logger::get("lorentzvectors")->debug("phis {}", phis);
        Logger::get("lorentzvectors")->debug("masses {}", masses);
        Logger::get("lorentzvectors")->debug("Index {}", index);

        const auto p4 = ROOT::Math::PtEtaPhiMVector(pts[index], etas[index],
                                                    phis[index], masses[index]);
        Logger::get("lorentzvectors")
            ->debug("P4 - Particle {} : {}", position, p4);
        return p4;
      },
      quantities);
  return df1;
}

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