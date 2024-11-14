#ifndef GUARD_EXAMPLE_H
#define GUARD_EXAMPLE_H

namespace example {
ROOT::RDF::RNode CutNMuon(ROOT::RDF::RNode df, const std::string &nMuons, const unsigned int &nMuons_req);

ROOT::RDF::RNode CutCMuons(ROOT::RDF::RNode df, const std::string &CMuons, const int &CMuons_req);

ROOT::RDF::RNode GetInvariantMass(
    ROOT::RDF::RNode df, 
    const std::string &new_column, 
    const std::string &particle_pts, 
    const std::string &particle_etas, 
    const std::string &particle_phis, 
    const std::string &particle_masses
);

ROOT::RDF::RNode MuonCSum(ROOT::RDF::RNode df, const std::string &new_column, const std::string &CMuons);

} // namespace example
#endif /* GUARD_EXAMPLE_H */