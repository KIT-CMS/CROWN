#ifndef GUARD_TOPRECO_H
#define GUARD_TOPRECO_H

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

ROOT::RDF::RNode LeptonSelection(ROOT::RDF::RNode df,
				 const std::string &str_n_loose_mu,
				 const std::string &str_n_loose_el,
				 const std::string &str_n_tight_mu,
				 const std::string &str_n_tight_el,
				 const std::string &tight_muons_mask,
				 const std::string &tight_electrons_mask,
				 const std::string &str_mu_pt,
				 const std::string &str_mu_eta,
				 const std::string &str_mu_phi,
				 const std::string &str_mu_mass,
				 const std::string &str_el_pt,
				 const std::string &str_el_eta,
				 const std::string &str_el_phi,
				 const std::string &str_el_mass,
				 const std::string &str_n_loose_lep,
				 const std::string &str_n_tight_lep,
				 const std::string &str_is_mu,
				 const std::string &str_is_el,
				 const std::string &str_lep_p4
				 );

double rad_py(double x, double lep_px);
double min_fplus(double *par);
double min_fminus(double *par);
void fcn_plus(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcn_minus(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

ROOT::RDF::RNode ReconstructLeptonicW(ROOT::RDF::RNode df,
				      const std::string &str_lep_p4,
				      const std::string &str_met_p4,
				      const std::string &str_lepw_p4
				      );

ROOT::RDF::RNode ReconstructLeptonicW_mt(ROOT::RDF::RNode df,
					 const std::string &outputname,
					 const std::string &particle_p4
					 );


#endif /* GUARD_TOPRECO_H */
