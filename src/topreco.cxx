#ifndef GUARD_TOPRECO_H
#define GUARD_TOPRECO_H

#include "../include/topreco.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include "TMinuit.h"
#include <Math/Boost.h>
#include "correction.h"

const float W_MASS = 80.377; // PDG value as of 10/22
const float TOP_MASS = 172.5; // gen mass

namespace topreco {

ROOT::RDF::RNode LeptonSelection(ROOT::RDF::RNode df,
				 const std::string &str_n_loose_mu,
				 const std::string &str_n_loose_el,
				 const std::string &str_n_tight_mu,
				 const std::string &str_n_tight_el,
				 const std::string &str_tight_muons_mask,
				 const std::string &str_tight_electrons_mask,
				 const std::string &str_n_antitight_mu,
				 const std::string &str_n_antitight_el,
				 const std::string &str_antitight_muons_mask,
				 const std::string &str_antitight_electrons_mask,
				 const std::string &str_mu_pt,
				 const std::string &str_mu_eta,
				 const std::string &str_mu_phi,
				 const std::string &str_mu_mass,
				 const std::string &str_mu_charge,
				 const std::string &str_el_pt,
				 const std::string &str_el_eta,
				 const std::string &str_el_detasc,
				 const std::string &str_el_phi,
				 const std::string &str_el_mass,
				 const std::string &str_el_charge,
				 const std::string &str_n_loose_lep,
				 const std::string &str_n_tight_lep,
				 const std::string &str_n_antitight_lep,
				 const std::string &str_is_mu,
				 const std::string &str_is_el,
				 const std::string &str_is_iso,
				 const std::string &str_lep_p4,
				 const std::string &str_lep_sceta,
				 const std::string &str_lep_charge
				 ) {


  auto df1 = df.Define(str_n_loose_lep,
		       [](const int &n_loose_mu,
			  const int &n_loose_el) {
			 Logger::get("lepsel")->debug("size of n_loose_mu and n_loose_el: {} {}",
						      n_loose_mu, n_loose_el);
			 return n_loose_mu + n_loose_el;
		       },
		       {str_n_loose_mu, str_n_loose_el}
		       );

  auto df2 = df1.Define(str_n_tight_lep,
			[](const int &n_tight_mu,
			   const int &n_tight_el) {
			  Logger::get("lepsel")->debug("size of n_tight_mu and n_tight_el: {} {}",
						       n_tight_mu, n_tight_el);
			  return n_tight_mu + n_tight_el;
			},
			{str_n_tight_mu, str_n_tight_el}
			);

  auto df3 = df2.Define(str_n_antitight_lep,
			[](const int &n_antitight_mu,
			   const int &n_antitight_el) {
			  Logger::get("lepsel")->debug("size of n_antitight_mu and n_antitight_el: {} {}",
						       n_antitight_mu, n_antitight_el);
			  return n_antitight_mu + n_antitight_el;
			},
			{str_n_antitight_mu, str_n_antitight_el}
			);

  auto df4 = df3.Define(str_is_mu,
			[](const int &n_tight_mu,
			   const int &n_antitight_mu,
			   const int &n_loose_lep,
			   const int &n_tight_lep,
			   const int &n_antitight_lep) {
			  return int(((n_tight_mu == 1) && (n_loose_lep == 1) && (n_tight_lep == 1)) || ((n_antitight_mu == 1) && (n_loose_lep == 0) && (n_antitight_lep == 1)));
			},
			{str_n_tight_mu, str_n_antitight_mu, str_n_loose_lep, str_n_tight_lep, str_n_antitight_lep}
			);

  auto df5 = df4.Define(str_is_el,
			[](const int &n_tight_el,
			   const int &n_antitight_el,
			   const int &n_loose_lep,
			   const int &n_tight_lep,
			   const int &n_antitight_lep) {
			  return int(((n_tight_el == 1) && (n_loose_lep == 1) && (n_tight_lep == 1)) || ((n_antitight_el == 1 )&& (n_loose_lep == 0) && (n_antitight_lep == 1)));
			},
			{str_n_tight_el, str_n_antitight_el, str_n_loose_lep, str_n_tight_lep, str_n_antitight_lep}
			);


  auto is_iso = [](const int &is_mu,
		   const int &is_el,
		   const int &n_tight_lep,
		   const int &n_antitight_lep){

    if (is_mu || is_el) {
      if (n_tight_lep == 1)
	return +1;
      if (n_antitight_lep == 1)
	return -1;
    }

    return 0;
  };

  auto df6 = df5.Define(str_is_iso,
			is_iso,
			{str_is_mu, str_is_el,
			    str_n_tight_lep, str_n_antitight_lep}
			);


  auto lep_p4 = [](const int is_mu,
		   const int is_el,
		   const int is_iso,
		   const ROOT::RVec<int> &tight_muons_mask,
		   const ROOT::RVec<int> &tight_electrons_mask,
		   const ROOT::RVec<int> &antitight_muons_mask,
		   const ROOT::RVec<int> &antitight_electrons_mask,
		   const ROOT::RVec<float> &mu_pt,
		   const ROOT::RVec<float> &mu_eta,
		   const ROOT::RVec<float> &mu_phi,
		   const ROOT::RVec<float> &mu_mass,
		   const ROOT::RVec<float> &el_pt,
		   const ROOT::RVec<float> &el_eta,
		   const ROOT::RVec<float> &el_phi,
		   const ROOT::RVec<float> &el_mass) {

    Logger::get("lep_p4")->debug("masks mu {}, el {}",
    				 tight_muons_mask, tight_electrons_mask);
    Logger::get("lep_p4")->debug("mask sizes mu {}, el {}",
    				 tight_muons_mask.size(), tight_electrons_mask.size());
    Logger::get("lep_p4")->debug("max in mu {}, el {}",
    				 ROOT::VecOps::Max(tight_muons_mask), ROOT::VecOps::Max(tight_electrons_mask));
    Logger::get("lep_p4")->debug("index of max in mu {}, el {}",
    				 ROOT::VecOps::ArgMax(tight_muons_mask), ROOT::VecOps::ArgMax(tight_electrons_mask));

    auto lep = ROOT::Math::PtEtaPhiMVector(-10,-10,-10,-10);

    if (is_iso == 0) {
      return lep;
    }
    else if (is_iso == +1) {
      if (is_mu) {
	Logger::get("lep_p4")->debug("---> should reco iso mu...");
	lep = ROOT::Math::PtEtaPhiMVector(mu_pt.at(ROOT::VecOps::ArgMax(tight_muons_mask),-2),
					  mu_eta.at(ROOT::VecOps::ArgMax(tight_muons_mask),-2),
					  mu_phi.at(ROOT::VecOps::ArgMax(tight_muons_mask),-2),
					  mu_mass.at(ROOT::VecOps::ArgMax(tight_muons_mask),-2));
      }
      else if (is_el) {
	Logger::get("lep_p4")->debug("---> should reco iso el...");
	lep = ROOT::Math::PtEtaPhiMVector(el_pt.at(ROOT::VecOps::ArgMax(tight_electrons_mask),-3),
					  el_eta.at(ROOT::VecOps::ArgMax(tight_electrons_mask),-3),
					  el_phi.at(ROOT::VecOps::ArgMax(tight_electrons_mask),-3),
					  el_mass.at(ROOT::VecOps::ArgMax(tight_electrons_mask),-3));
      }
    }
    else if (is_iso == -1) {
      if (is_mu) {
	Logger::get("lep_p4")->debug("---> should reco antiiso mu...");
	lep = ROOT::Math::PtEtaPhiMVector(mu_pt.at(ROOT::VecOps::ArgMax(antitight_muons_mask),-4),
					  mu_eta.at(ROOT::VecOps::ArgMax(antitight_muons_mask),-4),
					  mu_phi.at(ROOT::VecOps::ArgMax(antitight_muons_mask),-4),
					  mu_mass.at(ROOT::VecOps::ArgMax(antitight_muons_mask),-4));
      }
      else if (is_el) {
	Logger::get("lep_p4")->debug("---> should reco antiiso el...");
	lep = ROOT::Math::PtEtaPhiMVector(el_pt.at(ROOT::VecOps::ArgMax(antitight_electrons_mask),-5),
					  el_eta.at(ROOT::VecOps::ArgMax(antitight_electrons_mask),-5),
					  el_phi.at(ROOT::VecOps::ArgMax(antitight_electrons_mask),-5),
					  el_mass.at(ROOT::VecOps::ArgMax(antitight_electrons_mask),-5));
      }
    }
    else {
      lep = ROOT::Math::PtEtaPhiMVector(-6,-6,-6,-6);
    }

    Logger::get("final_lep")->debug("building p4 from lepton with {} {} {} {}",
				    lep.Pt(), lep.Eta(), lep.Phi(), lep.M());

    return lep;

  };

  auto df7 = df6.Define(str_lep_p4,
			lep_p4,
			{str_is_mu, str_is_el, str_is_iso,
			    str_tight_muons_mask, str_tight_electrons_mask,
			    str_antitight_muons_mask, str_antitight_electrons_mask,
			    str_mu_pt, str_mu_eta, str_mu_phi, str_mu_mass,
			    str_el_pt, str_el_eta, str_el_phi, str_el_mass}
			);

auto lep_sceta = [](const int is_el,
		    const int is_iso,
		    const ROOT::RVec<int> &tight_electrons_mask,
		    const ROOT::RVec<int> &antitight_electrons_mask,
		    const ROOT::RVec<float> &el_eta,
		    const ROOT::RVec<float> &el_detasc) {

  float lep_sceta = -10;

    if (is_iso == 0) {
      return lep_sceta;
    }
    else if (is_iso == +1) {
      if (is_el) {
	lep_sceta = el_eta.at(ROOT::VecOps::ArgMax(tight_electrons_mask),-5) + el_detasc.at(ROOT::VecOps::ArgMax(tight_electrons_mask),-5);
      }
    }
    else if (is_iso == -1) {
      if (is_el) {
	lep_sceta = el_eta.at(ROOT::VecOps::ArgMax(antitight_electrons_mask),-5) + el_detasc.at(ROOT::VecOps::ArgMax(antitight_electrons_mask),-5);
      }
    }
    return lep_sceta;
  };

  auto df7b = df7.Define(str_lep_sceta,
			lep_sceta,
			{str_is_el, str_is_iso,
			    str_tight_electrons_mask,
			    str_antitight_electrons_mask,
			    str_el_eta, str_el_detasc}
			);



  auto lep_charge= [](const int is_mu,
		      const int is_el,
		      const int is_iso,
		      const ROOT::RVec<int> &tight_muons_mask,
		      const ROOT::RVec<int> &tight_electrons_mask,
		      const ROOT::RVec<int> &antitight_muons_mask,
		      const ROOT::RVec<int> &antitight_electrons_mask,
		      const ROOT::RVec<int> &mu_charge,
		      const ROOT::RVec<int> &el_charge) {

    int charge = -10;

    if (is_iso == 0) {
      return charge;
    }
    else if (is_iso == +1) {
      if (is_mu) {
	charge = mu_charge.at(ROOT::VecOps::ArgMax(tight_muons_mask),-2);
      }
      else if (is_el) {
	charge = el_charge.at(ROOT::VecOps::ArgMax(tight_electrons_mask),-3);
      }
    }
    else if (is_iso == -1) {
      if (is_mu) {
	charge = mu_charge.at(ROOT::VecOps::ArgMax(antitight_muons_mask),-4);
      }
      else if (is_el) {
	charge = el_charge.at(ROOT::VecOps::ArgMax(antitight_electrons_mask),-5);
      }
    }
    else {
      charge = -6;
    }

    return charge;
  };


  auto df8 = df7b.Define(str_lep_charge,
			lep_charge,
			{str_is_mu, str_is_el, str_is_iso,
			    str_tight_muons_mask, str_tight_electrons_mask,
			    str_antitight_muons_mask, str_antitight_electrons_mask,
			    str_mu_charge,
			    str_el_charge}
			);


  return df8;
}





// helper function for minimizer constraint
double rad_py(double x, double lep_px) {
  return W_MASS*W_MASS + 4*lep_px*x;
}

// the delta plus function with the py nu plus solution
double min_fplus(double *par) {
  // par[0] = x, par[1] = lep_px, par[2] = lep_py, par[3] = lep_pt, par[4] = px_miss, par[5] = py_miss
  double r = rad_py(par[0],par[1]);
  double y = 0;
  //double res = 99999;
  //if (r>=0) {
  y = (W_MASS*W_MASS * par[2] + 2 * par[1] * par[2] * par[0] + W_MASS * par[3] * sqrt(r))/(2 * par[1]*par[1]);
  double res = sqrt((par[0]-par[4])*(par[0]-par[4]) + (y-par[5])*(y-par[5]));
  // }
  // else // FIXME: proper constraint in TMinuit?
  // res = 99999;
  return res;
}

// the delta minus function with the py nu minus solution
double min_fminus(double *par) {
  // par[0] = x, par[1] = lep_px, par[2] = lep_py, par[3] = lep_pt, par[4] = px_miss, par[5] = py_miss
  double r = rad_py(par[0],par[1]);
  double y = 0;
  double res = 99999;
  if (r>=0) {
    y = (W_MASS*W_MASS * par[2] + 2 * par[1] * par[2] * par[0] - W_MASS * par[3] * sqrt(r))/(2 * par[1]*par[1]);
    res = sqrt((par[0]-par[4])*(par[0]-par[4]) + (y-par[5])*(y-par[5]));
  }
  else
    res = 99999;
  return res;
}


// TMinuit fit function for the py nu plus solution
void fcn_plus(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  f = min_fplus(par);
}

// TMinuit fit function for the py nu minus solution
void fcn_minus(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
f = min_fminus(par);
}

ROOT::RDF::RNode ReconstructLeptonicW(ROOT::RDF::RNode df,
				      const std::string &str_lep_p4,
				      const std::string &str_met_p4,
				      const std::string &str_wlep_p4
				      ) {

  auto leptonicW = [](const ROOT::Math::PtEtaPhiMVector lep_p4,
		      const ROOT::Math::PtEtaPhiMVector met_p4) {

    auto wlep_p4 = ROOT::Math::PtEtaPhiMVector(-10,-10,-10,-10) ;

    if (lep_p4.Pt() < 0)
      return wlep_p4;

    double lep_e  = lep_p4.E();
    double lep_pt = lep_p4.Pt();
    double lep_px = lep_p4.Px();
    double lep_py = lep_p4.Py();
    double lep_pz = lep_p4.Pz();

    ROOT::Math::PtEtaPhiMVector nu_p4 = met_p4;
    double nu_e = met_p4.Pt();
    double nu_px = met_p4.Px();
    double nu_py = met_p4.Py();

    bool solution_is_real;


    Logger::get("wlep")->debug("building wlep p4 from lepton with E: {} px: {} py: {} pz: {}",
			       lep_e, lep_px, lep_py, lep_pz);
    Logger::get("wlep")->debug("building wlep p4 from pTmiss with E: {} px: {} py: {} pz: ???",
			       nu_e, nu_px, nu_py);

    // definition of the constant mu in Eq. 4.5 (here called alpha to not confuse mu and nu)
    // also converting p_T and cos dphi into p_x and p_y
    double alpha = (W_MASS*W_MASS)/2 + (lep_px*nu_px) + (lep_py*nu_py);

    // for p_z,nu there is a quadratic equation with up to two solutions as shown in Eq. 4.6 and A.7
    // (NOTE: there is a 'power of two' missing in the first denominator of Eq. 4.6)
    // first, check if we have complex solution, i.e. if the radicand is negative
    double rad = ((alpha*alpha * lep_pz*lep_pz)/(lep_pt*lep_pt*lep_pt*lep_pt)) - ((lep_e*lep_e * nu_e*nu_e - alpha*alpha)/(lep_pt*lep_pt));

    if(rad < 0){
      // complex solutions, in around 30% of all cases
      //    cout << "Complex neutrino p_z" << endl;

      // assumption: p_T^miss does not directly correspond to p_T,nu
      // set m_T^W to m^W, result is a quadratic equation for p_(x,y) that depends on p_(y,x)

      // save p_x^miss and p_y^miss as we need them later to determine the better solution

      Logger::get("wlep")->debug("complex solution");

      double px_miss = nu_px;
      double py_miss = nu_py;

      Logger::get("wlep")->debug("complex debug point 1");

      // initialize TMinuit with a maximum of 6 params for py nu plus and minus solution
      TMinuit *gMinuit_plus = new TMinuit(5);
      TMinuit *gMinuit_minus = new TMinuit(5);

      Logger::get("wlep")->debug("complex debug point 2");

      //    TMinuit *gMinuit_plus = new TMinuit(5);
      gMinuit_plus->SetFCN(fcn_plus);

      //    TMinuit *gMinuit_minus = new TMinuit(5);
      gMinuit_minus->SetFCN(fcn_minus);

      Logger::get("wlep")->debug("complex debug point 3");

      double arglist[10];
      int ierflg = 0;

      // no print
      arglist[0] = -1;
      gMinuit_plus->mnexcm("SET PRI", arglist, 1, ierflg);
      gMinuit_minus->mnexcm("SET PRI", arglist, 1, ierflg);

      // no warnings
      arglist[0] = -1;
      gMinuit_plus->mnexcm("SET NOW", arglist, 1, ierflg);
      gMinuit_minus->mnexcm("SET NOW", arglist, 1, ierflg);
      //    gMinuit_plus->mnexcm("SET WAR", arglist, 1, ierflg);
      //    gMinuit_minus->mnexcm("SET WAR", arglist, 1, ierflg);

      Logger::get("wlep")->debug("complex debug point 4");


      /*
      // set accuracy
      arglist[0] = 1.E-5L;
      gMinuit_plus->mnexcm("SET EPS", arglist, 1, ierflg);
      gMinuit_minus->mnexcm("SET EPS", arglist, 1, ierflg);
      */
      // set error
      arglist[0] = 1; // 0.5 for lnL, 1 for chi2
      // set the 1 sigma tolerance for the change in FCN
      // that determines when a function has been minimized
      gMinuit_plus->mnexcm("SET ERR", arglist, 1, ierflg);
      gMinuit_minus->mnexcm("SET ERR", arglist, 1, ierflg);

      // set strategy (0 is less accurate, 1 is default, 2 is more accurate)
      arglist[0] = 0;
      gMinuit_plus->mnexcm("SET STR", arglist, 1, ierflg);
      gMinuit_minus->mnexcm("SET STR", arglist, 1, ierflg);

      Logger::get("wlep")->debug("complex debug point 5");

      // set start value and range of x and fix all other params
      static double start_val = 0;
      static double step = 0.00001;
      double lower = -9999;
      double upper = 9999;
      if (lep_px > 0) {
	lower = -W_MASS*W_MASS/(4*lep_px) + 1e-5;
	start_val = lower + 1;
	//      cout << "lower: " << lower << endl;
      }
      if (lep_px < 0) {
	upper = -W_MASS*W_MASS/(4*lep_px) - 1e-5;
	start_val = upper - 1;
	//      cout << "upper: " << upper << endl;
      }
      gMinuit_plus->mnparm(0, "x", start_val, step, lower, upper, ierflg);
      gMinuit_plus->mnparm(1, "lep_px", lep_px, step, -999, 999, ierflg);
      gMinuit_plus->mnparm(2, "lep_py", lep_py, step, -999, 999, ierflg);
      gMinuit_plus->mnparm(3, "lep_pt", lep_pt, step, 0, 999, ierflg);
      gMinuit_plus->mnparm(4, "px_miss", px_miss, step, -999, 999, ierflg);
      gMinuit_plus->mnparm(5, "py_miss", py_miss, step, -999, 999, ierflg);
      gMinuit_plus->FixParameter(1);
      gMinuit_plus->FixParameter(2);
      gMinuit_plus->FixParameter(3);
      gMinuit_plus->FixParameter(4);
      gMinuit_plus->FixParameter(5);

      gMinuit_minus->mnparm(0, "x", start_val, step, lower, upper, ierflg);
      gMinuit_minus->mnparm(1, "lep_px", lep_px, step, -999, 999, ierflg);
      gMinuit_minus->mnparm(2, "lep_py", lep_py, step, -999, 999, ierflg);
      gMinuit_minus->mnparm(3, "lep_pt", lep_pt, step, 0, 999, ierflg);
      gMinuit_minus->mnparm(4, "px_miss", px_miss, step, -999, 999, ierflg);
      gMinuit_minus->mnparm(5, "py_miss", py_miss, step, -999, 999, ierflg);
      gMinuit_minus->FixParameter(1);
      gMinuit_minus->FixParameter(2);
      gMinuit_minus->FixParameter(3);
      gMinuit_minus->FixParameter(4);
      gMinuit_minus->FixParameter(5);

      Logger::get("wlep")->debug("complex debug point 6");

      // now ready for minimization step
      arglist[0] = 500000; // maximum number of iterations
      arglist[1] = 1.; // related to errors
      gMinuit_plus->mnexcm("MIGARD", arglist, 2, ierflg);
      gMinuit_minus->mnexcm("MIGRAD", arglist, 2, ierflg);

      // obtain fit results and calculate values of delta minus and delta plus functions
      // choose solution that leads to a smaller delta value
      double x_plus, x_pluserr;
      double d_plus;
      gMinuit_plus->GetParameter(0,x_plus, x_pluserr);
      double par_plus[6] = {x_plus,lep_px,lep_py,lep_pt,px_miss,py_miss};
      d_plus = min_fplus(par_plus);
      //    cout << "Fit result plus: x=" << x_plus << " " << "d(x)=" << d_plus << endl;

      double x_minus, x_minuserr;
      double d_minus;
      gMinuit_minus->GetParameter(0,x_minus, x_minuserr);
      double par_minus[6] = {x_minus,lep_px,lep_py,lep_pt,px_miss,py_miss};
      d_minus = min_fminus(par_minus);
      //    cout << "Fit result minus: x=" << x_minus << " d(x)=" << d_minus << endl;

      Logger::get("wlep")->debug("complex debug point 7");

      double nu_pxnew, nu_pynew, r_new;
      if (d_plus<d_minus){
	nu_pxnew = x_plus;
	r_new = rad_py(nu_pxnew,lep_px);
	nu_pynew = (W_MASS*W_MASS * lep_py + 2 * lep_px * lep_py * nu_pxnew + W_MASS * lep_pt * sqrt(r_new))/(2 * lep_px * lep_px);
      }
      else{
	nu_pxnew = x_minus;
	r_new = rad_py(nu_pxnew,lep_px);
	nu_pynew = (W_MASS*W_MASS * lep_py + 2 * lep_px * lep_py * nu_pxnew - W_MASS * lep_pt * sqrt(r_new))/(2 * lep_px * lep_px);
      }
      // calculate new nu pz (only one solution with fixed px and py)
      double nu_pznew = lep_pz / (lep_pt*lep_pt) * ((W_MASS*W_MASS / 2) + (lep_px * nu_pxnew) + (lep_py * nu_pynew));
      //    cout << "new nu px: " << nu_pxnew << ", new nu py: " << nu_pynew << ", new nu pz: " << nu_pznew << endl;

      // set 4 momenta of neutrino and W boson
      nu_p4.SetPxPyPzE(nu_pxnew,nu_pynew,nu_pznew,sqrt(nu_pxnew*nu_pxnew + nu_pynew*nu_pynew + nu_pznew*nu_pznew));

      Logger::get("wlep")->debug("complex debug point 8");

      delete gMinuit_plus;
      delete gMinuit_minus;

    }

    else { // two real solutions for pz nu
      //    cout << "Two neutrino pz solutions" << endl;
      Logger::get("wlep")->debug("real solution");
      double sol1, sol2, nu_pz;
      sol1 = lep_pz * alpha / (lep_pt * lep_pt) + sqrt(rad);
      sol2 = lep_pz * alpha / (lep_pt * lep_pt) - sqrt(rad);

      // choose the smaller pz solution
      if (abs(sol1) < abs(sol2)) {
	nu_pz = sol1;
      } else {
	nu_pz = sol2;
      }

      // set 4 momenta of neutrino and W boson
      nu_p4.SetPxPyPzE(nu_px, nu_py, nu_pz, sqrt(nu_px * nu_px + nu_py * nu_py + nu_pz * nu_pz));
      solution_is_real = true;
    }

    wlep_p4 = lep_p4 + nu_p4;

    Logger::get("wlep")->debug("final leptonic W boson pT: {} eta: {} phi: {} mass: {}",
			       wlep_p4.Pt(), wlep_p4.Eta(), wlep_p4.Phi(), wlep_p4.M());

    return wlep_p4;

  };



  return df.Define(str_wlep_p4,
		   leptonicW,
		   {str_lep_p4, str_met_p4}
		   );

}



ROOT::RDF::RNode JetSelection(ROOT::RDF::RNode df,
			      const int &njets,
			      const int &nbjets,
			      const std::string &str_good_jets_mask,
			      const std::string &str_good_bjets_mask
			      ) {

  auto df2 = df.Filter([njets, nbjets](const ROOT::RVec<int> &good_jets_mask,
				       const ROOT::RVec<int> &good_bjets_mask) {

			 int nj = ROOT::VecOps::Sum(good_jets_mask,0);
			 int nbj = ROOT::VecOps::Sum(good_bjets_mask,0);

			 return (nj == njets) && (nbj == nbjets);
		       },
		       {str_good_jets_mask, str_good_bjets_mask},
		       "jet selection and b jet selection"
			);

  return df2;
}



ROOT::RDF::RNode TopReco(ROOT::RDF::RNode df,
			 const std::string &str_wlep_p4,
			 const std::string &str_n_nonbjets,
			 const std::string &str_nonbjet_p4_1,
			 const std::string &str_nonbjet_btag_1,
			 const std::string &str_nonbjet_p4_2,
			 const std::string &str_nonbjet_btag_2,
			 const std::string &str_n_bjets,
			 const std::string &str_bjet_p4_1,
			 const std::string &str_bjet_btag_1,
			 const std::string &str_bjet_p4_2,
			 const std::string &str_bjet_btag_2,
			 const std::string &str_is_reco,
			 const std::string &str_is_jjb,
			 const std::string &str_is_jjbb,
			 const std::string &str_is_jjjb,
			 const std::string &str_is_jjjbb,
			 const std::string &str_reco_p4s,
			 const std::string &str_top_p4,
			 const std::string &str_tb_p4,
			 const std::string &str_sb_p4
			 ) {

  auto df2 = df.Define(str_is_jjb,
		       [](const int n_nonbjets,
			  const int n_bjets){
			 return int((n_nonbjets + n_bjets) == 2 && n_bjets == 1);
		       },
		       {str_n_nonbjets, str_n_bjets}
		       );

  auto df3 = df2.Define(str_is_jjbb,
			[](const int n_nonbjets,
			   const int n_bjets){
			  return int((n_nonbjets + n_bjets) == 2 && n_bjets == 2);
			},
			{str_n_nonbjets, str_n_bjets}
			);

  auto df4 = df3.Define(str_is_jjjb,
			[](const int n_nonbjets,
			   const int n_bjets){
			  return int((n_nonbjets + n_bjets) == 3 && n_bjets == 1);
			},
			{str_n_nonbjets, str_n_bjets}
			);

  auto df5 = df4.Define(str_is_jjjbb,
			[](const int n_nonbjets,
			   const int n_bjets){
			  return int((n_nonbjets + n_bjets) == 3 && n_bjets == 2);
			},
			{str_n_nonbjets, str_n_bjets}
			);

  auto df6 = df5.Define(str_is_reco,
		       [](const int is_jjb,
			  const int is_jjbb,
			  const int is_jjjb,
			  const int is_jjjbb
			  ){
			 return is_jjb + is_jjbb + is_jjjb + is_jjjbb;
		       },
		       {str_is_jjb, str_is_jjbb, str_is_jjjb, str_is_jjjbb}
		       );


  auto top_reco = [](const int is_reco,
		     const int is_jjb,
		     const int is_jjbb,
		     const int is_jjjb,
		     const int is_jjjbb,
		     const ROOT::Math::PtEtaPhiMVector wlep_p4,
		     const ROOT::Math::PtEtaPhiMVector nonbjet_p4_1,
		     const float nonbjet_btag_1,
		     const ROOT::Math::PtEtaPhiMVector nonbjet_p4_2,
		     const float nonbjet_btag_2,
		     const ROOT::Math::PtEtaPhiMVector bjet_p4_1,
		     const float bjet_btag_1,
		     const ROOT::Math::PtEtaPhiMVector bjet_p4_2,
		     const float bjet_btag_2
		     ){


    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> reco_vec {
      ROOT::Math::PtEtaPhiMVector(-10,-10,-10,-10), // top_p4
      ROOT::Math::PtEtaPhiMVector(-10,-10,-10,-10), // tb_p4
      ROOT::Math::PtEtaPhiMVector(-10,-10,-10,-10)  // sb_p4
    };


    if (wlep_p4.Pt() < 0)
      return reco_vec;

    if (!is_reco)
      return reco_vec;

    if (is_jjb) { // 2j1b
      reco_vec[0] = wlep_p4 + bjet_p4_1;
      reco_vec[1] = bjet_p4_1;
      reco_vec[2] = nonbjet_p4_1;
    }
    else if (is_jjbb) { // 2j2b
      auto cand1_p4 = wlep_p4 + bjet_p4_1;
      auto cand2_p4 = wlep_p4 + bjet_p4_2;
      if (abs(cand1_p4.M() - TOP_MASS) < abs(cand2_p4.M() - TOP_MASS)) {
	reco_vec[0] = cand1_p4;
	reco_vec[1] = bjet_p4_1;
	reco_vec[2] = bjet_p4_2;
      }
      else {
	reco_vec[0] = cand2_p4;
	reco_vec[1] = bjet_p4_2;
	reco_vec[2] = bjet_p4_1;
      }
    }
    else if (is_jjjb) { // 3j1b
      reco_vec[0] = wlep_p4 + bjet_p4_1;
      reco_vec[1] = bjet_p4_1;
      if (nonbjet_btag_1 > nonbjet_btag_2)
	reco_vec[2] = nonbjet_p4_1;
      else
	reco_vec[2] = nonbjet_p4_2;
    }
    else if (is_jjjbb) { // 3j2b
      auto cand1_p4 = wlep_p4 + bjet_p4_1;
      auto cand2_p4 = wlep_p4 + bjet_p4_2;
      if (abs(cand1_p4.M() - TOP_MASS) < abs(cand2_p4.M() - TOP_MASS)) {
	reco_vec[0] = cand1_p4;
	reco_vec[1] = bjet_p4_1;
	reco_vec[2] = bjet_p4_2;
      }
      else {
	reco_vec[0] = cand2_p4;
	reco_vec[1] = bjet_p4_2;
	reco_vec[2] = bjet_p4_1;
      }
    }

    return reco_vec;
  };


  auto df7 = df6.Define(str_reco_p4s,
			top_reco,
			{str_is_reco,
			    str_is_jjb,
			    str_is_jjbb,
			    str_is_jjjb,
			    str_is_jjjbb,
			    str_wlep_p4,
			    str_nonbjet_p4_1,
			    str_nonbjet_btag_1,
			    str_nonbjet_p4_2,
			    str_nonbjet_btag_2,
			    str_bjet_p4_1,
			    str_bjet_btag_1,
			    str_bjet_p4_2,
			    str_bjet_btag_2
			    }
			);

  auto df8 = df7.Define(str_top_p4,
			[](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &vec) {
			 return vec[0];
			},
			{str_reco_p4s}
			);

  auto df9 = df8.Define(str_tb_p4,
			[](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &vec) {
			 return vec[1];
			},
			{str_reco_p4s}
			);

  auto df10 = df9.Define(str_sb_p4,
			 [](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &vec) {
			 return vec[2];
			},
			{str_reco_p4s}
			);

  return df10;
}



ROOT::RDF::RNode DNNQuantities(ROOT::RDF::RNode df,
			       const std::string &str_is_reco,
			       const std::string &str_lep_p4,
			       const std::string &str_met_p4,
			       const std::string &str_wlep_p4,
			       const std::string &str_nonbjet_p4_1,
			       const std::string &str_nonbjet_btag_1,
			       const std::string &str_nonbjet_p4_2,
			       const std::string &str_nonbjet_btag_2,
			       const std::string &str_bjet_p4_1,
			       const std::string &str_bjet_btag_1,
			       const std::string &str_bjet_p4_2,
			       const std::string &str_bjet_btag_2,
			       const std::string &str_top_p4,
			       const std::string &str_tb_p4,
			       const std::string &str_sb_p4,
			       const std::string &str_good_jetslowpt_mask,
			       const std::string &str_jet_pt,
			       const std::string &str_jet_eta,
			       const std::string &str_jet_phi,
			       const std::string &str_jet_mass,
			       const std::string &str_dphi_top_tb,
			       const std::string &str_deta_top_sb,
			       const std::string &str_dphi_tb_sb,
			       const std::string &str_deta_lep_tb,
			       const std::string &str_m_lep_sb,
			       const std::string &str_pt_tb_sb,
			       const std::string &str_costhetastar,
			       const std::string &str_sumht,
			       const std::string &str_wolfram,
			       const std::string &str_deta_topsb_tb
			       ) {

  auto df2 = df.Define(str_dphi_top_tb,
		       [](const int reco,
			  const ROOT::Math::PtEtaPhiMVector &top,
			  const ROOT::Math::PtEtaPhiMVector &tb) {
			 if (!reco) {return -10.0;}
			 return ROOT::Math::VectorUtil::DeltaPhi(top, tb);
		       },
		       {str_is_reco, str_top_p4, str_tb_p4}
		       );

  auto df3 = df2.Define(str_deta_top_sb,
			[](const int reco,
			   const ROOT::Math::PtEtaPhiMVector &top,
			   const ROOT::Math::PtEtaPhiMVector &sb) {
			  if (!reco) {return -10.0;}
			  return abs(top.Eta() - sb.Eta());
			},
			{str_is_reco, str_top_p4, str_sb_p4}
			);

  auto df4 = df3.Define(str_dphi_tb_sb,
			[](const int reco,
			   const ROOT::Math::PtEtaPhiMVector &tb,
			   const ROOT::Math::PtEtaPhiMVector &sb) {
			 if (!reco) {return -10.0;}
			  return ROOT::Math::VectorUtil::DeltaPhi(tb, sb);
			},
			{str_is_reco, str_tb_p4, str_sb_p4}
			);

  auto df5 = df4.Define(str_deta_lep_tb,
			[](const int reco,
			   const ROOT::Math::PtEtaPhiMVector &lep,
			   const ROOT::Math::PtEtaPhiMVector &tb) {
			 if (!reco) {return -10.0;}
			  return abs(lep.Eta() - tb.Eta());
			},
			{str_is_reco, str_lep_p4, str_tb_p4}
			);

  auto df6 = df5.Define(str_m_lep_sb,
			[](const int reco,
			   const ROOT::Math::PtEtaPhiMVector &lep,
			   const ROOT::Math::PtEtaPhiMVector &sb) {
			 if (!reco) {return -10.0;}
			  return (lep + sb).M();
			},
			{str_is_reco, str_lep_p4, str_sb_p4}
			);

  auto df7 = df6.Define(str_pt_tb_sb,
			[](const int reco,
			   const ROOT::Math::PtEtaPhiMVector &tb,
			   const ROOT::Math::PtEtaPhiMVector &sb) {
			 if (!reco) {return -10.0;}
			 return (tb + sb).Pt();
			},
			{str_is_reco, str_tb_p4, str_sb_p4}
			);

  auto df8 = df7.Define(str_costhetastar,
			[](const int reco,
			   const ROOT::Math::PtEtaPhiMVector &top,
			   const ROOT::Math::PtEtaPhiMVector &sb,
			   const ROOT::Math::PtEtaPhiMVector &lep) {
			  if (!reco) {
			    return -10.0;}
			  double costhetastar = 0;

			  auto top_boost_vec = ROOT::Math::Cartesian3D(top.X()/top.T(), top.Y()/top.T(), top.Z()/top.T());
			  ROOT::Math::Boost top_boost(top_boost_vec);
			  top_boost.Invert();
			  ROOT::Math::PtEtaPhiMVector lep_boosted = top_boost(lep);
			  ROOT::Math::PtEtaPhiMVector sb_boosted = top_boost(sb);
			  costhetastar = lep_boosted.Vect().Dot(sb_boosted.Vect()) / sqrt(lep_boosted.Vect().Mag2()*sb_boosted.Vect().Mag2());

			  Logger::get("DNN_costhetastar")->debug("top_boost {} {} {}", top_boost_vec.X(), top_boost_vec.Y(), top_boost_vec.Z());
			  Logger::get("DNN_costhetastar")->debug("sb boosted {} {} {} {}", sb_boosted.Pt(), sb_boosted.Eta(), sb_boosted.Phi(), sb_boosted.M());
			  Logger::get("DNN_costhetastar")->debug("lep boosted {} {} {} {}", lep_boosted.Pt(), lep_boosted.Eta(), lep_boosted.Phi(), lep_boosted.M());
			  Logger::get("DNN_costhetastar")->debug("costhetastar {}", costhetastar);

			  return costhetastar;
			},
			{str_is_reco, str_top_p4, str_sb_p4, str_lep_p4}
			);


  auto df9 = df8.Define(str_sumht,
			[](const int reco,
			   const ROOT::Math::PtEtaPhiMVector &tb,
			   const ROOT::Math::PtEtaPhiMVector &sb,
			   const ROOT::Math::PtEtaPhiMVector &lep,
			   const ROOT::Math::PtEtaPhiMVector &met) {
			  if (!reco) {return -10.0;}
			  return tb.Pt() + sb.Pt() + lep.Pt() + met.Pt();
			},
			{str_is_reco, str_tb_p4, str_sb_p4, str_lep_p4, str_met_p4}
			);


  auto df10 = df9.Define(str_wolfram,
			[](const int reco,
			   const ROOT::RVec<int> &jet_mask,
			   const ROOT::RVec<float> &jet_pt,
			   const ROOT::RVec<float> &jet_eta,
			   const ROOT::RVec<float> &jet_phi,
			   const ROOT::RVec<float> &jet_mass) {
			  if (!reco) {return -10.0;}

			  Logger::get("DNN_wolfram")->debug("eval sum_e...");
			  Logger::get("DNN_wolfram")->debug("jet mask {}", jet_mask);

			  float sum_e = 0;
			  for (std::size_t index = 0; index < jet_pt.size(); index++) {
			      Logger::get("DNN_wolfram")->debug("checking jet {}", index);
			    if (jet_mask.at(index)) {
			      sum_e += (ROOT::Math::PtEtaPhiMVector(jet_pt.at(index,0),
								    jet_eta.at(index,0),
								    jet_phi.at(index,0),
								    jet_mass.at(index,0))
					).E();
			      Logger::get("DNN_wolfram")->debug("sum_e now {} after adding jet {}", sum_e, index);

			    }
			  }

			  Logger::get("DNN_wolfram")->debug("now momenta");

			  double h3 = 0;

			  for (std::size_t index = 0; index < jet_pt.size(); index++) {
			    for (std::size_t index2 = index + 1; index2 < jet_pt.size(); index2++) {
			      if (jet_mask.at(index) && jet_mask.at(index2)) {

				auto jet_a = ROOT::Math::PtEtaPhiMVector(jet_pt[index],jet_eta[index],jet_phi[index],jet_mass[index]);
				auto jet_b = ROOT::Math::PtEtaPhiMVector(jet_pt[index2],jet_eta[index2],jet_phi[index2],jet_mass[index2]);

				auto jet_aXYZT = ROOT::Math::XYZTVector(jet_a.X(),jet_a.Y(),jet_a.Z(),jet_a.T());
				auto jet_bXYZT = ROOT::Math::XYZTVector(jet_b.X(),jet_b.Y(),jet_b.Z(),jet_b.T());

				float costh = ROOT::Math::VectorUtil::CosTheta(jet_aXYZT,jet_bXYZT);
				float p3 = 0.5*(5.0*costh*costh - 3.0*costh);
				float pipj = jet_aXYZT.P()*jet_bXYZT.P();

				h3 += (pipj/(sum_e*sum_e))*p3;

			      }
			    }
			  }

			  return h3;

			},
			 {str_is_reco, str_good_jetslowpt_mask, str_jet_pt, str_jet_eta, str_jet_phi, str_jet_mass}
			 );

  auto df11 = df10.Define(str_deta_topsb_tb,
			  [](const int reco,
			     const ROOT::Math::PtEtaPhiMVector &tb,
			     const ROOT::Math::PtEtaPhiMVector &sb,
			     const ROOT::Math::PtEtaPhiMVector &wlep) {
			    if (!reco) {return -10.0;}
			    return abs((sb + wlep).Eta() - tb.Eta());
			  },
			  {str_is_reco, str_tb_p4, str_sb_p4, str_wlep_p4}
			  );


  return df11;



}





void sf_from_root_file(TH2D* h, int nbinsx, int nbinsy, float pt, float eta, int var, double* sf) {
// double sf_from_root_file(TH2D *h,float pt, float eta, int var, int nbinsx, int nbinsy, int ymax) {

  Logger::get("trigger_func")->debug("complex debug point 1");

  Logger::get("trigger_func")->debug("hist bins {} {}", h->GetXaxis()->GetNbins(), h->GetYaxis()->GetNbins());


  int xbin_index = h->GetXaxis()->FindBin(eta);
  int ybin_index = h->GetYaxis()->FindBin(pt);

  // int ybin_max_index = h->GetYaxis()->FindBin(ymax);
  Logger::get("trigger_func")->debug("complex debug point 2");

  if(xbin_index > nbinsx) xbin_index = nbinsx;
  if(ybin_index > nbinsy) ybin_index = nbinsy;
  // if(ybin_max_index>nbinsy) ybin_max_index=nbinsy;
  Logger::get("trigger_func")->debug("complex debug point 3");

  if(xbin_index == 0) xbin_index = 1;
  if(ybin_index == 0) ybin_index = 1;

  // int bin_index;
  Logger::get("trigger_func")->debug("complex debug point 4");

  const int bin_index = h->GetBin(xbin_index, ybin_index, 0);
  // const int bin_max_index = h->GetBin(xbin_index,ybin_max_index,0);
  Logger::get("trigger_func")->debug("complex debug point 5");

  *sf = 1.;

  Logger::get("trigger_func")->debug("complex debug point 6");

  if (var == 0)
    *sf = h->GetBinContent(bin_index);
  if (var == +1)
    *sf = h->GetBinContent(bin_index) + h->GetBinError(bin_index);
  if (var == -1)
    *sf = h->GetBinContent(bin_index) - h->GetBinError(bin_index);
}

ROOT::RDF::RNode LeptonScaleFactors(ROOT::RDF::RNode df,
				    const std::string &str_lep_pt,
				    const std::string &str_lep_eta,
				    const std::string &str_lep_sceta,
				    const std::string &str_lep_is_mu,
				    const std::string &str_lep_is_el,
				    const std::string &str_lep_is_iso,
				    const std::string &str_lep_sf_mu_trigger_nom,
				    const std::string &str_lep_sf_mu_trigger_up,
				    const std::string &str_lep_sf_mu_trigger_down,
				    const std::string &str_lep_sf_mu_iso_nom,
				    const std::string &str_lep_sf_mu_iso_up,
				    const std::string &str_lep_sf_mu_iso_down,
				    const std::string &str_lep_sf_mu_id_nom,
				    const std::string &str_lep_sf_mu_id_up,
				    const std::string &str_lep_sf_mu_id_down,
				    const std::string &str_lep_sf_el_trigger_nom,
				    const std::string &str_lep_sf_el_trigger_up,
				    const std::string &str_lep_sf_el_trigger_down,
				    const std::string &str_lep_sf_el_id_nom,
				    const std::string &str_lep_sf_el_id_up,
				    const std::string &str_lep_sf_el_id_down,
				    const std::string &str_lep_sf_el_reco_nom,
				    const std::string &str_lep_sf_el_reco_up,
				    const std::string &str_lep_sf_el_reco_down,
				    const std::string &mu_sf_era,
				    const std::string &mu_trigger_sf_file,
				    const std::string &mu_trigger_sf_file_syst,
				    const std::string &mu_trigger_sf_name,
				    const std::string &mu_trigger_sf_name_syst,
				    const std::string &mu_iso_sf_file,
				    const std::string &mu_iso_sf_file_syst,
				    const std::string &mu_iso_sf_name,
				    const std::string &mu_iso_sf_name_syst,
				    const std::string &mu_sf_file,
				    const std::string &mu_id_sf_name,
				    const std::string &el_sf_era,
				    const std::string &el_trigger_sf_file,
				    const std::string &el_trigger_sf_file_syst,
				    const std::string &el_trigger_sf_name,
				    const std::string &el_trigger_sf_name_syst,
				    const std::string &el_sf_file,
				    const std::string &el_id_sf_name
				    ) {


  Logger::get("lepsf_muonTriggerSF")->debug("Setting up functions for muon trigger sf");
  Logger::get("lepsf_muonTriggerSF")->debug("Trigger - File {}", mu_trigger_sf_file);
  Logger::get("lepsf_muonTriggerSF")->debug("Trigger - Name {}", mu_trigger_sf_name);
  auto evaluator_mu_trigger = correction::CorrectionSet::from_file(mu_trigger_sf_file)->at(mu_trigger_sf_name);
  auto mu_trigger_nom = [evaluator_mu_trigger](const float &pt,
					       const float &eta,
					       const int &is_mu,
					       const int &is_el,
					       const int &is_iso){


    double sf = 1.;
    if (is_mu == 0 || is_iso != +1)return sf;
    Logger::get("lepsf_muonTriggerSF")->debug("Muon - pt {}, eta {}", pt, eta);
    if (pt >= 0.0)
      sf = evaluator_mu_trigger->evaluate({std::abs(eta), pt});
    Logger::get("lepsf_muonTriggerSF")->debug("Trigger - sf {}", sf);
    return sf;

  };
  Logger::get("lepsf_muonTriggerSF")->debug("Trigger syst - File {}", mu_trigger_sf_file_syst);
  Logger::get("lepsf_muonTriggerSF")->debug("Trigger syst - Name {}", mu_trigger_sf_name_syst);
  auto evaluator_mu_trigger_syst = correction::CorrectionSet::from_file(mu_trigger_sf_file_syst)->at(mu_trigger_sf_name_syst);
  auto mu_trigger_up = [evaluator_mu_trigger, evaluator_mu_trigger_syst](const float &pt,
									 const float &eta,
									 const int &is_mu,
									 const int &is_el,
									 const int &is_iso){


    double sf = 1.;
    if (is_mu == 0 || is_iso != +1)return sf;
    if (pt >= 0.0)
      sf = evaluator_mu_trigger->evaluate({std::abs(eta), pt}) + evaluator_mu_trigger_syst->evaluate({std::abs(eta), pt});
    Logger::get("lepsf_muonTriggerSF")->debug("Trigger up - sf {}", sf);
    return sf;
  };
  auto mu_trigger_down = [evaluator_mu_trigger, evaluator_mu_trigger_syst](const float &pt,
									   const float &eta,
									   const int &is_mu,
									   const int &is_el,
									   const int &is_iso){


    double sf = 1.;
    if (is_mu == 0 || is_iso != +1)return sf;
    if (pt >= 0.0)
      sf = evaluator_mu_trigger->evaluate({std::abs(eta), pt}) - evaluator_mu_trigger_syst->evaluate({std::abs(eta), pt});
    Logger::get("lepsf_muonTriggerSF")->debug("Trigger down - sf {}", sf);
    return sf;
  };


  Logger::get("lepsf_muonIsoSF")->debug("Setting up functions for muon iso sf");
  Logger::get("lepsf_muonIsoSF")->debug("ISO - File {}", mu_iso_sf_file);
  Logger::get("lepsf_muonIsoSF")->debug("ISO - Name {}", mu_iso_sf_name);
  auto evaluator_mu_iso = correction::CorrectionSet::from_file(mu_iso_sf_file)->at(mu_iso_sf_name);
  auto mu_iso_nom = [evaluator_mu_iso](const float &pt,
					       const float &eta,
					       const int &is_mu,
					       const int &is_el,
					       const int &is_iso){


    double sf = 1.;
    if (is_mu == 0 || is_iso != +1)return sf;
    if (pt >= 0.0)
      sf = evaluator_mu_iso->evaluate({std::abs(eta), pt});
    Logger::get("lepsf_muonIsoSF")->debug("Iso - sf {}", sf);
    return sf;

  };
  Logger::get("lepsf_muonisoSF")->debug("iso syst - File {}", mu_iso_sf_file_syst);
  Logger::get("lepsf_muonisoSF")->debug("iso syst - Name {}", mu_iso_sf_name_syst);
  auto evaluator_mu_iso_syst = correction::CorrectionSet::from_file(mu_iso_sf_file_syst)->at(mu_iso_sf_name_syst);
  auto mu_iso_up = [evaluator_mu_iso, evaluator_mu_iso_syst](const float &pt,
									 const float &eta,
									 const int &is_mu,
									 const int &is_el,
									 const int &is_iso){


    double sf = 1.;
    if (is_mu == 0 || is_iso != +1)return sf;
    if (pt >= 0.0)
      sf = evaluator_mu_iso->evaluate({std::abs(eta), pt}) + evaluator_mu_iso_syst->evaluate({std::abs(eta), pt});
    Logger::get("lepsf_muonIsoSF")->debug("Iso up - sf {}", sf);
    return sf;
  };
  auto mu_iso_down = [evaluator_mu_iso, evaluator_mu_iso_syst](const float &pt,
									   const float &eta,
									   const int &is_mu,
									   const int &is_el,
									   const int &is_iso){


    double sf = 1.;
    if (is_mu == 0 || is_iso != +1)return sf;
    if (pt >= 0.0)
      sf = evaluator_mu_iso->evaluate({std::abs(eta), pt}) - evaluator_mu_iso_syst->evaluate({std::abs(eta), pt});
    Logger::get("lepsf_muonIsoSF")->debug("Iso down - sf {}", sf);
    return sf;
  };


  Logger::get("lepsf_muonIdSF")->debug("Setting up functions for muon id sf");
  Logger::get("lepsf_muonIdSF")->debug("ID - File {}", mu_sf_file);
  Logger::get("lepsf_muonIdSF")->debug("ID - Name {}", mu_id_sf_name);
  Logger::get("lepsf_muonIdSF")->debug("ID - Era {}", mu_sf_era);
  auto evaluator_mu_id = correction::CorrectionSet::from_file(mu_sf_file)->at(mu_id_sf_name);

  auto mu_id_nom = [evaluator_mu_id, mu_sf_era](const float &pt,
					     const float &eta,
					     const int &is_mu,
					     const int &is_el,
					     const int &is_iso
					     ) {
    double sf = 1.;
    if (is_mu == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_mu_id->evaluate({mu_sf_era, std::abs(eta), pt, "sf"});
    Logger::get("lepsf_muonIdSF")->debug("ID - sf {}", sf);
    return sf;
  };

  auto mu_id_up = [evaluator_mu_id, mu_sf_era](const float &pt,
					     const float &eta,
					     const int &is_mu,
					     const int &is_el,
					     const int &is_iso
					     ) {
    double sf = 1.;
    if (is_mu == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_mu_id->evaluate({mu_sf_era, std::abs(eta), pt, "systup"});
    Logger::get("lepsf_muonIdSF")->debug("ID up - sf {}", sf);
    return sf;
  };

  auto mu_id_down = [evaluator_mu_id, mu_sf_era](const float &pt,
					     const float &eta,
					     const int &is_mu,
					     const int &is_el,
					     const int &is_iso
					     ) {
    double sf = 1.;
    if (is_mu == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_mu_id->evaluate({mu_sf_era, std::abs(eta), pt, "systdown"});
    Logger::get("lepsf_muonIdSF")->debug("ID down - sf {}", sf);
    return sf;
  };


  auto df1 = df.Define(str_lep_sf_mu_trigger_nom,
  		       mu_trigger_nom,
  		       {str_lep_pt, str_lep_eta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
  		       );
  auto df2 = df1.Define(str_lep_sf_mu_trigger_up,
			mu_trigger_up,
			{str_lep_pt, str_lep_eta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);
  auto df3 = df2.Define(str_lep_sf_mu_trigger_down,
			mu_trigger_down,
			{str_lep_pt, str_lep_eta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);
  auto df4 = df3.Define(str_lep_sf_mu_iso_nom,
		       mu_iso_nom,
		       {str_lep_pt, str_lep_eta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
		       );
  auto df5 = df4.Define(str_lep_sf_mu_iso_up,
			mu_iso_up,
			{str_lep_pt, str_lep_eta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);
  auto df6 = df5.Define(str_lep_sf_mu_iso_down,
			mu_iso_down,
			{str_lep_pt, str_lep_eta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);
  auto df7 = df6.Define(str_lep_sf_mu_id_nom,
		       mu_id_nom,
		       {str_lep_pt, str_lep_eta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
		       );
  auto df8 = df7.Define(str_lep_sf_mu_id_up,
			mu_id_up,
			{str_lep_pt, str_lep_eta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);
  auto df9 = df8.Define(str_lep_sf_mu_id_down,
			mu_id_down,
			{str_lep_pt, str_lep_eta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);






  Logger::get("lepsf_electronTriggerSF")->debug("Setting up functions for electron trigger sf");
  Logger::get("lepsf_electronTriggerSF")->debug("Trigger - File {}", el_trigger_sf_file);
  Logger::get("lepsf_electronTriggerSF")->debug("Trigger - Name {}", el_trigger_sf_name);
  auto evaluator_el_trigger = correction::CorrectionSet::from_file(el_trigger_sf_file)->at(el_trigger_sf_name);
  auto el_trigger_nom = [evaluator_el_trigger](const float &pt,
					       const float &sceta,
					       const int &is_mu,
					       const int &is_el,
					       const int &is_iso){


    double sf = 1.;
    if (is_el == 0 || is_iso != +1) return sf;
    Logger::get("lepsf_electronTriggerSF")->debug("Electron - pt {}, sceta {}", pt, sceta);
    if (pt >= 0.0)
      sf = evaluator_el_trigger->evaluate({sceta, pt});
    Logger::get("lepsf_electronTriggerSF")->debug("Trigger - sf {}", sf);
    return sf;

  };


  auto evaluator_el_trigger_syst = correction::CorrectionSet::from_file(el_trigger_sf_file_syst)->at(el_trigger_sf_name_syst);
  auto el_trigger_up = [evaluator_el_trigger, evaluator_el_trigger_syst](const float &pt,
									 const float &sceta,
									 const int &is_mu,
									 const int &is_el,
									 const int &is_iso){


    double sf = 1.;
    if (is_el == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_el_trigger->evaluate({sceta, pt}) + evaluator_el_trigger_syst->evaluate({sceta, pt});
    Logger::get("lepsf_electronTriggerSF")->debug("Trigger up - sf {}", sf);
    return sf;
  };
  auto el_trigger_down = [evaluator_el_trigger, evaluator_el_trigger_syst](const float &pt,
									   const float &sceta,
									   const int &is_mu,
									   const int &is_el,
									   const int &is_iso){


    double sf = 1.;
    if (is_el == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_el_trigger->evaluate({sceta, pt}) - evaluator_el_trigger_syst->evaluate({sceta, pt});
    Logger::get("lepsf_electronTriggerSF")->debug("Trigger up - sf {}", sf);
    return sf;
  };


  Logger::get("lepsf_electronSF")->debug("Setting up functions for electron id+reco sf");
  Logger::get("lepsf_electronSF")->debug("SF - File {}", el_sf_file);
  Logger::get("lepsf_electronSF")->debug("SF - Name {}", el_id_sf_name);
  Logger::get("lepsf_electronSF")->debug("SF - Era {}", el_sf_era);
  auto evaluator_el_sf = correction::CorrectionSet::from_file(el_sf_file)->at(el_id_sf_name);

  auto el_id_nom = [evaluator_el_sf, el_sf_era](const float &pt,
					     const float &sceta,
					     const int &is_mu,
					     const int &is_el,
					     const int &is_iso
					     ) {
    double sf = 1.;
    if (is_el == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_el_sf->evaluate({el_sf_era, "sf", "Tight", sceta, pt});
    Logger::get("lepsf_electronIdSF")->debug("ID - sf {}", sf);
    return sf;
  };

  auto el_id_up = [evaluator_el_sf, el_sf_era](const float &pt,
					     const float &sceta,
					     const int &is_mu,
					     const int &is_el,
					     const int &is_iso
					     ) {
    double sf = 1.;
    if (is_el == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_el_sf->evaluate({el_sf_era, "sfup", "Tight", sceta, pt});
    Logger::get("lepsf_electronIdSF")->debug("ID up - sf {}", sf);
    return sf;
  };

  auto el_id_down = [evaluator_el_sf, el_sf_era](const float &pt,
					     const float &sceta,
					     const int &is_mu,
					     const int &is_el,
					     const int &is_iso
					     ) {
    double sf = 1.;
    if (is_el == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_el_sf->evaluate({el_sf_era, "sfdown", "Tight", sceta, pt});
    Logger::get("lepsf_electronIdSF")->debug("ID down - sf {}", sf);
    return sf;
  };

  auto el_reco_nom = [evaluator_el_sf, el_sf_era](const float &pt,
					     const float &sceta,
					     const int &is_mu,
					     const int &is_el,
					     const int &is_iso
					     ) {
    double sf = 1.;
    if (is_el == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_el_sf->evaluate({el_sf_era, "sf", "RecoAbove20", sceta, pt});
    Logger::get("lepsf_electronRecoSF")->debug("RECO - sf {}", sf);
    return sf;
  };

  auto el_reco_up = [evaluator_el_sf, el_sf_era](const float &pt,
					     const float &sceta,
					     const int &is_mu,
					     const int &is_el,
					     const int &is_iso
					     ) {
    double sf = 1.;
    if (is_el == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_el_sf->evaluate({el_sf_era, "sfup", "RecoAbove20", sceta, pt});
    Logger::get("lepsf_electronRecoSF")->debug("RECO up - sf {}", sf);
    return sf;
  };

  auto el_reco_down = [evaluator_el_sf, el_sf_era](const float &pt,
					     const float &sceta,
					     const int &is_mu,
					     const int &is_el,
					     const int &is_iso
					     ) {
    double sf = 1.;
    if (is_el == 0 || is_iso != +1) return sf;
    if (pt >= 0.0)
      sf = evaluator_el_sf->evaluate({el_sf_era, "sfdown", "RecoAbove20", sceta, pt});
    Logger::get("lepsf_electronRecoSF")->debug("RECO down - sf {}", sf);
    return sf;
  };



  auto df10 = df9.Define(str_lep_sf_el_trigger_nom,
		       el_trigger_nom,
		       {str_lep_pt, str_lep_sceta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
		       );
  auto df11 = df10.Define(str_lep_sf_el_trigger_up,
			el_trigger_up,
			{str_lep_pt, str_lep_sceta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);
  auto df12 = df11.Define(str_lep_sf_el_trigger_down,
			el_trigger_down,
			{str_lep_pt, str_lep_sceta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);


  auto df13 = df12.Define(str_lep_sf_el_id_nom,
		       el_id_nom,
		       {str_lep_pt, str_lep_sceta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
		       );
  auto df14 = df13.Define(str_lep_sf_el_id_up,
			el_id_up,
			{str_lep_pt, str_lep_sceta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);
  auto df15 = df14.Define(str_lep_sf_el_id_down,
			el_id_down,
			{str_lep_pt, str_lep_sceta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);


  auto df16 = df15.Define(str_lep_sf_el_reco_nom,
		       el_reco_nom,
		       {str_lep_pt, str_lep_sceta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
		       );
  auto df17 = df16.Define(str_lep_sf_el_reco_up,
			el_reco_up,
			{str_lep_pt, str_lep_sceta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);
  auto df18 = df17.Define(str_lep_sf_el_reco_down,
			el_reco_down,
			{str_lep_pt, str_lep_sceta, str_lep_is_mu, str_lep_is_el, str_lep_is_iso}
			);


  return df18;



}


} // end namespace topreco

#endif /* GUARD_TOPRECO_H */
