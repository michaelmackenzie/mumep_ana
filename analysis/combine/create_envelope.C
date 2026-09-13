// Create a data-driven envelope
#ifndef __CONVANA_ANALYSIS_CREATEENVELOPE__
#define __CONVANA_ANALYSIS_CREATEENVELOPE__

#include <functional>

#include "../defaults.C"
#include "../tools/types.C"
#include "../tools/families.C"

// Defaults for envelope construction
bool useFrameChiSq_         = false; //use a roo plot frame to evaluate chi^2
bool useManualChisq_        = true ; //calcuate chi^2 values by hand with histograms
bool useDataBinErrors_      = false; //use data bin errors when calculating chi^2

bool use_exp_family_        = false; //exp(x)
bool use_power_family_      = false; //x^p
bool use_laurent_family_    = false; //Laurent
bool use_inv_poly_family_   = false; //1/(polynomial)
bool use_poly_family_       = true ; //polynomial
bool use_gaus_poly_family_  = true ; //Gaussian + polynomial
bool use_gaus_expo_family_  = false; //Gaussian + exp(x)
bool use_gaus_power_family_ = false; //Gaussian + x^p
bool use_generic_bernstein_ = false;
bool use_fast_bernstein_    = false;

bool force_fit_order_       = false; //force only the inclusion of fixed orders of each family
bool force_best_fit_        = false; //force the nominal PDF to be a specific function
bool enforce_ftest_         = true ; //once a function fails the F-test, stop adding functions
bool add_all_fits_          = false; //add all fits, even if they fail

bool test_single_function_  = false; //only pass 1 function into the PDF ensemble

double chisq_p_min_         = 0.001; //minimum p(chi^2, ndof) for a fit to be acceptable
double ftest_p_max_         = 0.20 ; //maximum p(F) for a higher order to be considered an improvement

RooAbsPdf* additional_bkg_  = nullptr; //if some background component is modeled separately
double     additional_bkg_norm_ = -1.; //fixed additional background yield, in the fit range
double     sideband_data_yield_ = -1.; //data yield in the fit range (-1: use the full dataset)

//-------------------------------------------------------------------------------------------------------------------
// Structure for storing fit result information
struct fit_res_t {
  bool accept_;
  double chi_sq_;
  int ndof_;
  int nparams_;
  int nentries_;
  double chi_sq_per_dof_;
  double p_;
  fit_res_t(bool accept, double chi_sq, int ndof, int nparams, int nentries) : accept_(accept), chi_sq_(chi_sq), ndof_(ndof),
                                                                               nparams_(nparams), nentries_(nentries) {
    chi_sq_per_dof_ = (ndof > 0) ? chi_sq / ndof : 0.;
    p_ = (ndof > 0.) ? TMath::Prob(chi_sq, ndof) : 0.;
  }
};

//----------------------------------------------------------------------------------------------------------------
//Perform the F-test selection for a higher order function: function_1 = nominal, function_2 = higher order
bool perform_f_test(double chisq_1, int ndof_1, double chisq_2, int ndof_2, double* p_val = nullptr) {
  if(p_val) *p_val = -1.; //undefined unless the test below is actually evaluated
  if(chisq_1 < 0. || ndof_1 <= 0) return true; //pick the higher order function
  if(chisq_2 < 0. || ndof_2 <= 0) return false; //stick with the original function
  if(ndof_1 == ndof_2) return chisq_1 > chisq_2; //if equal degrees of freedom, pick the better chi^2 model
  if(ndof_1 < ndof_2) {
    cout << __func__ << ": Error, higher order has higher N(dof)!\n";
    return false;
  }

  const int Mode = 1; //how to perform the F-test: 0: use F-distribution; 1: use delta chi^2 distribution
  double p_of_f = 1.;
  if(Mode == 0) { //use CDF for the F-distribution
    const double ftest = (chisq_1 / ndof_1) / (chisq_2 / ndof_2); //(chisq_1 - chisq_2) / (ndof_1 - ndof_2) / (chisq_2 / ndof_2);
    p_of_f = 1. - ROOT::Math::fdistribution_cdf(ftest, ndof_1, ndof_2);
  } else if(Mode == 1) { //use the difference of chi^2 p-value given N(dof) = Delta N(dof)
    const double ftest = chisq_1 - chisq_2;
    if(ftest < 0.) p_of_f = 1.; //lower order is a better fit
    else           p_of_f = ROOT::Math::chisquared_cdf_c(ftest, ndof_1 - ndof_2);
  } else {
    cout << __func__ << ": Unkown F-test mode " << Mode << endl;
    return false;
  }
  if(p_val) *p_val = p_of_f;
  return (p_of_f < ftest_p_max_); //select a higher order if confident it's better
}

//----------------------------------------------------------------------------------------------------------------
//Perform a p(chi^2, ndof) test to determine if the function has an acceptable agreement with data
bool perform_chisq_test(double chisq, int ndof, double* p_val = nullptr) {
  if(chisq < 0. || ndof <= 0) return false; //must have well-defined chi^2 and N(dof)
  const double p_chi_sq = TMath::Prob(chisq, ndof);
  if(p_val) *p_val = p_chi_sq;
  return (p_chi_sq > chisq_p_min_); //accept if the fit is compatible with the data
}

//----------------------------------------------------------------------------------------------------------------
//Get the chi-squared between a model and data from the model
double get_hist_chisquare(TH1* model, TH1* data, int* nbins = nullptr) {
  if(model->GetNbinsX() != data->GetNbinsX()) {
    cout << __func__ << ": Model and data don't have the same number of bins ("
         << model->GetNbinsX() << " vs " << data->GetNbinsX() << ")!\n";
    return -1.;
  }

  double chisq = 0.;
  for(int ibin = 1; ibin <= model->GetNbinsX(); ++ibin) {
    const double x_data = data ->GetBinCenter(ibin);
    const double x_pdf  = model->GetBinCenter(ibin);
    if(std::fabs(x_data - x_pdf) > 1.e-6) {
      cout << __func__ << ": Warning! Data center = " << x_data << " but PDF center = " << x_pdf << endl;
    }
    const double npdf  = model->GetBinContent(ibin);
    const double ndata = data->GetBinContent(ibin);
    const double error = data->GetBinError  (ibin);
    const double val  = ndata - npdf;
    const double sigma = (error > 0.) ? std::pow(val / error, 2) : 0.;
    chisq += sigma;
    if(verbose_ > 2) {
      cout << "Bin " << ibin << " (" << data->GetBinLowEdge(ibin) << " - "
           << data->GetBinLowEdge(ibin) + data->GetBinWidth(ibin) << "): Data = " << ndata << " PDF = " << npdf
           << " --> sigma = " << sigma << endl;
    }
  }

  if(verbose_ > 1) cout << __func__ << ": Total chi^2 = " << chisq << " / " << model->GetNbinsX() << " bins\n";
  if(nbins) *nbins = model->GetNbinsX();
  return chisq;
}

//----------------------------------------------------------------------------------------------------------------
//Get the chi-squared using a by-hand calculation
double get_manual_subrange_chisquare(RooRealVar& obs, RooAbsPdf* pdf, RooDataHist& data, const char* range = nullptr,
                                     const char* norm_range = nullptr, bool norm_skip = true, int* nbins = nullptr) {
  if(!pdf) return 0.;
  const int nhist_bins = data.numEntries();
  TH1* htmp_pdf  = pdf->createHistogram("htmp_chisq_pdf" , obs, RooFit::Binning(nhist_bins, obs.getMin(), obs.getMax()));
  TH1* htmp_data = data.createHistogram("htmp_chisq_data", obs, RooFit::Binning(nhist_bins, obs.getMin(), obs.getMax()));
  if(htmp_pdf->GetNbinsX() != htmp_data->GetNbinsX()) {
    cout << __func__ << ": PDF and data don't have the same number of bins ("
         << htmp_pdf->GetNbinsX() << " vs " << htmp_data->GetNbinsX() << ")!\n";
    delete htmp_pdf;
    delete htmp_data;
    return -1.;
  }

  //Create the PDF normalization by matching it to the data, skipping the region requested
  const double xmin_norm = (norm_range) ? obs.getMin(norm_range) : (norm_skip) ?  1. : obs.getMin();
  const double xmax_norm = (norm_range) ? obs.getMax(norm_range) : (norm_skip) ? -1. : obs.getMax();
  const int bin_norm_lo = (xmin_norm < xmax_norm) ? htmp_data->GetXaxis()->FindBin(xmin_norm) : (norm_skip) ? -1 : 1;
  const int bin_norm_hi = (xmin_norm < xmax_norm) ? htmp_data->GetXaxis()->FindBin(xmax_norm) : (norm_skip) ? -1 : htmp_data->GetNbinsX();
  double pdf_norm = 0.;
  double data_norm = 0.;
  for(int ibin = 1; ibin <= htmp_data->GetNbinsX(); ++ibin) {
    const double x_center = htmp_data->GetBinCenter(ibin);
    bool in_norm = (xmin_norm < xmax_norm) ? (x_center > xmin_norm && x_center < xmax_norm) : (ibin >= bin_norm_lo && ibin <= bin_norm_hi);
    if(norm_skip && in_norm) continue; //if given a range to skip
    if(!norm_skip && !in_norm) continue; //if given a range for normalizing
    const double npdf = htmp_pdf->GetBinContent(ibin)*htmp_pdf->GetBinWidth(ibin); //expected N(events)
    const double ndata = htmp_data->GetBinContent(ibin); //observed N(events)
    pdf_norm  += npdf;
    data_norm += ndata;
    if(verbose_ > 2) {
      cout << "Norm Bin " << ibin << " (" << htmp_data->GetBinLowEdge(ibin) << " - "
           << htmp_data->GetBinLowEdge(ibin) + htmp_data->GetBinWidth(ibin) << "): Data = " << ndata << " PDF = " << npdf  << endl;
    }
  }
  if(pdf_norm <= 0.) {
    cout << __func__ << ": PDF normalization is non-positive!\n";
    delete htmp_pdf;
    delete htmp_data;
    return -1.;
  }
  htmp_pdf->Scale(data_norm / pdf_norm);
  if(verbose_ > 2) {
    cout << __func__ << ": Scaling PDF:\n"
         << "#### Data norm = " << data_norm << endl
         << "#### PDF norm = " << pdf_norm << endl
         << "#### Scaling the PDF histogram by " << data_norm / pdf_norm << endl;
  }

  const double xmin = obs.getMin(range);
  const double xmax = obs.getMax(range);

  double chisq = 0.;
  int n_used_bins = 0;
  for(int ibin = 1; ibin <= htmp_data->GetNbinsX(); ++ibin) {
    const double x_center = htmp_data->GetBinCenter(ibin);
    if(x_center <= xmin || x_center >= xmax) continue;
    const double x_data = htmp_data->GetBinCenter(ibin);
    const double x_pdf  = htmp_pdf ->GetBinCenter(ibin);
    if(std::fabs(x_data - x_pdf) > 1.e-6) {
      cout << __func__ << ": Warning! Data center = " << x_data << " but PDF center = " << x_pdf << endl;
    }
    const double npdf  = htmp_pdf->GetBinContent(ibin)*htmp_pdf->GetBinWidth(ibin);
    const double ndata = htmp_data->GetBinContent(ibin);
    const double error = htmp_data->GetBinError  (ibin);
    const double val   = ndata - npdf;
    const double sigma = val*val / ((useDataBinErrors_) ? error*error : (npdf <= 0. ? 1.e-5 : npdf));
    chisq += sigma;
    ++n_used_bins;
    if(verbose_ > 2) {
      cout << "Bin " << ibin << " (" << htmp_data->GetBinLowEdge(ibin) << " - "
           << htmp_data->GetBinLowEdge(ibin) + htmp_data->GetBinWidth(ibin) << "): Data = " << ndata << " PDF = " << npdf
           << " --> sigma = " << sigma << endl;
    }
  }

  if(verbose_ > 1) cout << __func__ << ": Total chi^2 = " << chisq << " / " << n_used_bins << " bins\n";
  if(nbins) *nbins = n_used_bins;
  delete htmp_pdf;
  delete htmp_data;
  return chisq;
}

//Get the chi-squared using a RooChi2Var
double get_subrange_chisquare(RooRealVar& obs, RooAbsPdf* pdf, RooDataHist& data, const char* range) {
  RooAbsReal* chi = pdf ? pdf->createChi2(data, RooFit::Range(range)) : nullptr;
  if(!chi) return -1.;
  const double val = chi->getVal();
  delete chi;
  return val;
}

//Evaluate the chi-squared
double get_chi_squared(RooRealVar& obs, RooAbsPdf* pdf, RooDataHist& data, bool useSideBands, int* nbins = nullptr) {
  if(useManualChisq_) {
    if(useSideBands || blind_data_) {
      const char* blind_region = "BlindRegion"; //skip the blinding region if fitting the sidebands
      int nbin_running = 0;
      double chi_sq = get_manual_subrange_chisquare(obs, pdf, data, "LowSideband", blind_region, true, &nbin_running);
      if(nbins) *nbins = nbin_running;
      chi_sq += get_manual_subrange_chisquare(obs, pdf, data, "HighSideband", blind_region, true, &nbin_running);
      if(nbins) *nbins = *nbins + nbin_running;
      return chi_sq;
    }
    return get_manual_subrange_chisquare(obs, pdf, data, "full");
  }
  if(useFrameChiSq_ || TString(pdf->GetTitle()).Contains("Exponential")) {
    auto xframe = obs.frame();
    data.plotOn(xframe);
    pdf->plotOn(xframe);
    const double chi_sq = xframe->chiSquare() * data.numEntries(); //returns chi squared / entries
    delete xframe;
    return chi_sq;
  }
  if(useSideBands) {
    double chi_sq = get_subrange_chisquare(obs, pdf, data, "LowSideband");
    chi_sq += get_subrange_chisquare(obs, pdf, data, "HighSideband");
    return chi_sq;
  }
  return get_subrange_chisquare(obs, pdf, data, "full");
}

//-----------------------------------------------------------------------------------------------------------------------------------
// Helper function to combine the data-driven fit with additional backgrounds
// The yields are expressed in the range the fit normalizes over, so additional_bkg_norm_ and
// sideband_data_yield_ must already be restricted to that range by the caller
RooAbsPdf* wrap_pdf(RooAbsPdf* pdf, double norm) {
    const double fixed_norm = (additional_bkg_ && additional_bkg_norm_ > 0.) ? additional_bkg_norm_ : 0.;
    const double target = (sideband_data_yield_ > 0.) ? sideband_data_yield_ : norm;
    const double float_norm = std::max(1.e-3, target - fixed_norm);
    RooRealVar* pdf_norm = new RooRealVar(Form("%s_pdf_norm", pdf->GetName()), Form("%s_pdf_norm", pdf->GetName()),
                                          float_norm, 0., std::max(2.*target, 10.));
    if(!additional_bkg_ || fixed_norm <= 0.) {
      return new RooAddPdf(Form("%s_total", pdf->GetName()), pdf->GetTitle(), RooArgList(*pdf), RooArgList(*pdf_norm));
    }
    // The explicitly modeled background is held at its expectation while the envelope floats
    RooRealVar* add_norm = new RooRealVar(Form("%s_add_norm", pdf->GetName()), Form("%s_add_norm", pdf->GetName()), fixed_norm);
    add_norm->setConstant(true);
    return new RooAddPdf(Form("%s_total", pdf->GetName()), pdf->GetTitle(),
                         RooArgList(*pdf, *additional_bkg_), RooArgList(*pdf_norm, *add_norm));
}

//-----------------------------------------------------------------------------------------------------------------------------------
// Fit a PDF to the dataset and determine whether or not to accept it
fit_res_t fit_pdf_to_data(RooAbsPdf* pdf, RooDataHist& data, RooRealVar& obs, bool useSideBands, int verbose = 0) {
  if(useSideBands)
    // Range() also fixes the normalization range, so no separate NormRange is needed (nor accepted)
    pdf->fitTo(data, RooFit::PrintLevel(-1 + max(0, verbose-2)), RooFit::Warnings(0), RooFit::PrintEvalErrors(-1),
               RooFit::Range("LowSideband,HighSideband"));
  else
    pdf->fitTo(data, RooFit::PrintLevel(-1 + max(0, verbose-2)), RooFit::Warnings(0), RooFit::PrintEvalErrors(-1));
  int nentries = data.numEntries(); //non-constant, as this is updated in the chi^2 calculation based on region used
  const double chi_sq = get_chi_squared(obs, pdf, data, useSideBands, &nentries);
  const int nparams = count_pdf_params(pdf, obs);
  const int dof = (nentries - nparams - 1);  //Assume normalization is not counted in the parameter list
  double p_chi_sq;
  const bool accept = perform_chisq_test(chi_sq, dof, &p_chi_sq);
  return fit_res_t(accept,chi_sq, dof, nparams, nentries);
}

///////////////////////////////////////////
// Add functions to the envelope
///////////////////////////////////////////

//----------------------------------------------------------------------------------------------------------------
// Record of a function accepted into the envelope, kept for the summary plot and printouts
struct env_fit_t {
  RooAbsPdf* pdf_       = nullptr; //shape added to the envelope
  RooAbsPdf* total_pdf_ = nullptr; //fitted shape plus the fixed explicit background, as used in the fit
  TString    label_     = ""     ;
  double     chi_sq_    = -1.    ;
  int        ndof_      =  0     ;
  int        nparams_   =  0     ;
};

//Function that builds a given order of a family
typedef std::function<RooAbsPdf*(RooRealVar&, int, TString)> pdf_creator_t;

//Fit each order of a family and add the passing ones
//  create       : builds the order being tested
//  min/max_order: orders to scan
//  forced_order : the only order kept when force_fit_order_ is set
std::pair<int,double> add_family(RooDataHist& data, RooRealVar& obs, RooArgList& list,
                                 std::vector<env_fit_t>& fits, bool useSideBands, TString name, int verbose,
                                 pdf_creator_t create, TString label,
                                 const int min_order, const int max_order, const int forced_order) {
  double min_chi = 1.e10;
  int min_index = -1;
  double prev_chi = -1.;
  int prev_ndof = -1;
  int prev_order = -1;
  bool has_prev = false;
  int n_tested = 0, n_accepted = 0;
  if(verbose > 0) {
    TString policy = (add_all_fits_)    ? "keeping every order (add_all_fits_)"
                   : (force_fit_order_) ? TString::Format("keeping only order %i (force_fit_order_)", forced_order)
                                        : TString::Format("keeping orders with p(chi^2) > %.3g", chisq_p_min_);
    if(!add_all_fits_ && !force_fit_order_ && enforce_ftest_)
      policy += TString::Format(", stopping once p(F) > %.3g", ftest_p_max_);
    printf("===============================================================================================\n");
    printf("%s: ===== %s family, orders %i - %i: %s =====\n",
           __func__, label.Data(), min_order, max_order, policy.Data());
    printf("===============================================================================================\n");
  }
  for(int order = min_order; order <= max_order; ++order) {
    RooAbsPdf* basePdf = create(obs, order, name);
    if(!basePdf) {
      if(verbose > 0) printf("%s:   %s order %i: not available\n", __func__, label.Data(), order);
      continue;
    }
    ++n_tested;
    RooAbsPdf* pdf = wrap_pdf(basePdf, data.sumEntries());
    auto res = fit_pdf_to_data(pdf, data, obs, useSideBands, verbose);
    bool accept = add_all_fits_ || (res.accept_ && !force_fit_order_) || (force_fit_order_ && order == forced_order);

    //Record why this order was or was not kept, so a scan can be followed from the log alone
    TString reason;
    if(add_all_fits_)         reason = "add_all_fits_ is set";
    else if(force_fit_order_) reason = (order == forced_order)
                                     ? "the forced order"
                                     : TString::Format("not the forced order (%i)", forced_order);
    else                      reason = (res.accept_)
                                     ? TString::Format("p(chi^2) = %.3g > %.3g", res.p_, chisq_p_min_)
                                     : TString::Format("p(chi^2) = %.3g < %.3g", res.p_, chisq_p_min_);

    bool stop_orders = false;
    TString ftest_str;
    //the F-test picks an order, so it only applies when the order isn't already forced
    if(enforce_ftest_ && !add_all_fits_ && !force_fit_order_ && has_prev) {
      double p_ftest = -1.;
      const bool pass_ftest = perform_f_test(prev_chi, prev_ndof, res.chi_sq_, res.ndof_, &p_ftest);
      ftest_str = (p_ftest >= 0.) ? TString::Format(" F-test vs order %i: p = %8.2e,", prev_order, p_ftest)
                                  : TString::Format(" F-test vs order %i: n/a,", prev_order);
      if(!pass_ftest) {
        accept = false;
        stop_orders = true;
        reason = TString::Format("no significant improvement over order %i", prev_order);
      }
    }

    if(verbose > 0) {
    printf("===============================================================================================\n");
      printf("%s:   %s order %i: %2i par, %3i bins, chi^2/dof = %8.2f / %3i = %6.3f, p = %8.2e,%s %s (%s)\n",
             __func__, label.Data(), order, res.nparams_, res.nentries_,
             res.chi_sq_, res.ndof_, res.chi_sq_per_dof_, res.p_, ftest_str.Data(),
             (accept) ? "ACCEPTED" : "rejected", reason.Data());
      pdf->Print("tree");
      printf("===============================================================================================\n");
    }

    if(accept) {
      ++n_accepted;
      list.add(*basePdf);
      env_fit_t fit;
      fit.pdf_       = basePdf;
      fit.total_pdf_ = pdf;
      fit.label_     = Form("%s order %i", label.Data(), order);
      fit.chi_sq_    = res.chi_sq_;
      fit.ndof_      = res.ndof_;
      fit.nparams_   = res.nparams_;
      fits.push_back(fit);
    } else {
      delete pdf;
    }
    has_prev = true;
    prev_chi = res.chi_sq_;
    prev_ndof = res.ndof_;
    prev_order = order;
    if(res.chi_sq_ < min_chi) {min_chi = res.chi_sq_; min_index = list.getSize() - 1;}
    if(stop_orders) {
      printf("===============================================================================================\n");
      if(verbose > 0) printf("%s:   %s: stopping the scan at order %i (higher orders not tested)\n",
                             __func__, label.Data(), order);
      printf("===============================================================================================\n");
      break;
    }
  }
  if(verbose > 0) {
    printf("===============================================================================================\n");
    printf("%s:   %s: %i / %i tested order(s) added to the envelope\n",
           __func__, label.Data(), n_accepted, n_tested);
    printf("===============================================================================================\n");
  }
  return std::pair<int, double>(min_index, min_chi);
}

//----------------------------------------------------------------------------------------------------------------
// Plot the accepted envelope functions and their agreement with the data they were fit to
// The comparison is done by hand, consistently with the chi^2 evaluation, so the blinded region is
// left out of the normalization and of the pulls
int plot_envelope_functions(RooRealVar& obs, RooDataHist& data, const std::vector<env_fit_t>& fits,
                            const TString figdir, const TString name,
                            const double blind_min, const double blind_max) {
  if(fits.empty()) {
    cout << __func__ << ": No accepted functions to plot for " << name.Data() << endl;
    return 1;
  }
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir.Data(), figdir.Data()));

  const int nbins   = data.numEntries();
  const double xmin = obs.getMin(), xmax = obs.getMax();
  auto in_sideband  = [&](const double x) { return !(x > blind_min && x < blind_max); };

  TH1* h_data = (TH1*) data.createHistogram(Form("%s_plot_data", name.Data()), obs,
                                            RooFit::Binning(nbins, xmin, xmax));
  double data_sideband = 0.;
  for(int ibin = 1; ibin <= h_data->GetNbinsX(); ++ibin) {
    const double n = h_data->GetBinContent(ibin);
    if(h_data->GetBinError(ibin) <= 0.) h_data->SetBinError(ibin, std::sqrt(std::max(0., n)));
    if(in_sideband(h_data->GetBinCenter(ibin))) data_sideband += n;
  }

  //Build each model as expected events per bin, matched to the data in the sidebands
  const int colors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7, kCyan+2, kViolet+1, kSpring-6, kAzure+7, kPink+7};
  const int ncolors = sizeof(colors)/sizeof(colors[0]);
  //families often land on top of each other, so vary the line style as well as the color
  const int styles[] = {kSolid, kDashed, kDotted, kDashDotted};
  const int nstyles = sizeof(styles)/sizeof(styles[0]);
  std::vector<TH1*> models, pulls;
  for(size_t ifit = 0; ifit < fits.size(); ++ifit) {
    TH1* h = (TH1*) fits[ifit].total_pdf_->createHistogram(Form("%s_plot_model_%i", name.Data(), (int) ifit), obs,
                                                          RooFit::Binning(nbins, xmin, xmax));
    double model_sideband = 0.;
    for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
      h->SetBinContent(ibin, h->GetBinContent(ibin)*h->GetBinWidth(ibin)); //density --> expected events
      h->SetBinError(ibin, 0.);
      if(in_sideband(h->GetBinCenter(ibin))) model_sideband += h->GetBinContent(ibin);
    }
    if(model_sideband > 0.) h->Scale(data_sideband / model_sideband);

    TH1* h_pull = (TH1*) h->Clone(Form("%s_plot_pull_%i", name.Data(), (int) ifit));
    h_pull->Reset();
    for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
      const double model = h->GetBinContent(ibin);
      const double ndata = h_data->GetBinContent(ibin);
      const double err   = std::sqrt(std::max(1.e-9, model));
      h_pull->SetBinContent(ibin, (ndata - model)/err);
    }

    const int color = colors[ifit % ncolors];
    for(TH1* hh : {h, h_pull}) {
      hh->SetLineColor(color);
      hh->SetMarkerColor(color);
      hh->SetLineWidth(2);
      hh->SetStats(0);
    }
    h->SetLineStyle(styles[(ifit/ncolors) % nstyles]);
    if(fits.size() > 1) h->SetLineStyle(styles[ifit % nstyles]);
    h_pull->SetMarkerStyle(20);
    h_pull->SetMarkerSize(0.7);
    models.push_back(h);
    pulls.push_back(h_pull);
  }

  //Canvas with the functions on top and the pulls below
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas(Form("c_%s_functions", name.Data()), "Envelope functions", 1200, 1000);
  TPad* pad1 = new TPad("pad1_env", "pad1_env", 0., 0.3, 1., 1.);
  TPad* pad2 = new TPad("pad2_env", "pad2_env", 0., 0. , 1., 0.3);
  pad1->SetBottomMargin(0.03); pad1->SetLeftMargin(0.12); pad1->SetRightMargin(0.04); pad1->SetTopMargin(0.08);
  pad2->SetTopMargin(0.04);    pad2->SetLeftMargin(0.12); pad2->SetRightMargin(0.04); pad2->SetBottomMargin(0.32);
  pad1->Draw(); pad2->Draw();

  //Mark the blinded region left out of the fits
  auto draw_blind_box = [&](const double ymin, const double ymax) {
    TBox* box = new TBox(blind_min, ymin, blind_max, ymax);
    box->SetFillColorAlpha(kGray+1, 0.25);
    box->SetLineColor(kGray+2);
    box->SetLineStyle(kDashed);
    box->Draw("same");
  };

  h_data->SetTitle("");
  h_data->SetLineColor(kBlack);
  h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.8);
  h_data->GetYaxis()->SetTitle(Form("Events / %.2g MeV/c", h_data->GetBinWidth(1)));
  h_data->GetYaxis()->SetTitleSize(0.055);
  h_data->GetYaxis()->SetTitleOffset(1.05);
  h_data->GetYaxis()->SetLabelSize(0.045);
  h_data->GetXaxis()->SetLabelSize(0.);
  const double ymax_data = h_data->GetMaximum();
  double ymin_log = 1.e30;
  for(int ibin = 1; ibin <= h_data->GetNbinsX(); ++ibin) {
    const double n = h_data->GetBinContent(ibin);
    if(n > 0. && n < ymin_log) ymin_log = n;
  }
  if(ymin_log > 1.e29) ymin_log = 1.;

  TLegend* leg = new TLegend(0.14, 0.63, 0.96, 0.92);
  leg->SetNColumns(2);
  leg->SetLineWidth(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.03);
  leg->AddEntry(h_data, "Data (sidebands fit)", "PE");
  for(size_t ifit = 0; ifit < fits.size(); ++ifit) {
    leg->AddEntry(models[ifit], Form("%s: #chi^{2}/dof = %.1f/%i", fits[ifit].label_.Data(),
                                     fits[ifit].chi_sq_, fits[ifit].ndof_), "L");
  }

  //Redraw the top pad from scratch for each scale, so the axis limits are set before it is painted
  auto draw_top = [&](const bool logy) {
    pad1->cd();
    pad1->Clear();
    pad1->SetLogy(logy);
    const double ymin = (logy) ? 0.5*ymin_log      : 0.;
    const double ymax = (logy) ? 100.*ymax_data    : 2.2*ymax_data;
    h_data->SetMinimum(ymin);
    h_data->SetMaximum(ymax);
    h_data->Draw("E1");
    draw_blind_box(std::max(ymin, 0.), ymax);
    for(auto h : models) h->Draw("hist same");
    h_data->Draw("E1 same");
    leg->Draw();
    gPad->RedrawAxis();
  };
  draw_top(false);

  pad2->cd();
  TH1* h_axis = (TH1*) pulls[0]->Clone(Form("%s_plot_pull_axis", name.Data()));
  h_axis->Reset();
  h_axis->SetTitle("");
  h_axis->GetYaxis()->SetTitle("#frac{Data - Fit}{#sigma}");
  h_axis->GetYaxis()->SetRangeUser(-4.5, 4.5);
  h_axis->GetYaxis()->SetNdivisions(505);
  h_axis->GetYaxis()->SetTitleSize(0.11);
  h_axis->GetYaxis()->SetTitleOffset(0.45);
  h_axis->GetYaxis()->SetLabelSize(0.10);
  h_axis->GetXaxis()->SetTitle("Momentum (MeV/c)");
  h_axis->GetXaxis()->SetTitleSize(0.13);
  h_axis->GetXaxis()->SetTitleOffset(1.05);
  h_axis->GetXaxis()->SetLabelSize(0.11);
  h_axis->Draw();
  draw_blind_box(-4.5, 4.5);
  for(int sigma : {-2, -1, 1, 2}) {
    TLine* line = new TLine(xmin, sigma, xmax, sigma);
    line->SetLineStyle(kDotted);
    line->SetLineColor(kGray+2);
    line->Draw("same");
  }
  TLine* zero = new TLine(xmin, 0., xmax, 0.);
  zero->SetLineColor(kBlack);
  zero->Draw("same");
  for(auto h : pulls) h->Draw("P same");
  gPad->RedrawAxis();

  const TString out_base = Form("%s/%s_functions", figdir.Data(), name.Data());
  c->SaveAs(Form("%s.png", out_base.Data()));
  draw_top(true);
  c->SaveAs(Form("%s_log.png", out_base.Data()));
  delete c;
  return 0;
}

//----------------------------------------------------------------------------------------------------------------
// Main function to create the envelope
RooMultiPdf* create_envelope(RooRealVar& obs, RooCategory& cat, RooDataHist& data, bool useSideBands, TString name,
                             std::vector<env_fit_t>& fits, const int verbose = 0) {
  RooArgList pdf_list;
  std::pair<int, double> result(-1, 1.e10);
  double chi_min = 1.e10;
  int best_index = -1;
  fits.clear();

  auto consider = [&](pdf_creator_t create, TString label, int min_order, int max_order, int forced_order) {
    result = add_family(data, obs, pdf_list, fits, useSideBands, name, verbose,
                        create, label, min_order, max_order, forced_order);
    if(result.second < chi_min) { chi_min = result.second; best_index = result.first; }
  };

  if(use_poly_family_) {
    // consider([](RooRealVar& o, int n, TString s) -> RooAbsPdf* {
    //            if(use_generic_bernstein_) return create_generic_bernstein(o, n, s);
    //            if(use_fast_bernstein_)    return create_fast_bernstein   (o, n, s);
    //            return create_bernstein(o, n, s);
    //          }, "Bernstein", 1, 4, 3);
    consider([](RooRealVar& o, int n, TString s) -> RooAbsPdf* { return create_chebychev(o, n, s); },
             "Chebychev", 1, 4, 3);
  }
  if(use_exp_family_)
    consider([](RooRealVar& o, int n, TString s) -> RooAbsPdf* { return create_exponential(o, n, s); },
             "Exponential", 1, 3, 2);
  if(use_power_family_)
    consider([](RooRealVar& o, int n, TString s) -> RooAbsPdf* { return create_powerlaw(o, n, s); },
             "Power law", 1, 3, 2);
  if(use_laurent_family_)
    consider([](RooRealVar& o, int n, TString s) -> RooAbsPdf* { return create_laurent(o, n, s); },
             "Laurent", 1, 6, 1);
  if(use_inv_poly_family_)
    consider([](RooRealVar& o, int n, TString s) -> RooAbsPdf* { return create_inv_polynomial(o, n, s); },
             "Inverse poly", 1, 1, 1);
  if(use_gaus_poly_family_)
    consider([](RooRealVar& o, int n, TString s) -> RooAbsPdf* { return create_gaus_poly_pdf(o, n, s); },
             "Gaussian poly", -1, 3, 1);
  if(use_gaus_expo_family_)
    consider([](RooRealVar& o, int n, TString s) -> RooAbsPdf* { return create_gaus_expo_pdf(o, n, s); },
             "Gaussian expo", 1, 2, 1);
  if(use_gaus_power_family_)
    consider([](RooRealVar& o, int n, TString s) -> RooAbsPdf* { return create_gaus_power_pdf(o, n, s); },
             "Gaussian power", 1, 2, 1);

  if(pdf_list.getSize() <= 0) {
    cout << __func__ << ": No PDFs were added to the envelope for " << name.Data() << endl;
    return nullptr;
  }

  if(test_single_function_) {
    if(best_index < 0 || best_index >= pdf_list.getSize()) best_index = 0;
    RooAbsPdf* single_pdf = (RooAbsPdf*) pdf_list.at(best_index);
    RooArgList single_list;
    single_list.add(*single_pdf);
    if(fits.size() > (size_t) best_index) fits = std::vector<env_fit_t>{fits[best_index]};
    return new RooMultiPdf(Form("%s_pdf", name.Data()), Form("%s PDF", name.Data()), cat, single_list);
  }

  return new RooMultiPdf(Form("%s_pdf", name.Data()), Form("%s PDF", name.Data()), cat, pdf_list);
}

///////////////////////////////////////////
// Envelope model construction
///////////////////////////////////////////

//----------------------------------------------------------------------------------------------------------------
// Everything the model builder needs to write the envelope into a workspace/data card
struct envelope_t {
  RooMultiPdf* pdf_           = nullptr; // discretely profiled envelope of background functions
  RooCategory* cat_           = nullptr; // index over the envelope functions, profiled by Combine
  RooRealVar*  norm_          = nullptr; // freely floating envelope yield, named <pdf>_norm
  RooAbsPdf*   explicit_pdf_  = nullptr; // frozen sum of the explicitly modeled backgrounds
  double       explicit_rate_ = 0.;      // expected explicit background yield, full range
  double       rate_          = 0.;      // envelope yield extrapolated to the full range
  TString      name_          = "";
  TString      title_         = "Background envelope";
  std::vector<env_fit_t> fits_;          // the functions accepted into the envelope
};

//----------------------------------------------------------------------------------------------------------------
// Define the ranges the envelope fits use: the two sidebands around the blinded signal region
void set_envelope_ranges(RooRealVar& obs, const double blind_min, const double blind_max) {
  const double xmin = obs.getMin(), xmax = obs.getMax();
  obs.setRange("full"        , xmin     , xmax     );
  obs.setRange("LowSideband" , xmin     , blind_min);
  obs.setRange("HighSideband", blind_max, xmax     );
  obs.setRange("BlindRegion" , blind_min, blind_max);
}

//----------------------------------------------------------------------------------------------------------------
// Fraction of a PDF falling into the two sidebands
// Evaluated from a fine histogram rather than a multi-range RooFit integral, which normalizes
// each range to itself instead of to the full range
double pdf_sideband_fraction(RooRealVar& obs, RooAbsPdf* pdf, const double blind_min, const double blind_max) {
  if(!pdf) return 0.;
  std::unique_ptr<TH1> h((TH1*) pdf->createHistogram(Form("h_sideband_frac_%s", pdf->GetName()), obs,
                                                     RooFit::Binning(2000, obs.getMin(), obs.getMax())));
  if(!h) return 0.;
  const double total = h->Integral();
  if(total <= 0.) return 0.;
  double sidebands = 0.;
  for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
    const double x = h->GetBinCenter(ibin);
    if(x > blind_min && x < blind_max) continue;
    sidebands += h->GetBinContent(ibin);
  }
  return sidebands / total;
}

//----------------------------------------------------------------------------------------------------------------
// Data yield outside of the blinded region
double data_sideband_yield(RooRealVar& obs, RooDataHist& data, const double blind_min, const double blind_max) {
  double total = 0.;
  for(int entry = 0; entry < data.numEntries(); ++entry) {
    const RooArgSet* vars = data.get(entry);
    auto var = (RooRealVar*) vars->find(obs.GetName());
    if(!var) continue;
    const double x = var->getVal();
    if(x > blind_min && x < blind_max) continue;
    total += data.weight();
  }
  return total;
}

//----------------------------------------------------------------------------------------------------------------
// Sum the explicitly modeled backgrounds into a single PDF with every parameter frozen, including
// the expected rates, so the sideband fits only float the envelope on top of them
RooAbsPdf* build_explicit_background_pdf(RooRealVar& obs,
                                         const std::vector<pdf_info>& models,
                                         const TString name,
                                         double& total_rate) {
  total_rate = 0.;
  RooArgList pdfs, rates;
  for(auto& info : models) {
    if(!info.pdf_) {
      cout << __func__ << ": Explicit background " << info.name_.Data() << " has no PDF!\n";
      continue;
    }
    // PDFs read back from a workspace carry their own copy of the observable, which has none of
    // the sideband ranges defined on it, so point them at the observable the model is built with
    info.pdf_->recursiveRedirectServers(RooArgSet(obs));
    // Freeze the shape
    std::unique_ptr<RooArgSet> params(info.pdf_->getParameters(RooArgSet(obs)));
    for(auto param : *params) {
      auto var = dynamic_cast<RooRealVar*>(param);
      if(var) var->setConstant(true);
    }
    // Freeze the expected rate
    auto rate = new RooRealVar(Form("%s_%s_rate", name.Data(), info.name_.Data()),
                               Form("%s expected yield", info.title_.Data()), info.rate_);
    rate->setConstant(true);
    pdfs.add(*info.pdf_);
    rates.add(*rate);
    total_rate += info.rate_;
  }
  if(pdfs.getSize() == 0) return nullptr;
  return new RooAddPdf(Form("%s_explicit_pdf", name.Data()), "Explicitly modeled background", pdfs, rates);
}

//----------------------------------------------------------------------------------------------------------------
// Build the data-driven background envelope from the data sidebands
//   explicit_models: backgrounds modeled outside of the envelope, held fixed during the fits
envelope_t build_background_envelope(RooRealVar& obs,
                                     RooDataHist& data,
                                     const std::vector<pdf_info>& explicit_models,
                                     const TString name,
                                     const double blind_min,
                                     const double blind_max,
                                     const int verbose = 0,
                                     const TString figdir = "") {
  envelope_t env;
  env.name_ = name;

  set_envelope_ranges(obs, blind_min, blind_max);

  env.explicit_pdf_ = build_explicit_background_pdf(obs, explicit_models, name, env.explicit_rate_);

  // The sideband fits normalize over the sidebands only, so the fixed explicit yield and the
  // envelope starting point have to be expressed in that range as well
  const double explicit_frac = pdf_sideband_fraction(obs, env.explicit_pdf_, blind_min, blind_max);
  const double data_yield    = data_sideband_yield(obs, data, blind_min, blind_max);

  additional_bkg_      = env.explicit_pdf_;
  additional_bkg_norm_ = env.explicit_rate_ * explicit_frac;
  sideband_data_yield_ = data_yield;

  if(verbose > 0) {
    printf("%s: Fitting the envelope to the sidebands of %s (blinding %.2f - %.2f)\n",
           __func__, name.Data(), blind_min, blind_max);
    printf("%s: Sideband data yield = %.2f, fixed explicit background = %.2f (%.1f%% of %.2f)\n",
           __func__, data_yield, additional_bkg_norm_, 100.*explicit_frac, env.explicit_rate_);
  }

  env.cat_ = new RooCategory(Form("%s_cat", name.Data()), Form("%s function index", name.Data()));
  env.pdf_ = create_envelope(obs, *env.cat_, data, true, name, env.fits_, verbose);

  // Leave the shared state clean for any later calls
  additional_bkg_      = nullptr;
  additional_bkg_norm_ = -1.;
  sideband_data_yield_ = -1.;

  if(!env.pdf_) {
    cout << __func__ << ": Failed to build the envelope for " << name.Data() << endl;
    return env;
  }

  // Extrapolate the fitted sideband yield back to the full range for the starting normalization
  const double env_sideband = std::max(0., data_yield - env.explicit_rate_*explicit_frac);
  const double env_frac = pdf_sideband_fraction(obs, env.pdf_->getCurrentPdf(), blind_min, blind_max);
  env.rate_ = (env_frac > 0.) ? env_sideband / env_frac : env_sideband;

  // Combine picks up <pdf name>_norm as the normalization of a parametric shape, so this is the
  // freely floating, data-driven background yield
  env.norm_ = new RooRealVar(Form("%s_norm", env.pdf_->GetName()), Form("%s yield", env.title_.Data()),
                             env.rate_, 0., std::max(10.*env.rate_, 10.));
  env.norm_->setConstant(false);

  if(verbose > 0) {
    printf("===============================================================================================\n");
    printf("%s: Envelope has %i function(s), sideband fraction = %.3f --> starting yield = %.2f\n",
           __func__, env.pdf_->getNumPdfs(), env_frac, env.rate_);
    for(int ipdf = 0; ipdf < env.pdf_->getNumPdfs(); ++ipdf) {
      const bool has_fit = ipdf < (int) env.fits_.size();
      printf("%s:   [%i] %-30s %s\n", __func__, ipdf, env.pdf_->getPdf(ipdf)->GetName(),
             (has_fit) ? Form("%-22s chi^2/dof = %.1f / %i (%i par)", env.fits_[ipdf].label_.Data(),
                              env.fits_[ipdf].chi_sq_, env.fits_[ipdf].ndof_, env.fits_[ipdf].nparams_)
                       : env.pdf_->getPdf(ipdf)->GetTitle());
    }
    printf("===============================================================================================\n");
  }

  // Show the accepted functions and how well each describes the data they were fit to
  if(figdir != "") plot_envelope_functions(obs, data, env.fits_, figdir, name, blind_min, blind_max);

  return env;
}

#endif
