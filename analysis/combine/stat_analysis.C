// Perform statistical analysis using a RooFit workspace built by build_model.C

#include "../tools/types.C"
#include "../physics.C"

using namespace RooStats;

bool include_systematics_ = false; // whether or not to consider systematics

//---------------------------------------------------------------------------
// Retrieve a variable from the workspace (read-only copy)
//---------------------------------------------------------------------------
RooRealVar* get_var(RooWorkspace* ws, const char* name) {
  RooRealVar* v = (RooRealVar*) ws->var(name);
  if(!v) std::cout << "stat_analysis: WARNING – variable '" << name << "' not in workspace\n";
  return v;
}

//---------------------------------------------------------------------------
// Tetrieve a PDF from the workspace
//---------------------------------------------------------------------------
RooAbsPdf* get_pdf(RooWorkspace* ws, const char* name) {
  RooAbsPdf* p = (RooAbsPdf*) ws->pdf(name);
  if(!p) std::cout << "stat_analysis: WARNING – PDF '" << name << "' not in workspace\n";
  return p;
}

//---------------------------------------------------------------------------
// Retrieve a dataset from the workspace
//---------------------------------------------------------------------------
RooAbsData* get_data_ws(RooWorkspace* ws, const char* name = "data_obs") {
  RooAbsData* d = ws->data(name);
  if(!d) std::cout << "stat_analysis: WARNING – dataset '" << name << "' not in workspace\n";
  return d;
}

//---------------------------------------------------------------------------
// Reset nuisances
//---------------------------------------------------------------------------
void reset_nuisances(RooArgSet& nuisanceParams) {
  // Reset nuisances
  TIter nit(nuisanceParams.createIterator());
  TObject* nobj;
  while((nobj = nit())) {
    if(auto* rv = dynamic_cast<RooRealVar*>(nobj)) {
      if(TString(rv->GetName()).BeginsWith("theta")) rv->setVal(0.); // thetas
      else                                           rv->setVal(1.); // kappas
      rv->setConstant(!include_systematics_);
    }
  }
}

//---------------------------------------------------------------------------
// Systematic uncertainty descriptors
//---------------------------------------------------------------------------
struct ShapeSysInfo_t {
  TString name_;
  std::map<TString,TString> pdf_up_;
  std::map<TString,TString> pdf_down_;
};

struct RateSysInfo_t {
  TString name_;
  std::map<TString,double> scales_;
};

//---------------------------------------------------------------------------
// Build the full extended PDF with systematic nuisance parameters
// PDF = p(data | model) * p(theta_0) * p(theta_1)...
//---------------------------------------------------------------------------
RooProdPdf* build_full_model(
    RooWorkspace* ws,
    const TString process,
    const int selection,
    RooRealVar* mu,                  // signal strength
    std::vector<ShapeSysInfo_t>& sys_infos,
    RooArgSet& nuisanceParams,
    RooArgSet& constraintPdfs,
    std::vector<RooAbsPdf*>& owned_pdfs,   // keep PDFs alive
    std::vector<RooAbsReal*>& owned_reals  // keep RooFormulaVar etc. alive
) {
  const bool is_mumem = process == "mumem";
  TString sel_str = Form("_%i_", selection);

  //-----------------------------------------------------------------
  // Collect component names from workspace
  //-----------------------------------------------------------------
  std::vector<TString> components;
  components.push_back("signal");
  if(is_mumem) components.push_back("dio");
  components.push_back("cosmic");
  components.push_back("rpc_ext");
  components.push_back("rpc_int");
  if(!use_evtana_) components.push_back("pbar");
  if(!is_mumem) {
    components.push_back("rmc_ext");
    components.push_back("rmc_int");
  }

  //-----------------------------------------------------------------
  // Scan workspace for systematic PDF variants
  //-----------------------------------------------------------------
  // Pattern: <process>_<sel>_<comp>_pdf_<SysName>Up / Down
  // Example: mumem_20_signal_pdfScaleUp
  std::map<TString, ShapeSysInfo_t> sys_map;
  {
    RooArgSet all_ws_pdfs = ws->allPdfs();
    TIter it(all_ws_pdfs.createIterator());
    TObject* obj;
    while((obj = it())) {
      TString pname = obj->GetName();
      for(auto& comp : components) {
        TString base = Form("%s_%i_%s_pdf", process.Data(), selection, comp.Data());
        if(!pname.BeginsWith(base)) continue;
        TString suffix = pname(base.Length(), pname.Length()-base.Length());
        if(suffix == "" || suffix == "_nom") continue;
        // suffix is like "ScaleUp" or "ScaleDown"
        bool is_up   = suffix.EndsWith("Up");
        bool is_down = suffix.EndsWith("Down");
        if(!is_up && !is_down) continue;
        TString sys_name = suffix(0, suffix.Length() - (is_up ? 2 : 4));
        if(!sys_map.count(sys_name)) { ShapeSysInfo_t s; s.name_ = sys_name; sys_map[sys_name] = s; }
        if(is_up)   sys_map[sys_name].pdf_up_  [comp] = pname;
        else        sys_map[sys_name].pdf_down_ [comp] = pname;
      }
    }
  }
  for(auto& kv : sys_map) sys_infos.push_back(kv.second);

  // Scan for rate uncertainties
  // Pattern: <process>_<sel>_<comp>_RateSys_<name>
  std::map<TString, RateSysInfo_t> rate_sys_map;
  {
    RooArgSet all_ws_vars = ws->allVars();
    TIter it(all_ws_vars.createIterator());
    TObject* obj;
    while((obj = it())) {
      TString pname = obj->GetName();
      for(auto& comp : components) {
        TString base = Form("%s_%i_%s_RateSys_", process.Data(), selection, comp.Data());
        if(!pname.BeginsWith(base)) continue;
        TString suffix = pname(base.Length(), pname.Length()-base.Length());
        if(suffix == "") continue;
        TString sys_name = suffix(0, suffix.Length());
        const double value = ((RooRealVar*) obj)->getVal();
        cout << "Found systematic " << sys_name << " for component " << comp
             << " with value " << value
             << endl;
        if(value <= -1.) cout << ">>> WARNING! 1 sigma shift < -1 for " << sys_name << " component "
                              << comp << " (sigma = " << value << ")\n";
        if(value != 0.) { // 0 = no effect
          if(!rate_sys_map.count(sys_name)) { // add this to the map if not already there
            RateSysInfo_t s; s.name_ = sys_name; rate_sys_map[sys_name] = s;
          }
          const double scale = (value <= -1.) ? 1. : (1. + value);
          rate_sys_map[sys_name].scales_[comp] = scale;
        }
      }
    }
  }
  std::vector<RateSysInfo_t> rate_sys_infos;
  for(auto& kv : rate_sys_map) rate_sys_infos.push_back(kv.second);


  //-----------------------------------------------------------------
  // Create nuisance parameters (one theta per systematic)
  //-----------------------------------------------------------------
  std::map<TString, RooRealVar*> theta_map;

  // Shape uncertainties
  for(auto& sinfo : sys_infos) {
    TString tname = Form("theta_%s", sinfo.name_.Data());
    RooRealVar* theta = new RooRealVar(tname, Form("Nuisance param for %s", sinfo.name_.Data()), 0., -5., 5.);
    theta_map[sinfo.name_] = theta;
    nuisanceParams.add(*theta);
    owned_reals.push_back(theta);
    theta->setConstant(!include_systematics_); // freeze if not using systematics

    // Gaussian constraint
    RooRealVar* theta_mean  = new RooRealVar(Form("%s_mean" , tname.Data()), "mean" , 0.); theta_mean->setConstant(true);
    RooRealVar* theta_sigma = new RooRealVar(Form("%s_sigma", tname.Data()), "sigma", 1., 0.01, 10.); theta_sigma->setConstant(true);
    owned_reals.push_back(theta_mean);
    owned_reals.push_back(theta_sigma);
    auto* gauss = new RooGaussian(Form("constraint_%s", sinfo.name_.Data()),
                                  Form("Constraint for %s", sinfo.name_.Data()),
                                  *theta, *theta_mean, *theta_sigma);
    constraintPdfs.add(*gauss);
    owned_pdfs.push_back(gauss);
  }

  // Rate uncertainties
  for(auto& sinfo : rate_sys_infos) {
    // Use a log-normal constraint
    TString tname = Form("kappa_%s", sinfo.name_.Data());
    RooRealVar* kappa = new RooRealVar(tname, Form("Nuisance param for %s", sinfo.name_.Data()), 1., 0.001, 10.);
    // theta_map[sinfo.name_] = kappa;
    // nuisanceParams.add(*kappa);
    // owned_reals.push_back(kappa);
    kappa->setConstant(!include_systematics_); // freeze if not using systematics

    // log-normal constraint
    auto median  = new RooConstVar(Form("%s_median", tname.Data()), "median" , 1.);
    for(auto comp : components) {
      if(sinfo.scales_.count(comp)) {
        auto shape   = new RooConstVar(Form("%s_%s_shape", tname.Data(), comp.Data()), "shape", TMath::Exp(sinfo.scales_[comp]));
        owned_reals.push_back(shape);
        auto* logn = new RooLognormal(Form("constraint_%s_%s", sinfo.name_.Data(), comp.Data()),
                                      Form("Constraint for %s %s", sinfo.name_.Data(), comp.Data()),
                                      *kappa, *median, *shape);
        // constraintPdfs.add(*logn);
        // owned_pdfs.push_back(logn);
      }
    }
    // owned_reals.push_back(theta_mean);
  }

  //-----------------------------------------------------------------
  // Observable and rates
  //-----------------------------------------------------------------
  RooRealVar* obs = get_var(ws, Form("obs_%i", selection));
  if(!obs) return nullptr;

  //-----------------------------------------------------------------
  // Build morphed PDFs and collect extended terms
  //-----------------------------------------------------------------
  RooArgList term_pdfs;
  RooArgList term_yields;

  for(auto& comp : components) {
    TString nom_name = Form("%s_%i_%s_pdf", process.Data(), selection, comp.Data());
    TString nrm_name = Form("%s_%i_%s_norm", process.Data(), selection, comp.Data());
    RooAbsPdf*  nom_pdf = get_pdf(ws, nom_name);
    RooRealVar* norm    = get_var(ws, nrm_name);
    if(!nom_pdf || !norm) {
      std::cout << "stat_analysis: skipping component " << comp << " (missing PDF or norm)\n";
      continue;
    }

    // For the signal component, yield = mu (free parameter)
    // For backgrounds, yield = norm value (fixed to MC prediction)
    RooAbsReal* yield = nullptr;
    if(comp == "signal") {
      yield = mu; // free parameter
    } else {
      // Background yield is fixed
      RooRealVar* bkg_yield = new RooRealVar(Form("yield_%s_%i_%s", process.Data(), selection, comp.Data()),
                                              Form("Yield for %s", comp.Data()),
                                              norm->getVal(), 0., norm->getVal() * 10.);
      bkg_yield->setConstant(true);
      owned_reals.push_back(bkg_yield);
      yield = bkg_yield;
    }

    // Check if this component has any systematics
    bool has_sys = false;
    for(auto& sinfo : sys_infos) {
      if(sinfo.pdf_up_.count(comp) || sinfo.pdf_down_.count(comp)) { has_sys = true; break; }
    }

    RooAbsPdf* effective_pdf = nom_pdf;

    if(has_sys && !sys_infos.empty()) {
      // Build a morphed PDF using histogram interpolation
      // We use a RooHistPdf for the nominal, and add shape variations via
      // a weighted sum approach:
      //   pdf_eff = sum_i w_i * pdf_i
      // where the weights are piecewise-linear functions of the theta parameters.
      //
      // For simplicity we implement it component-by-component by forming the
      // PDF as a RooAddPdf where the fractions are driven by the theta nuisances.
      // When multiple systematics exist we apply them multiplicatively.

      // We'll use a PiecewiseInterpolation-style sum.
      // For each systematic with up/down shapes:
      //   N_eff(theta) = N_nom * (1 + alpha_up * max(0,theta) - alpha_down * max(0,-theta))
      // where alpha_up  = (N_up  - N_nom) / N_nom
      //       alpha_down = (N_nom - N_down) / N_nom
      //
      // The PDF shape morphing is done by reweighting RooHistPdf bins.
      // For RooFit, the cleanest approach without external libraries is to use
      // a pair of RooAddPdf with theta-driven fractions.

      // Collect all (sys_name -> theta) pairs relevant to this component
      std::vector<std::pair<TString,RooRealVar*>> comp_thetas;
      for(auto& sinfo : sys_infos) {
        if(sinfo.pdf_up_.count(comp) && sinfo.pdf_down_.count(comp))
          comp_thetas.push_back({sinfo.name_, theta_map[sinfo.name_]});
      }

      if(comp_thetas.size() >= 1) {
        // Blend nominal with up/down PDFs using the first systematic.
        // pdf_eff = (1 - |theta|/5) * nom
        //         + max(0, theta/5)  * pdf_up
        //         + max(0,-theta/5)  * pdf_down
        // Implemented as a 3-component RooAddPdf where the fractions are
        // RooFormulaVars driven by the nuisance parameter theta.
        // For multiple systematics, chain the morphing sequentially.
        RooAbsPdf* morphed_pdf = nom_pdf;
        for(auto& ct : comp_thetas) {
          const TString& sys_name_ct = ct.first;
          RooRealVar*    theta_ct    = ct.second;
          // Find the SysInfo for this systematic
          const ShapeSysInfo_t* sinfo_ptr = nullptr;
          for(auto& si : sys_infos) {
            if(si.name_ == sys_name_ct) { sinfo_ptr = &si; break; }
          }
          if(!sinfo_ptr) continue;
          if(!sinfo_ptr->pdf_up_.count(comp) || !sinfo_ptr->pdf_down_.count(comp)) continue;
          RooAbsPdf* pdf_up_ct   = get_pdf(ws, sinfo_ptr->pdf_up_  .at(comp));
          RooAbsPdf* pdf_down_ct = get_pdf(ws, sinfo_ptr->pdf_down_.at(comp));
          if(!pdf_up_ct || !pdf_down_ct) continue;

          TString fn_base = Form("morph_%s_%i_%s_%s", process.Data(), selection, comp.Data(), sys_name_ct.Data());
          // f_up   = max(0,  theta/5)  (fraction of up shape)
          // f_down = max(0, -theta/5)  (fraction of down shape)
          // f_nom  = 1 - f_up - f_down (remainder; handled automatically by RooAddPdf)
          RooFormulaVar* fup_ct   = new RooFormulaVar(fn_base+"_fup",
                                                       "max(0., @0/5.)", RooArgList(*theta_ct));
          RooFormulaVar* fdown_ct = new RooFormulaVar(fn_base+"_fdown",
                                                       "max(0.,-@0/5.)", RooArgList(*theta_ct));
          // 3-component mixture: last PDF (nom) gets the remaining fraction automatically
          auto* step_morphed = new RooAddPdf(fn_base+"_pdf",
                                              Form("Morphed PDF for %s sys %s", comp.Data(), sys_name_ct.Data()),
                                              RooArgList(*pdf_up_ct, *pdf_down_ct, *morphed_pdf),
                                              RooArgList(*fup_ct, *fdown_ct));
          owned_reals.push_back(fup_ct);
          owned_reals.push_back(fdown_ct);
          owned_pdfs.push_back(step_morphed);
          morphed_pdf = step_morphed;
        }
        effective_pdf = morphed_pdf;
      }
      // For missing shapes, effective_pdf stays as nom_pdf
    }

    term_pdfs.add(*effective_pdf);
    term_yields.add(*yield);
  }

  if(term_pdfs.getSize() == 0) {
    std::cout << "stat_analysis: ERROR – no component PDFs built!\n";
    return nullptr;
  }

  //-----------------------------------------------------------------
  // Total PDF = sum of component PDFs (extended)
  //-----------------------------------------------------------------
  auto* total_pdf = new RooAddPdf("total_pdf", "Total signal+background PDF",
                                   term_pdfs, term_yields);
  owned_pdfs.push_back(total_pdf);

  //-----------------------------------------------------------------
  // Multiply by constraint terms to get full model.
  // If there are no systematics, the full model is just the total PDF.
  //-----------------------------------------------------------------
  if(constraintPdfs.getSize() == 0) {
    // No constraints: wrap in a trivial RooProdPdf so the return type is consistent
    auto* full_pdf_nosys = new RooProdPdf("full_pdf", "Full model (no systematics)", RooArgSet(*total_pdf));
    owned_pdfs.push_back(full_pdf_nosys);
    return full_pdf_nosys;
  }

  RooArgSet all_terms;
  all_terms.add(*total_pdf);
  all_terms.add(constraintPdfs);

  auto* full_pdf = new RooProdPdf("full_pdf", "Full model with constraints",
                                   all_terms);
  owned_pdfs.push_back(full_pdf);
  return full_pdf;
}

//---------------------------------------------------------------------------
// Compute total expected background yield
//---------------------------------------------------------------------------
double total_background(RooWorkspace* ws, const TString process, const int selection) {
  const bool is_mumem = (process == "mumem");
  std::vector<TString> bkg_comps = {"cosmic","rpc_ext","rpc_int"};
  if(!use_evtana_) bkg_comps.push_back("pbar");
  if(is_mumem) bkg_comps.insert(bkg_comps.begin(), "dio");
  else { bkg_comps.push_back("rmc_ext"); bkg_comps.push_back("rmc_int"); }
  double total = 0.;
  for(auto& comp : bkg_comps) {
    RooRealVar* n = get_var(ws, Form("%s_%i_%s_norm", process.Data(), selection, comp.Data()));
    if(n) total += n->getVal();
  }
  return total;
}

//---------------------------------------------------------------------------
// Print a divider line
//---------------------------------------------------------------------------
void print_divider(const char* title = "") {
  printf("\n%s\n", std::string(70,'=').c_str());
  if(title && strlen(title)) printf("  %s\n%s\n", title, std::string(70,'-').c_str());
}

//---------------------------------------------------------------------------
// Perform profile-likelihood fit for the signal rate
//---------------------------------------------------------------------------
int do_likelihood_fit(
    RooWorkspace* ws,
    const TString process,
    const int selection,
    const TString tag,
    RooProdPdf* full_pdf,
    RooAbsData* data,
    RooRealVar* mu,
    RooArgSet& nuisanceParams,
    const TString figdir
) {
  print_divider("Profile Likelihood Fit");

  // Keep mu >= 0 (physical signal rate); use a small negative floor only
  // to allow the fitter to find mu=0 without hitting the boundary.
  mu->setConstant(false);
  mu->setRange(0., 20.);
  mu->setVal(0.);

  RooFitResult* res = full_pdf->fitTo(*data,
                                       RooFit::Extended(true),
                                       RooFit::Save(true),
                                       RooFit::PrintLevel(-1),
                                       RooFit::Minimizer("Minuit2","Migrad"),
                                       RooFit::Strategy(2));
  if(!res) {
    std::cout << "stat_analysis: likelihood fit failed!\n";
    return 1;
  }

  printf("  Fit status    : %i (covQual=%i)\n", res->status(), res->covQual());
  printf("  Signal yield  : %.4f +/- %.4f events\n", mu->getVal(), mu->getError());

  // Retrieve physics reference info
  RooRealVar* ref_br   = get_var(ws, "ref_signal_br");
  RooRealVar* sig_norm = get_var(ws, Form("%s_%i_signal_norm", process.Data(), selection));

  if(ref_br && sig_norm && sig_norm->getVal() > 0.) {
    const double br_scale = ref_br->getVal() / sig_norm->getVal();
    const double br_best   = mu->getVal()  * br_scale;
    const double br_err    = mu->getError()* br_scale;
    printf("  Signal BR     : %.3e +/- %.3e\n", br_best, br_err);
  }

  printf("\n  Background yields:\n");
  const double bkg_total = total_background(ws, process, selection);
  printf("  Total background: %.4f events\n", bkg_total);

  //-----------------------------------------------------------------
  // Plot: data vs best-fit model
  //-----------------------------------------------------------------
  RooRealVar* obs = get_var(ws, Form("obs_%i", selection));
  if(obs) {
    TCanvas* c = new TCanvas("c_fit", "Profile likelihood fit", 800, 600);
    auto frame = obs->frame();
    frame->SetTitle("Best-fit model");
    frame->SetXTitle("Momentum (MeV/c)");
    frame->SetYTitle(Form("Events / %.1f MeV/c", (obs->getMax()-obs->getMin())/obs->getBins()));
    data->plotOn(frame, RooFit::Name("data"));
    full_pdf->plotOn(frame, RooFit::Name("total"), RooFit::LineColor(kBlue));
    // Signal-only component
    RooAbsPdf* sig_pdf = get_pdf(ws, Form("%s_%i_signal_pdf", process.Data(), selection));
    if(sig_pdf)
      sig_pdf->plotOn(frame, RooFit::Name("signal"),
                      RooFit::LineColor(kRed), RooFit::LineStyle(kDashed),
                      RooFit::Normalization(std::max(0., mu->getVal()), RooAbsReal::NumEvent));
    frame->Draw();
    TLegend* leg = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg->AddEntry("data"  , "Data"     , "PE");
    leg->AddEntry("total" , "Total fit", "L" );
    if(sig_pdf) leg->AddEntry("signal", "Signal", "L");
    leg->SetLineWidth(0); leg->Draw();
    c->SaveAs(Form("%s/likelihood_fit_%i.png", figdir.Data(), selection));
    delete c;
    delete frame;
  }

  delete res;
  return 0;
}

//---------------------------------------------------------------------------
// CLs asymptotic upper limit (RooStats AsymptoticCalculator)
//---------------------------------------------------------------------------
int do_cls_asymptotic(
    RooWorkspace* ws,
    const TString process,
    const int selection,
    const TString tag,
    RooProdPdf* full_pdf,
    RooAbsData* data,
    RooRealVar* mu,
    RooArgSet& nuisanceParams,
    const TString figdir
) {
  print_divider("CLs Asymptotic Upper Limit (90% CL)");

  RooRealVar* obs = get_var(ws, Form("obs_%i", selection));
  if(!obs) return 1;

  //-----------------------------------------------------------------
  // Estimate a sensible scan range.
  // The 90% CLs UL for a counting experiment with n_obs observed and
  // b expected background is approximately:
  //   mu_UL ~ (n_obs - b + 1.28*sqrt(n_obs)) / efficiency
  // We use 5x this as the upper bound, and fall back to 10 if it is
  // too small.  The lower bound is a small positive epsilon so that
  // the q_mu test statistic is well-defined.
  //-----------------------------------------------------------------
  const double n_obs  = data->sumEntries();
  const double n_bkg  = total_background(ws, process, selection);
  RooRealVar* sig_norm_var = get_var(ws, Form("%s_%i_signal_norm", process.Data(), selection));
  const double sig_eff = (sig_norm_var && sig_norm_var->getVal() > 0.) ? sig_norm_var->getVal() : 1.;
  // Conservative upper bound: Poisson 99.9% upper limit on (n_obs-n_bkg) / eff
  const double n_excess  = std::max(n_obs - n_bkg, 0.) + 3.*std::sqrt(std::max(n_obs, 1.));
  const double mu_max_est = 5. * n_excess / sig_eff;
  const double mu_scan_max = std::max(mu_max_est, 10. / sig_eff);
  const double mu_scan_min = 1.e-3;  // small positive: avoid boundary at 0
  printf("  Scan range: [%.4f, %.4f] events  (n_obs=%.0f, n_bkg=%.3f, sig_eff=%.4f)\n",
         mu_scan_min, mu_scan_max, n_obs, n_bkg, sig_eff);

  mu->setRange(0., mu_scan_max * 1.2);
  mu->setVal(1.);
  mu->setConstant(false);

  ModelConfig sb_model("S+B model", ws);
  sb_model.SetPdf(*full_pdf);
  sb_model.SetObservables(*obs);
  sb_model.SetParametersOfInterest(*mu);
  sb_model.SetNuisanceParameters(nuisanceParams);
  { RooArgSet s; s.add(*mu); sb_model.SetSnapshot(s); }

  // Background-only snapshot: mu=0 using the actual POI object
  mu->setVal(0.);
  ModelConfig bkg_model("B-only model", ws);
  bkg_model.SetPdf(*full_pdf);
  bkg_model.SetObservables(*obs);
  bkg_model.SetParametersOfInterest(*mu);
  bkg_model.SetNuisanceParameters(nuisanceParams);
  { RooArgSet s; s.add(*mu); bkg_model.SetSnapshot(s); }

  mu->setVal(0.); mu->setConstant(false);

  AsymptoticCalculator ac(*data, sb_model, bkg_model);
  ac.SetOneSidedDiscovery(false); // upper limit (not discovery)
  ac.SetQTilde(true);             // q~_mu: physical mu>=0

  HypoTestInverter inverter(ac);
  inverter.SetConfidenceLevel(0.90);
  inverter.UseCLs(true);
  inverter.SetVerbose(false);
  // 30 equidistant points from a small epsilon to mu_scan_max
  inverter.SetFixedScan(30, mu_scan_min, mu_scan_max);

  HypoTestInverterResult* result = inverter.GetInterval();
  if(!result) {
    std::cout << "stat_analysis: CLs asymptotic calculation failed!\n";
    return 1;
  }

  const double ul_obs  = result->UpperLimit();
  const double ul_exp  = result->GetExpectedUpperLimit( 0);
  const double ul_ep1  = result->GetExpectedUpperLimit(+1);
  const double ul_em1  = result->GetExpectedUpperLimit(-1);
  const double ul_ep2  = result->GetExpectedUpperLimit(+2);
  const double ul_em2  = result->GetExpectedUpperLimit(-2);

  printf("  Observed 90%% CL upper limit : %.4f events\n", ul_obs);
  printf("  Expected 90%% CL upper limit : %.4f events\n", ul_exp);
  printf("  Expected +1sigma band        : [%.4f, %.4f]\n", ul_em1, ul_ep1);
  printf("  Expected +2sigma band        : [%.4f, %.4f]\n", ul_em2, ul_ep2);

  // Convert to BR if possible
  RooRealVar* ref_br   = get_var(ws, "ref_signal_br");
  RooRealVar* ref_snorm = get_var(ws, Form("%s_%i_signal_norm", process.Data(), selection));
  double br_scale = 0.;
  if(ref_br && ref_snorm && ref_snorm->getVal() > 0.) {
    br_scale = ref_br->getVal() / ref_snorm->getVal();
    printf("  Observed 90%% CL UL on BR    : %.3e\n", ul_obs * br_scale);
    printf("  Expected 90%% CL UL on BR    : %.3e\n", ul_exp * br_scale);
  }

  // Linear y-axis plot via HypoTestInverterPlot
  {
    TCanvas* c = new TCanvas("c_cls_asym_lin", "CLs asymptotic (linear)", 800, 600);
    HypoTestInverterPlot hplot("cls_asym_plot", "", result);
    hplot.Draw("CLb 2CL");
    c->SaveAs(Form("%s/cls_asymptotic_%i.png", figdir.Data(), selection));
    delete c;
  }
  // Log y-axis: redraw on a fresh canvas and set log after drawing
  {
    TCanvas* c = new TCanvas("c_cls_asym_log", "CLs asymptotic (log)", 800, 600);
    HypoTestInverterPlot hplot_log("cls_asym_plot_log", "", result);
    hplot_log.Draw("CLb 2CL");
    c->SetLogy();
    c->SaveAs(Form("%s/cls_asymptotic_%i_log.png", figdir.Data(), selection));
    delete c;
  }

  delete result;
  return 0;
}

//---------------------------------------------------------------------------
// DEBUG STUB – manual CLs band reconstruction (not yet working correctly).
// The asymptotic HypoTestInverterResult does not expose a per-point expected
// CLs API (GetExpectedCLs does not exist in this ROOT version).  This
// function is kept here for future development and is not called by default.
//---------------------------------------------------------------------------
void debug_cls_bands_plot(HypoTestInverterResult* result,
                          double ul_obs, double ul_exp, double br_scale,
                          const TString& figdir, int selection)
{
  // TODO: find the correct per-point API for expected CLs bands from
  // HypoTestInverterResult in this ROOT version, or rerun the asymptotic
  // calculator on the Asimov dataset to obtain the expected CLs curve.
  // Options to investigate:
  //   - result->GetExpectedPValueDist(i)  (returns null for asymptotic)
  //   - Re-running AsymptoticCalculator with SetUseAsymptoticFormula and an
  //     Asimov dataset generated at mu=0 as the "observed" data
  //   - Manually computing 1-Phi(mu/sigma_A) from sigma_A = ul_exp/1.28155
  //     per scan point (approximate, ignores CLb denominator shape)
  std::cout << "debug_cls_bands_plot: not yet implemented\n";
  (void)result; (void)ul_obs; (void)ul_exp; (void)br_scale;
  (void)figdir; (void)selection;
}

//---------------------------------------------------------------------------
// Helper: fit full_pdf to a dataset with mu constrained to [mu_lo, mu_hi].
// Returns -2*log(L) at the best-fit point, or 1e30 on failure.
//---------------------------------------------------------------------------
double fit_nll(RooProdPdf* full_pdf, RooAbsData* d,
               RooRealVar* mu, RooArgSet& nuisanceParams,
               double mu_lo, double mu_hi)
{
  mu->setRange(mu_lo, mu_hi);
  mu->setVal(0.5*(mu_lo + mu_hi));
  mu->setConstant(false);
  reset_nuisances(nuisanceParams);
  RooFitResult* r = full_pdf->fitTo(*d,
                                     RooFit::Extended(true),
                                     RooFit::Save(true),
                                     RooFit::PrintLevel(-1),
                                     RooFit::Minimizer("Minuit2","Migrad"),
                                     RooFit::Strategy(1));
  const double nll = (r && r->status() == 0) ? r->minNll() : 1.e30;
  delete r;
  return nll;
}

//---------------------------------------------------------------------------
// CLs toy-based upper limit – manual loop to avoid ToyMCSampler's
// internal setData() incompatibility with cached PDFs (RooNumConvPdf etc.)
//
// For each scan point mu_test we:
//   1. Generate ntoys datasets from the B-only PDF (mu=0)  -> null toys
//   2. Generate ntoys datasets from the S+B PDF (mu=mu_test) -> alt toys
//   3. For each toy compute q_mu = -2*ln(L(mu_test)/L(mu_hat))
//      with the one-sided constraint mu_hat >= 0
//   4. CLs(mu_test) = ps+b / pb  where ps+b (pb) is the fraction of
//      S+B (B-only) toys with q_mu >= q_mu_obs
//   5. 90% UL is the mu_test where CLs = 0.10
//---------------------------------------------------------------------------
int do_cls_toys(
    RooWorkspace* ws,
    const TString process,
    const int selection,
    const TString tag,
    RooProdPdf* full_pdf,
    RooAbsData* data,
    RooRealVar* mu,
    RooArgSet& nuisanceParams,
    const TString figdir,
    const int ntoys = 1000
) {
  print_divider("CLs Toy-Based Upper Limit (90% CL)");
  printf("  Using %i toys per scan point\n", ntoys);

  RooRealVar* obs = get_var(ws, Form("obs_%i", selection));
  if(!obs) return 1;

  // Scan range
  const double n_obs_t  = data->sumEntries();
  const double n_bkg_t  = total_background(ws, process, selection);
  RooRealVar* sig_norm_t = get_var(ws, Form("%s_%i_signal_norm", process.Data(), selection));
  const double sig_eff_t = (sig_norm_t && sig_norm_t->getVal() > 0.) ? sig_norm_t->getVal() : 1.;
  const double n_excess_t    = std::max(n_obs_t - n_bkg_t, 0.) + 3.*std::sqrt(std::max(n_obs_t, 1.));
  const double mu_scan_max_t = std::max(5.*n_excess_t / sig_eff_t, 10./sig_eff_t);
  const int    n_scan        = 10;
  printf("  Scan range: [0, %.4f] events\n", mu_scan_max_t);

  mu->setRange(0., mu_scan_max_t * 1.2);

  // Compute q_mu_obs on the real data for each scan point
  // q_mu = max(0, -2*(NLL(mu_test) - NLL(mu_hat_free)))  [one-sided]
  mu->setRange(0., mu_scan_max_t * 1.2);
  const double nll_free_obs = fit_nll(full_pdf, data, mu, nuisanceParams,
                                       0., mu_scan_max_t * 1.2);

  // Vectors to accumulate CLs scan
  std::vector<double> mu_scan_pts, cls_obs_pts;

  // For plots at the last scan point (most diagnostic)
  std::vector<double> q_bonly_last, q_spb_last, muhat_bonly_last, muhat_spb_last;
  double q_obs_last = 0., mu_test_last = 0.;

  for(int iscan = 0; iscan < n_scan; ++iscan) {
    const double mu_test = mu_scan_max_t * (iscan + 1.) / n_scan;

    // Observed test statistic q_mu on real data
    reset_nuisances(nuisanceParams);
    const double nll_fixed_obs = fit_nll(full_pdf, data, mu, nuisanceParams,
                                          mu_test, mu_test);
    reset_nuisances(nuisanceParams);
    const double nll_free_obs  = fit_nll(full_pdf, data, mu, nuisanceParams,
                                          0., mu_scan_max_t * 1.2);
    const double q_obs = std::max(0., 2.*(nll_fixed_obs - nll_free_obs));

    // B-only toys
    mu->setVal(0.); mu->setConstant(true);
    reset_nuisances(nuisanceParams);
    int n_above_bonly = 0;
    std::vector<double> q_b_vec, muhat_b_vec;
    for(int itoy = 0; itoy < ntoys; ++itoy) {
      RooDataHist* toy = (RooDataHist*) full_pdf->generateBinned(*obs, RooFit::Extended(true));
      if(!toy) continue;
      mu->setConstant(false);
      reset_nuisances(nuisanceParams);
      const double nll_free  = fit_nll(full_pdf, toy, mu, nuisanceParams, 0., mu_scan_max_t*1.2);
      const double muhat     = mu->getVal();
      reset_nuisances(nuisanceParams);
      const double nll_fixed = fit_nll(full_pdf, toy, mu, nuisanceParams, mu_test, mu_test);
      const double q_toy     = std::max(0., 2.*(nll_fixed - nll_free));
      if(q_toy >= q_obs) ++n_above_bonly;
      q_b_vec.push_back(q_toy);
      muhat_b_vec.push_back(muhat);
      delete toy;
    }
    const double p_b = (double)n_above_bonly / ntoys;

    // S+B toys at mu_test
    mu->setVal(mu_test); mu->setConstant(true);
    reset_nuisances(nuisanceParams);
    int n_above_spb = 0;
    std::vector<double> q_spb_vec, muhat_spb_vec;
    for(int itoy = 0; itoy < ntoys; ++itoy) {
      RooDataHist* toy = (RooDataHist*) full_pdf->generateBinned(*obs, RooFit::Extended(true));
      if(!toy) continue;
      mu->setConstant(false);
      reset_nuisances(nuisanceParams);
      const double nll_free  = fit_nll(full_pdf, toy, mu, nuisanceParams, 0., mu_scan_max_t*1.2);
      const double muhat     = mu->getVal();
      reset_nuisances(nuisanceParams);
      const double nll_fixed = fit_nll(full_pdf, toy, mu, nuisanceParams, mu_test, mu_test);
      const double q_toy     = std::max(0., 2.*(nll_fixed - nll_free));
      if(q_toy >= q_obs) ++n_above_spb;
      q_spb_vec.push_back(q_toy);
      muhat_spb_vec.push_back(muhat);
      delete toy;
    }
    const double p_spb = (double)n_above_spb / ntoys;

    mu->setConstant(false); reset_nuisances(nuisanceParams);

    const double cls = (p_b > 0.) ? std::min(p_spb / p_b, 1.) : 1.;
    mu_scan_pts.push_back(mu_test);
    cls_obs_pts.push_back(cls);

    printf("  mu_test=%.4f  q_obs=%.4f  p_b=%.3f  p_spb=%.3f  CLs=%.3f\n",
           mu_test, q_obs, p_b, p_spb, cls);

    // Save distributions at last point for plots
    q_bonly_last   = q_b_vec;    muhat_bonly_last = muhat_b_vec;
    q_spb_last     = q_spb_vec;  muhat_spb_last   = muhat_spb_vec;
    q_obs_last     = q_obs;      mu_test_last      = mu_test;
  }

  // Find 90% UL by linear interpolation
  double ul_obs_toys = mu_scan_max_t;
  for(int i = 1; i < (int)mu_scan_pts.size(); ++i) {
    if(cls_obs_pts[i-1] >= 0.10 && cls_obs_pts[i] < 0.10) {
      double frac = (cls_obs_pts[i-1] - 0.10) / (cls_obs_pts[i-1] - cls_obs_pts[i]);
      ul_obs_toys = mu_scan_pts[i-1] + frac*(mu_scan_pts[i] - mu_scan_pts[i-1]);
      break;
    }
  }

  printf("  Observed 90%% CL upper limit (toys): %.4f events\n", ul_obs_toys);
  RooRealVar* ref_br   = get_var(ws, "ref_signal_br");
  RooRealVar* ref_sn   = get_var(ws, Form("%s_%i_signal_norm", process.Data(), selection));
  if(ref_br && ref_sn && ref_sn->getVal() > 0.)
    printf("  Observed 90%% CL UL on BR (toys): %.3e\n",
           ul_obs_toys * ref_br->getVal() / ref_sn->getVal());

  //----- Plot 1: CLs vs mu (linear + log) ----
  {
    int ng = (int)mu_scan_pts.size();
    TGraph g_cls(ng, mu_scan_pts.data(), cls_obs_pts.data());
    g_cls.SetLineColor(kBlack); g_cls.SetLineWidth(2); g_cls.SetMarkerStyle(20);
    TLine cl_line(mu_scan_pts.front(), 0.10, mu_scan_pts.back(), 0.10);
    cl_line.SetLineColor(kRed); cl_line.SetLineStyle(kDashed); cl_line.SetLineWidth(2);
    auto save_cls = [&](bool logy, const char* fname) {
      TCanvas* c = new TCanvas("ctmp","",800,600);
      double ymin = logy ? 1.e-3 : 0.;
      TH1* fr = c->DrawFrame(mu_scan_pts.front(), ymin, mu_scan_pts.back()*1.02, 1.2);
      fr->SetXTitle("Signal yield #mu (events)"); fr->SetYTitle("CL_{s}"); fr->SetTitle("");
      g_cls.Draw("LP same"); cl_line.Draw("same");
      if(logy) c->SetLogy();
      TLatex lat; lat.SetNDC(); lat.SetTextSize(0.032); lat.SetTextFont(42);
      lat.DrawLatex(0.14, 0.85, Form("Obs. 90%% UL: #mu < %.3g", ul_obs_toys));
      c->SaveAs(fname); delete c;
    };
    save_cls(false, Form("%s/cls_toys_%i.png",     figdir.Data(), selection));
    save_cls(true,  Form("%s/cls_toys_%i_log.png", figdir.Data(), selection));
  }

  //----- Plot 2: best-fit mu_hat distribution (B-only and S+B toys) ----
  if(!muhat_bonly_last.empty() && !muhat_spb_last.empty()) {
    double muhat_lo = -mu_test_last*0.5, muhat_hi = mu_test_last*2.;
    TH1F h_mub ("h_mub" , Form(";#hat{#mu} (events);Toys / %.2f", (muhat_hi-muhat_lo)/50.),
                50, muhat_lo, muhat_hi);
    TH1F h_mus ("h_mus" , ";#hat{#mu} (events);Toys", 50, muhat_lo, muhat_hi);
    for(auto v : muhat_bonly_last) h_mub.Fill(v);
    for(auto v : muhat_spb_last)   h_mus.Fill(v);
    h_mub.SetLineColor(kBlue);  h_mub.SetFillColorAlpha(kBlue, 0.3);
    h_mus.SetLineColor(kRed);   h_mus.SetFillColorAlpha(kRed,  0.3);
    TCanvas* c = new TCanvas("c_muhat_toys","",800,600);
    h_mub.Draw("hist"); h_mus.Draw("hist same");
    TLine lmu(mu_test_last, 0., mu_test_last, std::max(h_mub.GetMaximum(),h_mus.GetMaximum()));
    lmu.SetLineColor(kBlack); lmu.SetLineStyle(kDashed); lmu.SetLineWidth(2); lmu.Draw();
    TLegend leg(0.55,0.70,0.88,0.88); leg.SetFillStyle(0); leg.SetBorderSize(0);
    leg.AddEntry(&h_mub, "B-only toys",          "F");
    leg.AddEntry(&h_mus, Form("S+B toys (#mu=%.2f)",mu_test_last), "F");
    leg.AddEntry(&lmu,   "#mu_{test}",            "L");
    leg.Draw();
    c->SaveAs(Form("%s/cls_toys_muhat_%i.png", figdir.Data(), selection));
    delete c;
  }

  //----- Plot 3: signal rate pull (B-only toys: (mu_hat - 0)/sigma_mu) ----
  // We don't have sigma_mu from fit_nll, so use RMS as denominator
  if(!muhat_bonly_last.empty()) {
    double sum=0.; for(auto v:muhat_bonly_last) sum+=v;
    double mean_b = sum/muhat_bonly_last.size();
    double rms2=0.; for(auto v:muhat_bonly_last) rms2+=(v-mean_b)*(v-mean_b);
    double rms_b = (muhat_bonly_last.size()>1) ? std::sqrt(rms2/(muhat_bonly_last.size()-1)) : 1.;
    TH1F h_pull("h_pull_toys",";(#hat{#mu} - 0) / RMS;Toys",50,-5.,5.);
    for(auto v:muhat_bonly_last) h_pull.Fill((rms_b>0.) ? v/rms_b : 0.);
    TCanvas* c = new TCanvas("c_pull_toys","",800,600);
    h_pull.Fit("gaus","Q"); h_pull.Draw();
    TLine lz(0.,0.,0.,h_pull.GetMaximum());
    lz.SetLineColor(kRed); lz.SetLineWidth(2); lz.Draw();
    c->SaveAs(Form("%s/cls_toys_pull_%i.png", figdir.Data(), selection));
    delete c;
  }

  //----- Plot 4: q_mu test-statistic distributions at last scan point ----
  if(!q_bonly_last.empty() && !q_spb_last.empty()) {
    double qmax = *std::max_element(q_bonly_last.begin(), q_bonly_last.end());
    qmax = std::max(qmax, *std::max_element(q_spb_last.begin(), q_spb_last.end()));
    qmax = std::max(qmax, q_obs_last * 1.5);
    TH1F h_qb("h_qb", Form(";q_{#mu} (#mu_{test}=%.2f);Toys",mu_test_last), 50,0.,qmax);
    TH1F h_qs("h_qs", "", 50, 0., qmax);
    for(auto v:q_bonly_last) h_qb.Fill(v);
    for(auto v:q_spb_last)   h_qs.Fill(v);
    h_qb.SetLineColor(kBlue); h_qb.SetFillColorAlpha(kBlue, 0.3);
    h_qs.SetLineColor(kRed);  h_qs.SetFillColorAlpha(kRed,  0.3);
    TCanvas* c = new TCanvas("c_qmu_toys","",800,600);
    c->SetLogy();
    h_qb.Draw("hist"); h_qs.Draw("hist same");
    TLine lq(q_obs_last, 0., q_obs_last, std::max(h_qb.GetMaximum(),h_qs.GetMaximum()));
    lq.SetLineColor(kBlack); lq.SetLineWidth(2); lq.Draw();
    TLegend leg(0.50,0.70,0.88,0.88); leg.SetFillStyle(0); leg.SetBorderSize(0);
    leg.AddEntry(&h_qb, "B-only toys",          "F");
    leg.AddEntry(&h_qs, Form("S+B toys (#mu=%.2f)",mu_test_last), "F");
    leg.AddEntry(&lq,   "q_{#mu}^{obs}",         "L");
    leg.Draw();
    c->SaveAs(Form("%s/cls_toys_qmu_%i.png", figdir.Data(), selection));
    delete c;
  }

  mu->setVal(0.); mu->setConstant(false);
  return 0;
}
//---------------------------------------------------------------------------
// Toy likelihood fits
//---------------------------------------------------------------------------
int do_fit_toys(
    RooWorkspace* ws,
    const TString process,
    const int selection,
    const TString tag,
    RooProdPdf* full_pdf,
    RooAbsData* data,
    RooRealVar* mu,
    RooArgSet& nuisanceParams,
    const TString figdir,
    const int ntoys = 1000,
    const double mu_inject = 0.
) {
  print_divider("Toy likelihood fits");
  printf("  Using %i toys\n", ntoys);

  RooRealVar* obs = get_var(ws, Form("obs_%i", selection));
  if(!obs) return 1;

  // Scan range
  const double n_bkg_t  = total_background(ws, process, selection);
  RooRealVar* sig_norm_t = get_var(ws, Form("%s_%i_signal_norm", process.Data(), selection));
  const double sig_eff_t = (sig_norm_t && sig_norm_t->getVal() > 0.) ? sig_norm_t->getVal() : 1.;

  vector<pair<double,double>> fit_results;
  double min_val(1.e10), max_val(-1.e10);
  for(int itoy = 0; itoy < ntoys; ++itoy) {
    reset_nuisances(nuisanceParams);
    mu->setVal(mu_inject); mu->setConstant(true);
    RooDataHist* toy = (RooDataHist*) full_pdf->generateBinned(*obs, RooFit::Extended(true));
    if(!toy) continue;
    mu->setConstant(false);
    reset_nuisances(nuisanceParams);
    const double nll_free  = fit_nll(full_pdf, toy, mu, nuisanceParams, -1000., 1000.);
    fit_results.push_back(pair<double,double>(mu->getVal(), mu->getError()));
    min_val = min(mu->getVal(), min_val);
    max_val = max(mu->getVal(), max_val);
    delete toy;
  }

  // Create the histograms
  TH1* h_mus = new TH1D("mus", "Toy fit signal yields;#mu;", 40, min_val - 0.1, max_val + 0.1);
  TH1* h_pulls = new TH1D("pulls", "Toy fit signal yield pulls;(#mu_{fit} - #mu_{true})/#sigma_{fit};", 40, -5., 5.);
  for(auto & [val, err] : fit_results) {
    const double pull = (err > 0.) ? (val - mu_inject) / err : 10.;
    h_mus->Fill(val);
    h_pulls->Fill(pull);
  }

  // Plot the results
  TCanvas* c = new TCanvas("c_toys","",1600,600);
  c->Divide(2,1);
  c->cd(1);
  h_mus->SetLineColor(kBlack);
  h_mus->SetFillColor(kAzure-4);
  h_mus->Draw("hist");

  c->cd(2);
  h_pulls->SetLineColor(kBlack);
  h_pulls->SetFillColor(kAzure-4);
  h_pulls->Draw("hist");
  h_pulls->Fit("gaus", "Q");
  h_pulls->GetListOfFunctions()->At(0)->Draw("same");
  c->SaveAs(Form("%s/fit_toys_%i.png", figdir.Data(), selection));
  delete c;
  delete h_mus;
  delete h_pulls;

  return 0;
}

//---------------------------------------------------------------------------
// Feldman-Cousins 90% CI – manual Neyman belt construction.
//
// For each mu_true on a grid, generate ntoys toys and find the central
// 90% interval of the profile-likelihood test statistic q_mu.  The FC
// confidence interval for the observed data is the set of mu_true values
// whose belt contains the observed q_mu(mu_true).
//
// Uses fit_nll() directly to avoid ToyMCSampler::setData() crashes with
// cached PDFs (RooNumConvPdf etc.).
//---------------------------------------------------------------------------
int do_feldman_cousins(
    RooWorkspace* ws,
    const TString process,
    const int selection,
    const TString tag,
    RooProdPdf* full_pdf,
    RooAbsData* data,
    RooRealVar* mu,
    RooArgSet& nuisanceParams,
    const TString figdir,
    const int ntoys = 500
) {
  print_divider("Feldman-Cousins 90% Confidence Interval");
  printf("  Using %i toys per point\n", ntoys);

  RooRealVar* obs = get_var(ws, Form("obs_%i", selection));
  if(!obs) return 1;

  // Scan range
  const double n_obs_fc  = data->sumEntries();
  const double n_bkg_fc  = total_background(ws, process, selection);
  RooRealVar* sig_norm_fc = get_var(ws, Form("%s_%i_signal_norm", process.Data(), selection));
  const double sig_eff_fc = (sig_norm_fc && sig_norm_fc->getVal()>0.) ? sig_norm_fc->getVal() : 1.;
  const double n_excess_fc   = std::max(n_obs_fc - n_bkg_fc, 0.) + 3.*std::sqrt(std::max(n_obs_fc,1.));
  const double mu_max_fc     = std::max(5.*n_excess_fc / sig_eff_fc, 10./sig_eff_fc);
  const int    n_scan_fc     = 10;
  printf("  Scan range: [0, %.4f] events\n", mu_max_fc);

  mu->setRange(0., mu_max_fc * 1.2);

  // Observed free-fit NLL (denominator)
  reset_nuisances(nuisanceParams);
  const double nll_free_obs = fit_nll(full_pdf, data, mu, nuisanceParams, 0., mu_max_fc*1.2);

  // For each mu_true: compute q_mu_obs and the toy belt quantiles
  // FC interval = {mu_true : q_mu_obs <= q_mu_90pct(mu_true)}
  std::vector<double> mu_grid, fc_lo_v, fc_hi_v;
  bool in_interval = false;
  double fc_lo_val = 0., fc_hi_val = 0.;

  for(int iscan = 0; iscan <= n_scan_fc; ++iscan) {
    const double mu_true = mu_max_fc * iscan / n_scan_fc;

    // Observed q_mu at this mu_true
    reset_nuisances(nuisanceParams);
    const double nll_fixed_obs = fit_nll(full_pdf, data, mu, nuisanceParams, mu_true, mu_true);
    const double q_obs_fc = std::max(0., 2.*(nll_fixed_obs - nll_free_obs));

    // Generate toys at mu_true and collect q_mu values
    mu->setVal(mu_true); mu->setConstant(true);
    reset_nuisances(nuisanceParams);
    std::vector<double> q_toys;
    q_toys.reserve(ntoys);
    for(int itoy = 0; itoy < ntoys; ++itoy) {
      RooDataHist* toy = (RooDataHist*) full_pdf->generateBinned(*obs, RooFit::Extended(true));
      if(!toy) continue;
      mu->setConstant(false); reset_nuisances(nuisanceParams);
      const double nll_free_t  = fit_nll(full_pdf, toy, mu, nuisanceParams, 0., mu_max_fc*1.2);
      reset_nuisances(nuisanceParams);
      const double nll_fixed_t = fit_nll(full_pdf, toy, mu, nuisanceParams, mu_true, mu_true);
      q_toys.push_back(std::max(0., 2.*(nll_fixed_t - nll_free_t)));
      delete toy;
    }
    mu->setConstant(false); reset_nuisances(nuisanceParams);

    if(q_toys.empty()) { printf("  mu_true=%.4f  no toys converged\n", mu_true); continue; }

    // 90% quantile of q_mu distribution (upper edge of acceptance region)
    std::sort(q_toys.begin(), q_toys.end());
    const double q90 = q_toys[(int)(0.90 * q_toys.size())];

    // Data point is inside belt if q_obs <= q90
    const bool inside = (q_obs_fc <= q90);
    printf("  mu_true=%.4f  q_obs=%.4f  q90=%.4f  %s\n",
           mu_true, q_obs_fc, q90, inside ? "IN" : "out");

    mu_grid.push_back(mu_true);
    if(inside && !in_interval) { fc_lo_val = mu_true; in_interval = true; }
    if(!inside && in_interval) { fc_hi_val = mu_true; in_interval = false; }
  }
  // If still inside at the end of scan, close the interval
  if(in_interval) fc_hi_val = mu_max_fc;

  printf("  FC 90%% CI: [%.4f, %.4f] events\n", fc_lo_val, fc_hi_val);

  RooRealVar* ref_br   = get_var(ws, "ref_signal_br");
  RooRealVar* ref_sn   = get_var(ws, Form("%s_%i_signal_norm", process.Data(), selection));
  if(ref_br && ref_sn && ref_sn->getVal() > 0.) {
    const double scale = ref_br->getVal() / ref_sn->getVal();
    printf("  FC 90%% CI on BR: [%.3e, %.3e]\n", fc_lo_val*scale, fc_hi_val*scale);
  }

  mu->setVal(0.); mu->setConstant(false);
  return 0;
}

//---------------------------------------------------------------------------
// Toy MC closure / pull study
//---------------------------------------------------------------------------
int do_toy_study(
    RooWorkspace* ws,
    const TString process,
    const int selection,
    const TString tag,
    RooProdPdf* full_pdf,
    RooRealVar* mu,
    RooArgSet& nuisanceParams,
    const TString figdir,
    const int ntoys = 500,
    const double mu_inject = 0.  // injected signal strength
) {
  print_divider("Toy MC Closure / Pull Study");
  printf("  %i toys with mu_inject = %.3f events\n", ntoys, mu_inject);

  RooRealVar* obs = get_var(ws, Form("obs_%i", selection));
  if(!obs) return 1;

  // Total expected events for Poisson fluctuation
  const double bkg_total = total_background(ws, process, selection);
  const double n_exp = bkg_total + mu_inject;
  printf("  Expected events: %.2f (signal=%.3f + bkg=%.3f)\n", n_exp, mu_inject, bkg_total);

  // Store fit results
  std::vector<double> mu_fit_vals, mu_fit_errs, pulls;
  std::vector<double> ul_vals; // toy upper limits (asymptotic)

  RooRealVar* sig_norm_ts = get_var(ws, Form("%s_%i_signal_norm", process.Data(), selection));
  const double sig_eff_ts  = (sig_norm_ts && sig_norm_ts->getVal() > 0.) ? sig_norm_ts->getVal() : 1.;
  const double mu_range_ts = std::max(mu_inject * 3., 10. / sig_eff_ts);
  TH1F h_mu  ("h_mu_toys"  , Form(";#hat{#mu} (events);Toys / %.2f", 2.*mu_range_ts/50.),
              50, -mu_range_ts, mu_range_ts);
  TH1F h_pull("h_pull_toys", ";(#hat{#mu} - #mu_{inject}) / #sigma_{#mu};Toys", 50, -5., 5.);
  TH1F h_ul  ("h_ul_toys"  , ";Approx. 90% UL (events);Toys", 50, 0., mu_range_ts*2.);

  // Set mu to injected value for generation
  mu->setRange(0., mu_range_ts * 1.2);
  mu->setVal(mu_inject);
  mu->setConstant(true);

  for(int itoy = 0; itoy < ntoys; ++itoy) {
    if(itoy % 100 == 0) printf("    toy %i / %i\r", itoy, ntoys); fflush(stdout);

    // Generate toy dataset
    RooDataHist* toy_data = (RooDataHist*) full_pdf->generateBinned(*obs, RooFit::Extended(true));
    if(!toy_data) continue;

    // Fit to toy (mu >= 0: physical constraint)
    mu->setConstant(false);
    mu->setVal(std::max(mu_inject, 0.1));
    mu->setRange(0., mu_range_ts * 1.2);
    RooFitResult* r = full_pdf->fitTo(*toy_data,
                                       RooFit::Extended(true),
                                       RooFit::Save(true),
                                       RooFit::PrintLevel(-1),
                                       RooFit::Minimizer("Minuit2","Migrad"),
                                       RooFit::Strategy(1));
    if(r && r->status() == 0) {
      const double mu_val = mu->getVal();
      const double mu_err = mu->getError();
      mu_fit_vals.push_back(mu_val);
      mu_fit_errs.push_back(mu_err);
      h_mu.Fill(mu_val);
      if(mu_err > 0.) {
        const double pull = (mu_val - mu_inject) / mu_err;
        pulls.push_back(pull);
        h_pull.Fill(pull);
      }

      // Quick asymptotic UL for this toy
      // Use profile likelihood: UL = mu_hat + 1.28*sigma (one-sided 90% for mu_hat+UL)
      // More rigorously: use the Asimov dataset or rerun the inverter
      // For speed in toy studies we use the Wald approximation:
      // q_mu = 2*(NLL(mu) - NLL(mu_hat)), UL where q_mu = 1.642^2
      // Here we just store mu_hat + 1.645*sigma as an approximation
      const double ul_approx = std::max(0., mu_val) + 1.645 * mu_err;
      ul_vals.push_back(ul_approx);
      h_ul.Fill(ul_approx);
    }
    if(r) delete r;
    delete toy_data;

    // Reset nuisances to 0 between toys
    reset_nuisances(nuisanceParams);
  }
  printf("\n");

  // Summary statistics
  double mu_mean(0.), mu_rms(0.), pull_mean(0.), pull_rms(0.), ul_mean(0.);
  if(!mu_fit_vals.empty()) {
    for(auto v : mu_fit_vals) mu_mean += v;
    mu_mean /= mu_fit_vals.size();
    for(auto v : mu_fit_vals) mu_rms += (v-mu_mean)*(v-mu_mean);
    mu_rms = sqrt(mu_rms / mu_fit_vals.size());
  }
  if(!pulls.empty()) {
    for(auto v : pulls) pull_mean += v;
    pull_mean /= pulls.size();
    for(auto v : pulls) pull_rms += (v-pull_mean)*(v-pull_mean);
    pull_rms = sqrt(pull_rms / pulls.size());
  }
  if(!ul_vals.empty()) {
    for(auto v : ul_vals) ul_mean += v;
    ul_mean /= ul_vals.size();
  }

  printf("  Converged toys: %lu / %i\n"  , mu_fit_vals.size(), ntoys);
  printf("  Fitted mu: mean=%.3f  rms=%.3f  (inject=%.3f)\n", mu_mean, mu_rms, mu_inject);
  printf("  Pull:      mean=%.3f  rms=%.3f  (ideal: 0, 1)\n", pull_mean, pull_rms);
  printf("  Approx 90%% UL (mean over toys): %.4f events\n", ul_mean);

  // Plot 1: best-fit signal rate
  {
    TCanvas* c = new TCanvas("c_toy_mu","",800,600);
    h_mu.Fit("gaus","Q");
    h_mu.Draw();
    TLine lmu(mu_inject, 0., mu_inject, h_mu.GetMaximum());
    lmu.SetLineColor(kRed); lmu.SetLineStyle(kDashed); lmu.SetLineWidth(2); lmu.Draw();
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.032); lat.SetTextFont(42);
    lat.DrawLatex(0.14, 0.85, Form("#mu_{inject} = %.2f,  mean = %.2f,  rms = %.2f",
                                    mu_inject, mu_mean, mu_rms));
    c->SaveAs(Form("%s/toy_study_muhat_mu%g_%i.png", figdir.Data(), mu_inject, selection));
    delete c;
  }
  // Plot 2: pull distribution
  {
    TCanvas* c = new TCanvas("c_toy_pull","",800,600);
    h_pull.Fit("gaus","Q");
    h_pull.Draw();
    TLine lpu(0., 0., 0., h_pull.GetMaximum());
    lpu.SetLineColor(kRed); lpu.SetLineStyle(kDashed); lpu.SetLineWidth(2); lpu.Draw();
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.032); lat.SetTextFont(42);
    lat.DrawLatex(0.14, 0.85, Form("mean = %.2f,  rms = %.2f  (ideal: 0, 1)", pull_mean, pull_rms));
    c->SaveAs(Form("%s/toy_study_pull_mu%g_%i.png", figdir.Data(), mu_inject, selection));
    delete c;
  }
  // Plot 3: approximate 90% UL distribution
  {
    TCanvas* c = new TCanvas("c_toy_ul","",800,600);
    h_ul.Draw();
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.032); lat.SetTextFont(42);
    lat.DrawLatex(0.14, 0.85, Form("Mean approx UL = %.2f events", ul_mean));
    c->SaveAs(Form("%s/toy_study_ul_mu%g_%i.png", figdir.Data(), mu_inject, selection));
    delete c;
  }

  // Reset mu
  mu->setVal(0.);
  mu->setConstant(false);
  return 0;
}

//---------------------------------------------------------------------------
// Main entry point
//---------------------------------------------------------------------------
int stat_analysis(
    TString process   = "mumem",   // process: "mumem" or "mumep"
    int     selection = 20,        // selection number
    TString tag       = "",        // optional tag matching build_model.C
    int     analyses  = 0x1f,      // bitmask: 0x1=fit, 0x2=CLs_asym, 0x4=CLs_toys, 0x8=FC, 0x10=toys
    int     ntoys     = 500,       // toys for CLs/FC/pull studies
    double  mu_inject = 0.         // injected signal for toy study
) {
  process.ToLower();

  //-----------------------------------------------------------------
  // Output directory
  //-----------------------------------------------------------------
  TString figdir = Form("figures/%s%s", process.Data(), (tag != "") ? ("_"+tag).Data() : "");
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir.Data(), figdir.Data()));

  //-----------------------------------------------------------------
  // Load workspace
  //-----------------------------------------------------------------
  TString ws_file = Form("workspaces/workspace_%s_%i%s.root",
                          process.Data(), selection,
                          (tag != "") ? ("_"+tag).Data() : "");
  TFile* f = TFile::Open(ws_file, "READ");
  if(!f || f->IsZombie()) {
    printf("stat_analysis: ERROR – cannot open workspace file %s\n", ws_file.Data());
    return 1;
  }
  RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
  if(!ws) {
    printf("stat_analysis: ERROR – workspace object not found in %s\n", ws_file.Data());
    f->Close(); return 1;
  }
  printf("stat_analysis: Loaded workspace from %s\n", ws_file.Data());
  ws->Print("v");

  //-----------------------------------------------------------------
  // Retrieve data
  //-----------------------------------------------------------------
  RooAbsData* data = get_data_ws(ws);
  if(!data) {
    printf("stat_analysis: ERROR – dataset 'data_obs' not found\n");
    f->Close(); return 1;
  }
  printf("stat_analysis: Dataset has %.0f entries\n", data->sumEntries());

  //-----------------------------------------------------------------
  // Normalization information
  //-----------------------------------------------------------------
  RooRealVar* npot     = get_var(ws, "npot");
  RooRealVar* livetime = get_var(ws, "livetime");
  RooRealVar* nmuons   = get_var(ws, "nmuons");
  npot_     = (npot    ) ? npot    ->getVal() : -1.;
  livetime_ = (livetime) ? livetime->getVal() : -1.;
  nmuons_   = (nmuons  ) ? nmuons  ->getVal() : -1.;

  //-----------------------------------------------------------------
  // Signal parameter of interest
  //-----------------------------------------------------------------
  RooRealVar* sig_norm = get_var(ws, Form("%s_%i_signal_norm", process.Data(), selection));
  const double sig_norm_val = sig_norm ? sig_norm->getVal() : 1.;
  printf("stat_analysis: Nominal signal yield = %.4f events\n", sig_norm_val);

  RooRealVar* mu = new RooRealVar("mu", "Signal yield (events)", 0., 0., 20.);

  //-----------------------------------------------------------------
  // Build full model with systematic nuisances
  //-----------------------------------------------------------------
  std::vector<ShapeSysInfo_t> sys_infos;
  RooArgSet nuisanceParams;
  RooArgSet constraintPdfs;
  std::vector<RooAbsPdf*>  owned_pdfs;
  std::vector<RooAbsReal*> owned_reals;

  RooProdPdf* full_pdf = build_full_model(ws, process, selection, mu,
                                            sys_infos, nuisanceParams, constraintPdfs,
                                            owned_pdfs, owned_reals);
  if(!full_pdf) {
    printf("stat_analysis: ERROR – model construction failed!\n");
    f->Close(); return 1;
  }

  if(sys_infos.empty()) {
    printf("stat_analysis: No systematics found in workspace (using nominal model only)\n");
  } else {
    printf("stat_analysis: Found %lu systematic(s):\n", sys_infos.size());
    for(auto& s : sys_infos) printf("  - %s\n", s.name_.Data());
  }

  int status = 0;

  //-----------------------------------------------------------------
  // Run requested analyses
  //-----------------------------------------------------------------
  if(analyses & 0x1)
    status += do_likelihood_fit(ws, process, selection, tag, full_pdf, data, mu, nuisanceParams, figdir);

  if(analyses & 0x2)
    status += do_cls_asymptotic(ws, process, selection, tag, full_pdf, data, mu, nuisanceParams, figdir);

  if(analyses & 0x4)
    status += do_cls_toys(ws, process, selection, tag, full_pdf, data, mu, nuisanceParams, figdir, ntoys);

  if(analyses & 0x8)
    status += do_feldman_cousins(ws, process, selection, tag, full_pdf, data, mu, nuisanceParams, figdir, ntoys);

  if(analyses & 0x10)
    status += do_toy_study(ws, process, selection, tag, full_pdf, mu, nuisanceParams, figdir, ntoys, mu_inject);

  if(analyses & 0x20)
    status += do_fit_toys(ws, process, selection, tag, full_pdf, data, mu, nuisanceParams, figdir, 4*ntoys, mu_inject);

  //-----------------------------------------------------------------
  // Cleanup
  //-----------------------------------------------------------------
  for(auto* p : owned_pdfs)  delete p;
  for(auto* r : owned_reals) delete r;
  delete mu;
  f->Close();

  print_divider();
  printf("stat_analysis: Done. Status = %i\n", status);
  printf("  Output figures in: %s/\n", figdir.Data());
  return status;
}
