#ifndef __CONVANA_ANALYSIS_BACKGROUNDMODEL__
#define __CONVANA_ANALYSIS_BACKGROUNDMODEL__

#include "../tools/types.C"
#include "../defaults.C"
#include "../datasets.C"
#include "../physics.C"
#include "model_io_utils.C"
#include "RooStitchedPdf.cxx"
#include "../tools/utilities.C"

TString data_dataset_key_from_tag(TString tag) {
  tag.ToLower();
  if(tag.Contains("mds1d")) return "data_mds1d";
  if(tag.Contains("mds1f")) return "data_mds1f";
  if(tag.Contains("mds1g")) return "data_mds1g";
  if(tag.Contains("mds2a")) return "data_mds2a";
  if(tag.Contains("mds2b")) return "data_mds2b";
  if(tag.Contains("mds2c")) return "data_mds2c";
  if(tag.Contains("mds3c")) return "data_mds3c";
  return "";
}

//---------------------------------------------------------------------------------------------------------------------------
TH1* get_data_hist(const TString process, const int selection, TString tag = "") {
  const TString key = data_dataset_key_from_tag(tag);
  if(key == "") {
    cout << __func__ << ": No data dataset key mapped for tag " << tag.Data() << endl;
    return nullptr;
  }
  TH1* h = load_component_hist_from_dataset(key,
                                            selection,
                                            Form("%s_data%s", process.Data(), tag.Data()),
                                            -1,
                                            var_,
                                            false);
  return h;
}

//---------------------------------------------------------------------------------------------------------------------------
RooDataHist* get_data(RooRealVar& obs, const TString process, const int selection, TString tag = "") {
  TH1* data = get_data_hist(process, selection, tag);
  if(!data) return nullptr;
  const bool is_mumem = process == "mumem";
  const float xmin(is_mumem ? xmin_em_ : xmin_ep_), xmax(is_mumem ? xmax_em_ : xmax_ep_);
  data = trim_hist(data, xmin, xmax); // restrict to the relevant range
  RooDataHist* data_hist = new RooDataHist("data_obs", "Data histograms", obs, data);
  return data_hist;
}

//---------------------------------------------------------------------------------------------------------------------------
TH1* get_background_hist(const TString process, const int selection, const TString name = "bkg", const int isys = -1, TString var_name = var_) {
  return load_component_hist_from_dataset(process, selection, name, isys, var_name);
}

//---------------------------------------------------------------------------------------------------------------------------
TH1* get_background_hist_multi(const TString process,
                               const std::vector<int>& selections,
                               const TString name = "bkg",
                               const int isys = -1,
                               TString var_name = var_) {
  return load_component_hist_from_sets(process, selections, name, isys, var_name);
}

//---------------------------------------------------------------------------------------------------------------------------
pdf_info read_model(const TString name, const TString process, const int selection, TString tag = "", const int isys = -1) {

  if(tag != "") tag = "_" + tag;
  pdf_info res;
  const char* file_name = Form("workspaces/%s%s/%s_fit_%i%s.root", process.Data(), tag.Data(), name.Data(), selection, (isys > 0) ? Form("_sys_%i", isys) : "");
  TFile* f = TFile::Open(file_name, "READ");
  if(!f) {
    cout << __func__ << ": Unable to retrieve " << file_name << " workspace file!\n";
    return res;
  }
  RooWorkspace* workspace = (RooWorkspace*) f->Get("workspace");
  if(!workspace) {
    cout << __func__ << ": Unable to retrieve the workspace from file " << file_name << endl;
    return res;
  }

  RooAbsPdf* pdf = (RooAbsPdf*) workspace->pdf(Form("%s_%i_%s_pdf", process.Data(), selection, name.Data()));
  if(!pdf) {
    cout << __func__ << ": Unable to retrieve the PDF from file " << file_name << endl;
    return res;
  }

  RooRealVar* norm = (RooRealVar*) workspace->var(Form("%s_%i_%s_norm", process.Data(), selection, name.Data()));
  if(!norm) {
    cout << __func__ << ": Unable to retrieve the normalization from file " << file_name << endl;
    return res;
  }

  RooRealVar* obs = (RooRealVar*) workspace->var(Form("obs_%i", selection));
  if(!obs) {
    cout << __func__ << ": Unable to retrieve the observable from file " << file_name << endl;
    return res;
  }

  TH1* h = (TH1*) f->Get(Form("%s_%i_%s_hist", process.Data(), selection, name.Data()));
  if(!h) {
    cout << __func__ << ": Unable to retrieve the base histogram from file " << file_name << endl;
    return res;
  }
  h->SetDirectory(0);

  TH1* h_raw = (TH1*) f->Get(Form("%s_%i_%s_raw_hist", process.Data(), selection, name.Data()));
  if(h_raw) h_raw->SetDirectory(0);

  TH1* h_smoothed = (TH1*) f->Get(Form("%s_%i_%s_smoothed_hist", process.Data(), selection, name.Data()));
  if(h_smoothed) h_smoothed->SetDirectory(0);

  TH1* h_t0_raw = (TH1*) f->Get(Form("%s_%i_%s_t0_raw_hist", process.Data(), selection, name.Data()));
  if(h_t0_raw) h_t0_raw->SetDirectory(0);

  TH1* h_t0_smoothed = (TH1*) f->Get(Form("%s_%i_%s_t0_smoothed_hist", process.Data(), selection, name.Data()));
  if(h_t0_smoothed) h_t0_smoothed->SetDirectory(0);

  res.pdf_   = pdf;
  res.hist_  = h;
  res.raw_hist_ = h_raw;
  res.smoothed_hist_ = h_smoothed;
  res.t0_raw_hist_ = h_t0_raw;
  res.t0_smoothed_hist_ = h_t0_smoothed;
  res.obs_   = obs;
  res.norm_  = norm;
  res.rate_  = norm->getVal();
  res.name_  = name;
  res.sys_   = isys;
  set_style(name, res.title_, res.color_);

  return res;
}

//---------------------------------------------------------------------------------------------------------------------------
pdf_info get_dio_model(RooRealVar& obs, const TString process, const int selection, const bool freeze = true) {

  pdf_info res;

  // Convolve the theory PDF with the resolution function
  const char* name = Form("%s_%i_dio", process.Data(), selection);

  const int fit_version = 1; // 0: Convolve with response; 1: approx model

  RooAbsPdf* pdf = nullptr;
  if(fit_version == 0) {
    RooAbsPdf* resolution = nullptr;
    const int resolution_version = 0;

    if(resolution_version == 0) { // double-sided Crystal Ball
      RooRealVar* alpha1    = new RooRealVar(Form("%s_res_alpha1", name), "alpha1", 1.4, 0.1, 10.);
      RooRealVar* alpha2    = new RooRealVar(Form("%s_res_alpha2", name), "alpha2", 1.2, 0.1, 10.);
      RooRealVar* mean      = new RooRealVar(Form("%s_res_mean"  , name), "mean", 0.0, -0.3, 0.3);// mean->setConstant(true);
      RooRealVar* n1        = new RooRealVar(Form("%s_res_n1"    , name), "enne1", 7.3, 0.1, 30.);
      RooRealVar* n2        = new RooRealVar(Form("%s_res_n2"    , name), "enne2", 6.7, 0.1, 30.);
      RooRealVar* sigma     = new RooRealVar(Form("%s_res_sigma" , name), "sigma", 0.1, 0.001, 1.);
      resolution = new RooCrystalBall(Form("%s_res_pdf"  , name), "DIO resolution", obs, *mean, *sigma, *alpha1, *n1, *alpha2, *n2);

      if(selection == 20) {
        alpha1   ->setVal( 0.101064 ); // +/- 0.748468
        alpha2   ->setVal( 3.65888  ); // +/- 875.439
        mean     ->setVal( -0.299851); // +/- 0.0237112
        n1       ->setVal( 29.6829  ); // +/- 101.184
        n2       ->setVal( 24.8761  ); // +/- 1.53563e+07
        sigma    ->setVal( 0.0127344); // +/- 0.111337
      }

      alpha1 ->setConstant(freeze);
      alpha2 ->setConstant(freeze);
      mean   ->setConstant(freeze);
      n1     ->setConstant(freeze);
      n2     ->setConstant(freeze);
      sigma  ->setConstant(freeze);
      // } else if(resolution_version == 1) { // Landau + power-law tails
      //   RooRealVar* alpha1    = new RooRealVar(Form("%s_res_alpha1", name), "alpha1", 1.4, 0.1, 3.);
      //   RooRealVar* alpha2    = new RooRealVar(Form("%s_res_alpha2", name), "alpha2", 1.0, 0.1, 3.);
      //   RooRealVar* mean      = new RooRealVar(Form("%s_res_mean"  , name), "mean", -0.5, -1.5, 1.5);// mean->setConstant(true);
      //   RooRealVar* n1        = new RooRealVar(Form("%s_res_n1"    , name), "n1", 1.5, 0.1, 10.);
      //   RooRealVar* n2        = new RooRealVar(Form("%s_res_n2"    , name), "n2", 5.0, 0.1, 20.);
      //   RooRealVar* a         = new RooRealVar(Form("%s_res_a"     , name), "a", 0.3, 0.1, 1.);
      //   RooRealVar* b         = new RooRealVar(Form("%s_res_b"     , name), "b", 4.7, 0.1, 10.);
      //   resolution = new RooLandauCB(Form("%s_res_pdf"  , name), "DIO resolution", obs, *mean, *a, *b, *alpha1, *n1, *alpha2, *n2);

      //   alpha1 ->setConstant(freeze);
      //   alpha2 ->setConstant(freeze);
      //   mean   ->setConstant(freeze);
      //   n1     ->setConstant(freeze);
      //   n2     ->setConstant(freeze);
      //   a      ->setConstant(freeze);
      //   b      ->setConstant(freeze);
    }

    RooRealVar* ep        = new RooRealVar(Form("%s_theory_ep", name), "ep",  104.975   ); ep->setConstant(true);
    RooRealVar* a5        = new RooRealVar(Form("%s_theory_a5", name), "a5",  8.9       ); a5->setConstant(true);
    RooRealVar* a6        = new RooRealVar(Form("%s_theory_a6", name), "a6",  1.17169   ); a6->setConstant(true);
    RooRealVar* a7        = new RooRealVar(Form("%s_theory_a7", name), "a7", -1.06599e-2); a7->setConstant(true);
    RooRealVar* a8        = new RooRealVar(Form("%s_theory_a8", name), "a8",  8.14251e-3); a8->setConstant(true);
    RooRealVar* mmu       = new RooRealVar(Form("%s_theory_mmu", name), "mmu",  105.658 ); mmu->setConstant(true);
    RooRealVar* fl        = new RooRealVar(Form("%s_theory_fl", name), "fl",  0.023343129 ); fl->setConstant(true);

    RooAbsPdf* theory     = new RooGenericPdf(Form("%s_theory", name), "(@0 < @1)*(@2*pow(@1-@0,5) + @3*pow(@1-@0,6) + @4*pow(@1-@0,7) + @5*pow(@1-@0,8))*pow(max(0., (@1-@0))/@6,@7)",
                                              RooArgList(obs, *ep, *a5, *a6, *a7, *a8, *mmu, *fl));
    // RooAbsPdf* theory     = new RooGenericPdf(Form("%s_theory", name), "(@0 < @1)*(@2*pow(@1-@0,5) + @3*pow(@1-@0,6) + @4*pow(@1-@0,7) + @5*pow(@1-@0,8))",
    //                                           RooArgList(obs, *ep, *a5, *a6, *a7, *a8));

    // obs.setBinning(RooBinning(10000, obs.getMin(), obs.getMax()), "cache"); // for discrete convolution
    // auto pdf        = new RooFFTConvPdf(Form("%s_pdf", name), "DIO PDF", obs, *theory, *resolution, 2);
    // pdf->setBufferFraction(1.);
    // pdf->setBufferStrategy(RooFFTConvPdf::Flat); //Extend, Flat, or Mirror
    pdf        = new RooNumConvPdf(Form("%s_pdf", name), "DIO PDF", obs, *theory, *resolution);
  } else if(fit_version == 1) {
    RooRealVar* x0     = new RooRealVar(Form("%s_x0", name), "Low Threshold Edge", 96.0, 80., 96.);
    RooRealVar* alpha  = new RooRealVar(Form("%s_alpha", name), "Low-edge turn-on power", 1.915, 0.1, 4.0);
    RooRealVar* beta   = new RooRealVar(Form("%s_beta", name), "Bulk curvature power", 0.457, 0.01, 20.0);
    RooRealVar* lambda = new RooRealVar(Form("%s_lambda", name), "Tail exponential decay", 1.932, 0.01, 5.0);
    pdf = new RooGenericPdf(Form("%s_pdf", name), "Gamma-Poly Hybrid",
                                 Form("(pow(max(0., %s - %s), %s) * pow(120.0 - %s, %s) * exp(-%s * %s))",
                                      obs.GetName(), x0->GetName(), alpha->GetName(),
                                      obs.GetName(), beta->GetName(),
                                      lambda->GetName(), obs.GetName()),
                                 RooArgList(obs, *x0, *alpha, *beta, *lambda));
  }

  const double rate_per_run1 = 10000.; //FIXME

  res.pdf_   = pdf;
  res.rate_  = rate_per_run1;
  res.color_ = kRed;
  res.name_  = "dio";
  res.title_ = "DIO";

  return res;
}

//---------------------------------------------------------------------------------------------------------------------------
pdf_info get_cosmic_model(RooRealVar& obs, TString process, int selection, const bool freeze = true) {

  pdf_info res;

  // Simply approximate with a flat distribution
  const char* name = Form("%s_%i_cosmic", process.Data(), selection);

  RooRealVar* p0 = new RooRealVar(Form("%s_p0", name), "p0", 0.1, -1., 1.);
  RooRealVar* p1 = new RooRealVar(Form("%s_p1", name), "p1", 0.1, -1., 1.);

  auto pdf = new RooChebychev(Form("%s_pdf", name), "Cosmic background", obs, RooArgList(*p0, *p1));
  // auto pdf = new RooUniform(Form("%s_pdf", name), "Cosmic background", obs);

  const double rate_per_run1 = 0.048 * (obs.getMax() - obs.getMin())/1.3;

  if(selection == 20) {
    p0->setVal( -0.0491244 ); // +/- 0.0106895
    p1->setVal( -0.00400783); // +/- 0.0103816
  } else if(selection == 40) {
    p0->setVal( -0.00529238); // +/- 0.0202569
    p1->setVal( 0.00156699 ); // +/- 0.0196947
  }
  p0->setConstant(freeze);
  p1->setConstant(freeze);

  res.pdf_   = pdf;
  res.rate_  = rate_per_run1;
  res.color_ = kOrange;
  res.name_  = "cosmic";
  res.title_ = "Cosmic ray";

  return res;
}

//---------------------------------------------------------------------------------------------------------------------------
pdf_info get_pbar_model(RooRealVar& obs, TString process, int selection, const bool freeze = true) {

  pdf_info res;

  const char* name = Form("%s_%i_pbar", process.Data(), selection);

  RooRealVar* p0 = new RooRealVar(Form("%s_p0", name), "p0", -0.05 , -1., 1.);
  RooRealVar* p1 = new RooRealVar(Form("%s_p1", name), "p1", -0.004, -1., 1.);

  auto pdf = new RooChebychev(Form("%s_pdf", name), "Antiproton background", obs, RooArgList(*p0));

  const double rate_per_run1 = 0.01 * (obs.getMax() - obs.getMin())/1.3;

  p0->setConstant(freeze);
  p1->setConstant(freeze);

  res.pdf_   = pdf;
  res.rate_  = rate_per_run1;
  res.color_ = kGreen-6;
  res.name_  = "pbar";
  res.title_ = "Antiproton";

  return res;
}

//---------------------------------------------------------------------------------------------------------------------------
pdf_info get_rpc_ext_model(RooRealVar& obs, TString process, int selection, const bool freeze = true) {

  pdf_info res;

  const char* name = Form("%s_%i_rpc_ext", process.Data(), selection);

  RooRealVar* p0 = new RooRealVar(Form("%s_p0", name), "p0", 0.1, -1., 1.);
  RooRealVar* p1 = new RooRealVar(Form("%s_p1", name), "p1", 0.1, -1., 1.);

  auto pdf = new RooChebychev(Form("%s_pdf", name), "RPC (external) background", obs, RooArgList(*p0));
  // auto pdf = new RooChebychev(Form("%s_pdf", name), "RPC (external) background", obs, RooArgList(*p0, *p1));

  const double rate_per_run1 = 0.01 * (obs.getMax() - obs.getMin())/1.3;

  if(selection == 20) {
    p0->setVal( -0.0491244 ); // +/- 0.0106895
    p1->setVal( -0.00400783); // +/- 0.0103816
  } else if(selection == 40) {
    p0->setVal( -0.00529238); // +/- 0.0202569
    p1->setVal( 0.00156699 ); // +/- 0.0196947
  }
  p0->setConstant(freeze);
  p1->setConstant(freeze);

  res.pdf_   = pdf;
  res.rate_  = rate_per_run1;
  res.color_ = kMagenta;
  res.name_  = "rpc_ext";
  res.title_ = "RPC (external)";

  return res;
}

//---------------------------------------------------------------------------------------------------------------------------
pdf_info get_rpc_int_model(RooRealVar& obs, TString process, int selection, const bool freeze = true) {

  pdf_info res;

  const char* name = Form("%s_%i_rpc_int", process.Data(), selection);

  RooRealVar* p0 = new RooRealVar(Form("%s_p0", name), "p0", 0.1, -1., 1.);
  RooRealVar* p1 = new RooRealVar(Form("%s_p1", name), "p1", 0.1, -1., 1.);

  auto pdf = new RooChebychev(Form("%s_pdf", name), "RPC (internal) background", obs, RooArgList(*p0, *p1));

  const double rate_per_run1 = 0.01 * (obs.getMax() - obs.getMin())/1.3;

  if(selection == 20) {
    p0->setVal( -0.0491244 ); // +/- 0.0106895
    p1->setVal( -0.00400783); // +/- 0.0103816
  } else if(selection == 40) {
    p0->setVal( -0.00529238); // +/- 0.0202569
    p1->setVal( 0.00156699 ); // +/- 0.0196947
  }
  p0->setConstant(freeze);
  p1->setConstant(freeze);

  res.pdf_   = pdf;
  res.rate_  = rate_per_run1;
  res.color_ = kMagenta;
  res.name_  = "rpc_int";
  res.title_ = "RPC (internal)";

  return res;
}

//---------------------------------------------------------------------------------------------------------------------------
pdf_info get_rmc_ext_model(RooRealVar& obs, TString process, int selection, const bool freeze = true) {

  pdf_info res;

  const char* name = Form("%s_%i_rmc_ext", process.Data(), selection);

  RooRealVar* x0     = new RooRealVar(Form("%s_x0", name), "Low Threshold Edge", 96.0, 80., 96.);
  RooRealVar* alpha  = new RooRealVar(Form("%s_alpha", name), "Low-edge turn-on power", 1.915, 0.1, 4.0);
  RooRealVar* beta   = new RooRealVar(Form("%s_beta", name), "Bulk curvature power", 0.457, 0.01, 20.0);
  RooRealVar* lambda = new RooRealVar(Form("%s_lambda", name), "Tail exponential decay", 1.932, 0.01, 5.0);
  auto pdf = new RooGenericPdf(Form("%s_pdf", name), "Gamma-Poly Hybrid",
                               Form("(pow(max(0., %s - %s), %s) * pow(120.0 - %s, %s) * exp(-%s * %s))",
                                    obs.GetName(), x0->GetName(), alpha->GetName(),
                                    obs.GetName(), beta->GetName(),
                                    lambda->GetName(), obs.GetName()),
                               RooArgList(obs, *x0, *alpha, *beta, *lambda));


  res.pdf_   = pdf;
  res.rate_  = 50.; // rough starting point
  res.color_ = kAtlantic+2;
  res.name_  = "rmc_ext";
  res.title_ = "RMC (external)";

  return res;
}

//---------------------------------------------------------------------------------------------------------------------------
pdf_info get_rmc_int_model(RooRealVar& obs, TString process, int selection, const bool freeze = true) {

  pdf_info res;

  const char* name = Form("%s_%i_rmc_int", process.Data(), selection);

  RooRealVar* x0     = new RooRealVar(Form("%s_x0", name), "Low Threshold Edge", 96.0, 80., 96.);
  RooRealVar* alpha  = new RooRealVar(Form("%s_alpha", name), "Low-edge turn-on power", 1.915, 0.1, 4.0);
  RooRealVar* beta   = new RooRealVar(Form("%s_beta", name), "Bulk curvature power", 0.457, 0.01, 20.0);
  RooRealVar* lambda = new RooRealVar(Form("%s_lambda", name), "Tail exponential decay", 1.932, 0.01, 5.0);
  auto pdf = new RooGenericPdf(Form("%s_pdf", name), "Gamma-Poly Hybrid",
                               Form("(pow(max(0., %s - %s), %s) * pow(120.0 - %s, %s) * exp(-%s * %s))",
                                    obs.GetName(), x0->GetName(), alpha->GetName(),
                                    obs.GetName(), beta->GetName(),
                                    lambda->GetName(), obs.GetName()),
                               RooArgList(obs, *x0, *alpha, *beta, *lambda));


  res.pdf_   = pdf;
  res.rate_  = 50.; // rough starting point
  res.color_ = kAtlantic;
  res.name_  = "rmc_int";
  res.title_ = "RMC (internal)";

  return res;
}

//---------------------------------------------------------------------------------------------------------------------------
pdf_infos get_background_model(RooRealVar& obs, TString process, const int selection, const TString tag = "", const bool hist_pdf = false) {

  // to keep separate labels
  combine_rpc_ = false;
  combine_rmc_ = false;

  pdf_infos pdfs;
  const bool is_mumem = process == "mumem";
  if(is_mumem)
    pdfs.push_back(read_model("dio", process, selection, tag));
  pdfs.push_back(read_model("cosmic", process, selection, tag));
  pdfs.push_back(read_model("rpc_ext", process, selection, tag));
  pdfs.push_back(read_model("rpc_int", process, selection, tag));
  if(!use_evtana_) {
    pdfs.push_back(read_model("pbar", process, selection, tag));
  } else {
    pdfs.push_back(read_model("rmc_ext", process, selection, tag));
    pdfs.push_back(read_model("rmc_int", process, selection, tag));
  }

  return pdfs;
}

#endif
