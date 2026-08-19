// Shared helpers for exporting fitted PDFs/histograms into workspaces
#ifndef __CONVANA_ANALYSIS_FITWORKSPACEUTILS__
#define __CONVANA_ANALYSIS_FITWORKSPACEUTILS__

#include "../tools/utilities.C"

int enforce_uniform_if_sparse(TH1* h,
                              RooRealVar& obs,
                              RooAbsPdf*& pdf,
                              const TString pdf_name,
                              const TString pdf_title,
                              const double min_filled_fraction = 0.75);

TString normalize_pdf_type(TString pdf_type) {
  pdf_type.ToLower();
  pdf_type.ReplaceAll(" ", "");
  return pdf_type;
}

TString default_pdf_type_for_component(const TString& component) {
  TString c = component;
  c.ToLower();
  if(c == "signal")  return "auto";
  if(c == "dio")     return "auto";
  if(c == "cosmic")  return "auto";
  if(c == "rpc_ext") return "auto";
  if(c == "rpc_int") return "auto";
  if(c == "pbar")    return "auto";
  if(c == "rmc_ext") return "auto";
  if(c == "rmc_int") return "auto";
  return "auto";
}

TString default_tail_model_for_component(const TString& component) {
  TString c = component;
  c.ToLower();
  if(c == "signal")  return "exp";
  if(c == "dio")     return "convolution";
  if(c == "rmc_ext") return "exp";
  if(c == "rmc_int") return "exp";
  return "none";
}

TString resolve_pdf_type(const TString& component, TString requested_pdf_type) {
  requested_pdf_type = normalize_pdf_type(requested_pdf_type);
  if(requested_pdf_type == "" || requested_pdf_type == "default") {
    return default_pdf_type_for_component(component);
  }
  return requested_pdf_type;
}

TString resolve_tail_model(const TString& component, TString requested_tail_model) {
  requested_tail_model = normalize_pdf_type(requested_tail_model);
  if(requested_tail_model == "" || requested_tail_model == "default") {
    return default_tail_model_for_component(component);
  }
  return requested_tail_model;
}

std::vector<int> default_control_region_sets(const TString& component, const int selection) {
  TString c = component;
  c.ToLower();
  if(c == "cosmic") return {selection + 1000};
  return {};
}

int parse_poly_degree(const TString& pdf_type, const int default_degree = 1) {
  TString lower = normalize_pdf_type(pdf_type);
  if(!lower.BeginsWith("poly")) return default_degree;
  TString degree_str = lower;
  degree_str.ReplaceAll("poly", "");
  if(degree_str.Length() == 0) return default_degree;
  const int degree = degree_str.Atoi();
  return std::max(0, degree);
}

RooAbsPdf* make_poly_pdf(RooRealVar& obs,
                         const TString& pdf_name,
                         const TString& pdf_title,
                         const int degree) {
  RooArgList coeffs;
  for(int i = 0; i <= degree; ++i) {
    auto coeff = new RooRealVar(Form("%s_p%i", pdf_name.Data(), i),
                                Form("poly coeff %i", i),
                                0.0, -1.0, 1.0);
    coeffs.add(*coeff);
  }
  return new RooChebychev(pdf_name, pdf_title, obs, coeffs);
}

RooAbsPdf* make_crystal_ball_pdf(RooRealVar& obs,
                                 const TString& pdf_name,
                                 const TString& pdf_title) {
  const double xmin = obs.getMin();
  const double xmax = obs.getMax();
  const double center = 0.5*(xmin + xmax);
  const double width = std::max(1.e-3, 0.1*(xmax - xmin));

  auto mean   = new RooRealVar(Form("%s_cb_mean", pdf_name.Data()), "CB mean", center, xmin, xmax);
  auto sigma  = new RooRealVar(Form("%s_cb_sigma", pdf_name.Data()), "CB sigma", width, 1.e-3, std::max(0.2, xmax - xmin));
  auto alphaL = new RooRealVar(Form("%s_cb_alphaL", pdf_name.Data()), "CB alphaL", 1.5, 0.1, 10.0);
  auto nL     = new RooRealVar(Form("%s_cb_nL", pdf_name.Data()), "CB nL", 3.0, 0.1, 50.0);
  auto alphaR = new RooRealVar(Form("%s_cb_alphaR", pdf_name.Data()), "CB alphaR", 1.5, 0.1, 10.0);
  auto nR     = new RooRealVar(Form("%s_cb_nR", pdf_name.Data()), "CB nR", 3.0, 0.1, 50.0);
  return new RooCrystalBall(pdf_name, pdf_title, obs, *mean, *sigma, *alphaL, *nL, *alphaR, *nR);
}

RooAbsPdf* choose_pdf_model(const TString& requested_pdf_type,
                            const bool default_hist_pdf,
                            RooRealVar& obs,
                            RooDataHist& data_hist,
                            RooAbsPdf* analytic_pdf,
                            const TString& pdf_name,
                            const TString& pdf_title,
                            TH1* sparse_hist = nullptr,
                            const double sparse_uniform_threshold = 0.75) {
  const TString mode = normalize_pdf_type(requested_pdf_type);
  RooAbsPdf* pdf = nullptr;

  if(mode == "auto") {
    pdf = (default_hist_pdf || !analytic_pdf)
      ? static_cast<RooAbsPdf*>(new RooHistPdf(pdf_name, pdf_title, obs, data_hist))
      : analytic_pdf;
  } else if(mode == "hist" || mode == "histogram") {
    pdf = new RooHistPdf(pdf_name, pdf_title, obs, data_hist);
  } else if(mode == "uniform") {
    pdf = new RooUniform(pdf_name, pdf_title, obs);
  } else if(mode.BeginsWith("poly")) {
    pdf = make_poly_pdf(obs, pdf_name, pdf_title, parse_poly_degree(mode, 1));
  } else if(mode == "cb" || mode == "crystalball") {
    pdf = make_crystal_ball_pdf(obs, pdf_name, pdf_title);
  } else if(mode == "analytic" || mode == "model") {
    pdf = analytic_pdf;
  }

  if((mode == "analytic" || mode == "model") && !pdf) {
    cout << __func__ << ": Analytic model requested for " << pdf_name.Data()
         << " but no analytic model is available, using histogram PDF instead." << endl;
    pdf = new RooHistPdf(pdf_name, pdf_title, obs, data_hist);
  }

  if(!pdf) {
    cout << __func__ << ": Unknown PDF type \"" << requested_pdf_type.Data()
         << "\", using auto mode" << endl;
    pdf = (default_hist_pdf)
      ? static_cast<RooAbsPdf*>(new RooHistPdf(pdf_name, pdf_title, obs, data_hist))
      : analytic_pdf;
  }

  if(mode == "auto" && sparse_hist) {
    enforce_uniform_if_sparse(sparse_hist, obs, pdf, pdf_name, pdf_title, sparse_uniform_threshold);
  }
  return pdf;
}

bool should_fit_pdf(const TString& requested_pdf_type, const bool default_hist_pdf) {
  const TString mode = normalize_pdf_type(requested_pdf_type);
  if(mode == "uniform") return false;
  if(mode == "hist" || mode == "histogram") return false;
  if(mode == "auto" && default_hist_pdf) return false;
  return true;
}

int run_component_fit(RooAbsPdf* pdf,
                      RooDataHist& data_hist,
                      const TString& requested_pdf_type,
                      const bool default_hist_pdf,
                      const bool use_sumw2 = true,
                      const TString& fit_range = "") {
  if(!pdf) return 1;
  if(pdf->InheritsFrom("RooHistPdf") || pdf->InheritsFrom("RooUniform")) return 0;
  if(!should_fit_pdf(requested_pdf_type, default_hist_pdf)) return 0;

  if(fit_range != "") {
    if(use_sumw2) pdf->fitTo(data_hist, RooFit::Range(fit_range.Data()), RooFit::SumW2Error(true));
    else          pdf->fitTo(data_hist, RooFit::Range(fit_range.Data()));
  } else {
    if(use_sumw2) pdf->fitTo(data_hist, RooFit::SumW2Error(true));
    else          pdf->fitTo(data_hist);
  }
  return 0;
}

void smooth_tail_from_reference(TH1* h,
                                const TH1* reference,
                                const double fit_xmin,
                                const double xmax,
                                const double rel_err = 0.2) {
  if(!h || !reference) return;
  for(int ibin = h->FindBin(fit_xmin + 1.e-6); ibin <= h->FindBin(xmax - 1.e-6); ++ibin) {
    const double x = h->GetBinCenter(ibin);
    const double val = reference->Interpolate(x) * h->GetBinWidth(ibin) / reference->GetBinWidth(ibin);
    h->SetBinContent(ibin, val);
    h->SetBinError(ibin, std::max(0., rel_err*val));
  }
}

int smooth_right_tail(TH1* h,
                      TH1* h_orig,
                      const TString& fig_path,
                      const TString& model,
                      const double fit_xmin,
                      const double fit_xmax,
                      const int nsmooth = 0,
                      const double rel_err = 0.2,
                      const bool use_weighted_fit = true,
                      const std::vector<double>& init_params = {}) {
  if(!h) return 1;
  if(nsmooth > 0) h->Smooth(nsmooth);

  TString lower = normalize_pdf_type(model);
  if(lower == "none" || lower == "") return 0;

  TF1* tail_func = nullptr;
  if(lower == "exp") {
    tail_func = new TF1("tail_func", "exp([0] + [1]*x)", fit_xmin, fit_xmax);
  } else if(lower == "power") {
    tail_func = new TF1("tail_func", "pow(max(1.e-9, (x-[0])/[1]), [2])", fit_xmin, fit_xmax);
  } else {
    cout << __func__ << ": Unknown smoothing model \"" << model.Data() << "\"" << endl;
    return 2;
  }

  for(size_t i = 0; i < init_params.size(); ++i) tail_func->SetParameter(i, init_params[i]);

  const TString fit_opt = (use_weighted_fit) ? "wR" : "R";
  h->Fit(tail_func, fit_opt.Data());
  for(int ibin = h->FindBin(fit_xmin + 1.e-6); ibin <= h->FindBin(fit_xmax - 1.e-6); ++ibin) {
    const double x = h->GetBinCenter(ibin);
    const double val = tail_func->Eval(x);
    h->SetBinContent(ibin, val);
    h->SetBinError(ibin, std::max(0., rel_err*val));
  }

  if(fig_path != "") {
    TCanvas c;
    if(h_orig) {
      h_orig->Draw("E1");
      h->SetLineColor(kRed);
      h->Draw("E1 same");
    } else {
      h->Draw("E1");
    }
    tail_func->Draw("same");
    c.SetLogy();
    c.SaveAs(fig_path.Data());
  }
  return 0;
}

TH1* convolve_with_resolution(TH1* theory, TH1* resolution) {
  TH1* result = (TH1*) theory->Clone("result");
  result->Reset();
  const double estep = 0.01;
  const double e_low  = result->GetXaxis()->GetBinLowEdge(1);
  const double e_high = result->GetXaxis()->GetBinUpEdge(result->GetNbinsX());
  const double de_low  = resolution->GetXaxis()->GetBinLowEdge(1);
  const double de_high = resolution->GetXaxis()->GetBinUpEdge(resolution->GetNbinsX());
  for(double energy = e_low; energy <= e_high; energy += estep) {
    const double p_theory = estep*theory->Interpolate(energy);
    for(double de = de_low; de <= de_high; de += estep) {
      const double p_reco = estep*resolution->Interpolate(de);
      result->Fill(energy + de, p_theory*p_reco);
    }
  }
  if(result->Integral() > 0.) result->Scale(1./result->Integral());
  return result;
}

TH1* convolve_with_resolution(TH1* theory, TF1* resolution) {
  TH1* result = (TH1*) theory->Clone("result");
  result->Reset();
  const double estep = 0.01;
  const double e_low  = result->GetXaxis()->GetBinLowEdge(1);
  const double e_high = result->GetXaxis()->GetBinUpEdge(result->GetNbinsX());
  const double de_low  = resolution->GetXmin();
  const double de_high = resolution->GetXmax();
  for(double energy = e_low; energy <= e_high; energy += estep) {
    const double p_theory = estep*theory->Interpolate(energy);
    for(double de = de_low; de <= de_high; de += estep) {
      const double p_reco = estep*resolution->Eval(de);
      result->Fill(energy + de, p_theory*p_reco);
    }
  }
  if(result->Integral() > 0.) result->Scale(1./result->Integral());
  return result;
}

int save_fit_plot_pair(TCanvas* c,
                       RooPlot* frame,
                       const TString out_base,
                       const double log_ymin = -1.,
                       const double log_ymax = -1.) {
  if(!c || !frame) return 1;

  c->SaveAs(Form("%s.png", out_base.Data()));

  if(log_ymin >= 0. && log_ymax > log_ymin) {
    frame->GetYaxis()->SetRangeUser(log_ymin, log_ymax);
  }

  TPad* pad1 = (TPad*) c->GetPrimitive("pad1");
  if(pad1) pad1->SetLogy();
  c->SaveAs(Form("%s_log.png", out_base.Data()));

  delete frame;
  delete c;
  return 0;
}

int enforce_uniform_if_sparse(TH1* h,
                              RooRealVar& obs,
                              RooAbsPdf*& pdf,
                              const TString pdf_name,
                              const TString pdf_title,
                              const double min_filled_fraction = 0.75) {
  if(!h || !pdf) return -1;

  int filled = 0;
  const int nbins = h->GetNbinsX();
  for(int bin = 1; bin <= nbins; ++bin) {
    if(h->GetBinContent(bin) > 0.) ++filled;
  }

  if(filled < min_filled_fraction*nbins) {
    cout << "Only " << filled << " / " << nbins << " bins have content --> using a uniform model!\n";
    delete pdf;
    pdf = new RooUniform(pdf_name, pdf_title, obs);
  }
  return filled;
}

int save_fit_workspace(TString process,
                       int selection,
                       TString tag,
                       TString component,
                       TString component_title,
                       RooAbsPdf*& pdf,
                       RooRealVar& obs,
                       RooRealVar& norm,
                       const bool hist_pdfs,
                       TH1* raw_hist = nullptr,
                       TH1* normalized_hist = nullptr,
                       TH1* smoothed_hist = nullptr,
                       TH1* t0_raw_hist = nullptr,
                       TH1* t0_normalized_hist = nullptr,
                       TH1* t0_smoothed_hist = nullptr,
                       const TString out_suffix = "") {
  const char* hist_name = Form("%s_%i_%s_hist", process.Data(), selection, component.Data());
  auto h_fit = pdf->createHistogram(hist_name, obs);
  h_fit->SetName(hist_name);

  if(hist_pdfs) {
    pdf->SetName("tmp_pdf");
    RooDataHist* fit_data = new RooDataHist(Form("%s_%i_%s_data_hist", process.Data(), selection, component.Data()),
                                            Form("%s data hist", component_title.Data()),
                                            obs, h_fit);
    pdf = new RooHistPdf(Form("%s_%i_%s_pdf", process.Data(), selection, component.Data()),
                         Form("%s PDF", component_title.Data()),
                         obs, *fit_data);
  }

  const char* fitdir = Form("workspaces/%s%s", process.Data(), (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", fitdir, fitdir));
  TFile* fout = new TFile(Form("%s/%s_fit_%i%s.root", fitdir, component.Data(), selection, out_suffix.Data()), "RECREATE");
  RooWorkspace ws("workspace", "workspace");

  auto write_hist_to_workspace = [&](TH1* src,
                                     const TString& hist_name,
                                     const TString& data_name,
                                     const TString& data_title) {
    if(!src) return;
    TH1* h_copy = (TH1*) src->Clone(hist_name.Data());
    h_copy->SetName(hist_name.Data());
    RooDataHist data(data_name.Data(), data_title.Data(), obs, h_copy);
    ws.import(data);
    h_copy->Write();
  };

  ws.import(*pdf);
  ws.import(norm);
  write_hist_to_workspace(raw_hist,
                          Form("%s_%i_%s_raw_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_raw_data_hist", process.Data(), selection, component.Data()),
                          Form("%s raw data hist", component_title.Data()));
  write_hist_to_workspace(normalized_hist,
                          Form("%s_%i_%s_normalized_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_normalized_data_hist", process.Data(), selection, component.Data()),
                          Form("%s normalized data hist", component_title.Data()));
  write_hist_to_workspace(smoothed_hist,
                          Form("%s_%i_%s_smoothed_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_smoothed_data_hist", process.Data(), selection, component.Data()),
                          Form("%s smoothed data hist", component_title.Data()));
  write_hist_to_workspace(t0_raw_hist,
                          Form("%s_%i_%s_t0_raw_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_t0_raw_data_hist", process.Data(), selection, component.Data()),
                          Form("%s t0 raw data hist", component_title.Data()));
  write_hist_to_workspace(t0_normalized_hist,
                          Form("%s_%i_%s_t0_normalized_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_t0_normalized_data_hist", process.Data(), selection, component.Data()),
                          Form("%s t0 normalized data hist", component_title.Data()));
  write_hist_to_workspace(t0_smoothed_hist,
                          Form("%s_%i_%s_t0_smoothed_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_t0_smoothed_data_hist", process.Data(), selection, component.Data()),
                          Form("%s t0 smoothed data hist", component_title.Data()));
  ws.Write();
  h_fit->Write();
  fout->Close();

  print_pdf(pdf);
  return 0;
}

int save_fit_workspace_with_hist(TString process,
                                 int selection,
                                 TString tag,
                                 TString component,
                                 RooAbsPdf* pdf,
                                 RooRealVar& obs,
                                 RooRealVar& norm,
                                 TH1* hist,
                                 TH1* raw_hist = nullptr,
                                 TH1* normalized_hist = nullptr,
                                 TH1* smoothed_hist = nullptr,
                                 TH1* t0_raw_hist = nullptr,
                                 TH1* t0_normalized_hist = nullptr,
                                 TH1* t0_smoothed_hist = nullptr,
                                 const TString out_suffix = "") {
  if(!pdf || !hist) return 1;

  const char* hist_name = Form("%s_%i_%s_hist", process.Data(), selection, component.Data());
  auto h_fit = (TH1*) hist->Clone(hist_name);
  h_fit->SetName(hist_name);

  pdf->SetName("tmp_pdf");
  RooDataHist fit_data(Form("%s_%i_%s_data_hist", process.Data(), selection, component.Data()),
                       Form("%s data hist", component.Data()),
                       obs, h_fit);
  RooHistPdf fit_pdf(Form("%s_%i_%s_pdf", process.Data(), selection, component.Data()),
                     Form("%s PDF", component.Data()),
                     obs, fit_data);

  const char* fitdir = Form("workspaces/%s%s", process.Data(), (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", fitdir, fitdir));
  TFile* fout = new TFile(Form("%s/%s_fit_%i%s.root", fitdir, component.Data(), selection, out_suffix.Data()), "RECREATE");
  RooWorkspace ws("workspace", "workspace");

  auto write_hist_to_workspace = [&](TH1* src,
                                     const TString& hist_name,
                                     const TString& data_name,
                                     const TString& data_title) {
    if(!src) return;
    TH1* h_copy = (TH1*) src->Clone(hist_name.Data());
    h_copy->SetName(hist_name.Data());
    RooDataHist data(data_name.Data(), data_title.Data(), obs, h_copy);
    ws.import(data);
    h_copy->Write();
  };

  ws.import(fit_pdf);
  ws.import(norm);
  write_hist_to_workspace(raw_hist,
                          Form("%s_%i_%s_raw_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_raw_data_hist", process.Data(), selection, component.Data()),
                          Form("%s raw data hist", component.Data()));
  write_hist_to_workspace(normalized_hist,
                          Form("%s_%i_%s_normalized_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_normalized_data_hist", process.Data(), selection, component.Data()),
                          Form("%s normalized data hist", component.Data()));
  write_hist_to_workspace(smoothed_hist,
                          Form("%s_%i_%s_smoothed_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_smoothed_data_hist", process.Data(), selection, component.Data()),
                          Form("%s smoothed data hist", component.Data()));
  write_hist_to_workspace(t0_raw_hist,
                          Form("%s_%i_%s_t0_raw_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_t0_raw_data_hist", process.Data(), selection, component.Data()),
                          Form("%s t0 raw data hist", component.Data()));
  write_hist_to_workspace(t0_normalized_hist,
                          Form("%s_%i_%s_t0_normalized_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_t0_normalized_data_hist", process.Data(), selection, component.Data()),
                          Form("%s t0 normalized data hist", component.Data()));
  write_hist_to_workspace(t0_smoothed_hist,
                          Form("%s_%i_%s_t0_smoothed_hist", process.Data(), selection, component.Data()),
                          Form("%s_%i_%s_t0_smoothed_data_hist", process.Data(), selection, component.Data()),
                          Form("%s t0 smoothed data hist", component.Data()));
  ws.Write();
  h_fit->Write();
  fout->Close();

  print_pdf(&fit_pdf);
  return 0;
}

#endif
