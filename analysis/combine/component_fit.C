#ifndef __CONVANA_ANALYSIS_COMPONENT_FIT__
#define __CONVANA_ANALYSIS_COMPONENT_FIT__

#include "../tools/utilities.C"
#include "../defaults.C"
#include "../physics.C"
#include "background_model.C"
#include "signal_model.C"
#include "fit_workspace_utils.C"

TString component_dataset_key(const TString& component, const TString& process) {
  TString c = component;
  c.ToLower();
  if(c == "signal") return process;
  if(c == "rmc_ext") return "rmc_ext_0n";
  return c;
}

RooAbsPdf* component_analytic_pdf(const TString& component,
                                  RooRealVar& obs,
                                  const TString& process,
                                  const int selection) {
  TString c = component;
  c.ToLower();
  if(c == "signal")  return get_signal_model (obs, process, selection, false).pdf_;
  if(c == "dio")     return get_dio_model    (obs, process, selection, false).pdf_;
  if(c == "cosmic")  return get_cosmic_model (obs, process, selection, false).pdf_;
  if(c == "pbar")    return get_pbar_model   (obs, process, selection, false).pdf_;
  if(c == "rpc_ext") return get_rpc_ext_model(obs, process, selection, false).pdf_;
  if(c == "rpc_int") return get_rpc_int_model(obs, process, selection, false).pdf_;
  if(c == "rmc_ext") return get_rmc_ext_model(obs, process, selection, false).pdf_;
  if(c == "rmc_int") return get_rmc_int_model(obs, process, selection, false).pdf_;
  return nullptr;
}

void get_tail_defaults(const TString& component,
                       const TString& process,
                       double& fit_xmin,
                       int& nsmooth,
                       std::vector<double>& init_params) {
  fit_xmin = -1.;
  nsmooth = 0;
  init_params.clear();

  TString c = component;
  c.ToLower();
  const bool is_mumem = process == "mumem";

  if(c == "signal") {
    fit_xmin = (is_mumem) ? 105.5 : 93.;
    init_params = (is_mumem) ? std::vector<double>{353., -3.4} : std::vector<double>{316., -3.5};
  } else if(c == "dio") {
    fit_xmin = 102.0;
    init_params = {473., -5.};
  // } else if(c == "rmc_ext") {
  //   fit_xmin = 99.5;
  //   nsmooth = 0;
  //   init_params = {473., -4.6};
  // } else if(c == "rmc_int") {
  //   fit_xmin = 99.5;
  //   nsmooth = 0;
  //   init_params = {473., -4.6};
  }
}

int fit_component_model(TString process,
                        int selection,
                        TString tag,
                        TString component,
                        TString component_title,
                        TString pdf_type = "default",
                        int isys = -1,
                        TString tail_model = "default",
                        std::vector<int> shape_sets = {},
                        std::vector<int> control_region_sets = {},
                        TString fit_input = "primary") {
  if(use_evtana_) set_evtana_defaults();
  init_physics(tag);

  component.ToLower();
  pdf_type = resolve_pdf_type(component, pdf_type);
  tail_model = resolve_tail_model(component, tail_model);

  const TString dataset_key = component_dataset_key(component, process);
  const bool is_mumem = process == "mumem";
  const double xmin = (is_mumem) ? xmin_em_ : xmin_ep_;
  const double xmax = (is_mumem) ? xmax_em_ : xmax_ep_;

  TString figdir;
  if(component == "signal") figdir = Form("figures/signal_%s%s", process.Data(), (tag == "") ? "" : ("_" + tag).Data());
  else                       figdir = Form("figures/%s%s", component.Data(), (tag == "") ? "" : ("_" + tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir.Data(), figdir.Data()));

  TH1* h = nullptr;
  if(component == "signal") {
    if(shape_sets.empty()) h = get_signal_hist(process, selection, Form("%s_fit", component.Data()), isys);
    else                   h = load_component_hist_from_sets(process, shape_sets, Form("%s_fit_shape", component.Data()), isys, var_);
  } else {
    if(shape_sets.empty()) h = get_background_hist(dataset_key, selection, Form("%s_fit", component.Data()), isys);
    else                   h = get_background_hist_multi(dataset_key, shape_sets, Form("%s_fit_shape", component.Data()), isys);
  }
  if(!h) return 2;

  TH1* h_control = nullptr;
  if(component == "cosmic" && control_region_sets.empty()) {
    control_region_sets = default_control_region_sets("cosmic", selection);
  }
  if(!control_region_sets.empty() && component != "signal") {
    h_control = get_background_hist_multi(dataset_key, control_region_sets, Form("%s_fit_ctrl", component.Data()), isys);
  }

  TH1* h_orig = h;
  h = trim_hist(h, xmin, xmax);
  const double scale = (component == "signal")
    ? get_dataset_info(process).norm(npot_*signal_br_)
    : get_dataset_info(dataset_key).norm((dataset_key.BeginsWith("cosmic")) ? livetime_ : npot_);
  h->Scale(scale);
  if(h_control) {
    h_control = trim_hist(h_control, xmin, xmax);
    h_control->Scale(scale);
  }

  TH1* h_raw_for_ws = (TH1*) h->Clone(Form("%s_raw_for_ws", component.Data()));
  bool did_smooth = false;

  if(tail_model == "convolution" && component == "dio") {
    const int conv_set = 10;
    TH1* theory = get_background_hist("dio", conv_set, "dio_fit", -1, "MC_GenE");
    TH1* response = get_background_hist("dio", conv_set, "dio_fit", -1, "dP");
    if(theory && response) {
      TH1* convolution = convolve_with_resolution(theory, response);
      convolution = trim_hist(convolution, xmin, xmax);
      const double denom = convolution->Integral(convolution->FindBin(102.), convolution->GetNbinsX());
      const double numer = h->Integral(h->FindBin(102.), h->GetNbinsX());
      if(denom > 0.) convolution->Scale(numer / denom);
      TCanvas c;
      h_orig->Draw("E1");
      h->Draw("E1 same");
      h->SetLineColor(kRed);
      h_orig->SetAxisRange(xmin, xmax, "X");
      convolution->SetLineColor(kOrange);
      convolution->Draw("E1 same");
      c.SetLogy();
      c.SaveAs(Form("%s/%s_smoothing_%i%s.png", figdir.Data(), component.Data(), selection,
                    (isys > 0) ? Form("_sys_%i", isys) : ""));
      smooth_tail_from_reference(h, convolution, 102.0, xmax, 0.2);
      did_smooth = true;
    }
  } else if(tail_model != "none") {
    double fit_xmin(-1.);
    int nsmooth(0);
    std::vector<double> init_params;
    get_tail_defaults(component, process, fit_xmin, nsmooth, init_params);
    if(fit_xmin > 0.) {
       if(smooth_right_tail(h,
                h_orig,
                Form("%s/%s_smoothing_%i%s.png", figdir.Data(), component.Data(), selection,
                  (isys > 0) ? Form("_sys_%i", isys) : ""),
                tail_model,
                fit_xmin,
                xmax,
                nsmooth,
                0.2,
                false,
                init_params) == 0) {
         did_smooth = true;
       }
    }
  }

  TH1* h_t0_raw_for_ws = nullptr;
  TH1* h_t0_smoothed_for_ws = nullptr;
  if(include_t0_) {
    TH1* h_t0 = nullptr;
    if(shape_sets.empty()) {
      h_t0 = load_component_hist_from_dataset(dataset_key,
                                              selection,
                                              Form("%s_fit_t0", component.Data()),
                                              isys,
                                              "t0");
    } else {
      h_t0 = load_component_hist_from_sets(dataset_key,
                                           shape_sets,
                                           Form("%s_fit_t0_shape", component.Data()),
                                           isys,
                                           "t0");
    }

    if(h_t0) {
      const int trebin = (t_bin_width_ > 0.) ? int(t_bin_width_/h_t0->GetBinWidth(1) + 1.e-3) : 1;
      if(trebin > 1) h_t0->Rebin(trebin);

      const double target_raw_integral = (h_raw_for_ws) ? h_raw_for_ws->Integral() : h->Integral();
      const double t0_raw_integral = h_t0->Integral();
      if(target_raw_integral > 0. && t0_raw_integral > 0.) h_t0->Scale(target_raw_integral / t0_raw_integral);

      h_t0_raw_for_ws = (TH1*) h_t0->Clone(Form("%s_t0_raw_for_ws", component.Data()));

      TH1* h_t0_fit = (TH1*) h_t0->Clone(Form("%s_t0_smoothed_for_ws", component.Data()));
      const double txmin = h_t0_fit->GetXaxis()->GetBinLowEdge(1);
      const double txmax = h_t0_fit->GetXaxis()->GetBinUpEdge(h_t0_fit->GetNbinsX());
      TF1 t0_exp_fit(Form("%s_t0_exp_fit", component.Data()), "exp([0] + [1]*x)", txmin, txmax);
      t0_exp_fit.SetParameters(std::log(std::max(1.e-9, h_t0_fit->GetMaximum())), -1.e-3);
      const int fit_status = h_t0->Fit(&t0_exp_fit, "Q0R");
      if(fit_status == 0) {
        for(int ibin = 1; ibin <= h_t0_fit->GetNbinsX(); ++ibin) {
          const double x = h_t0_fit->GetBinCenter(ibin);
          const double val = std::max(0., t0_exp_fit.Eval(x));
          h_t0_fit->SetBinContent(ibin, val);
          h_t0_fit->SetBinError(ibin, std::sqrt(val));
        }
        const double target_smoothed_integral = (did_smooth) ? h->Integral() : target_raw_integral;
        const double t0_smoothed_integral = h_t0_fit->Integral();
        if(target_smoothed_integral > 0. && t0_smoothed_integral > 0.) {
          h_t0_fit->Scale(target_smoothed_integral / t0_smoothed_integral);
        }
        h_t0_smoothed_for_ws = h_t0_fit;
      } else {
        delete h_t0_fit;
      }
    }
  }

  RooRealVar obs(Form("obs_%i", selection), "p", (xmin+xmax)/2., xmin, xmax, "MeV/c");
  obs.setBins(h->FindBin(xmax-1.e-6) - h->FindBin(xmin+1.e-6) + 1);

  RooDataHist data_primary(Form("%s_data_hist", component.Data()), Form("%s input", component_title.Data()), obs, h);

  std::unique_ptr<RooDataHist> data_control;
  std::unique_ptr<RooDataHist> data_total;
  if(h_control) {
    data_control.reset(new RooDataHist(Form("%s_ctrl_data_hist", component.Data()),
                                       Form("%s control input", component_title.Data()), obs, h_control));
    TH1* h_total = (TH1*) h->Clone(Form("%s_total", component.Data()));
    h_total->Add(h_control);
    data_total.reset(new RooDataHist(Form("%s_total_data_hist", component.Data()),
                                     Form("%s total input", component_title.Data()), obs, h_total));
  }

  RooDataHist* fit_data = &data_primary;
  fit_input.ToLower();
  if(fit_input == "control" && data_control) fit_data = data_control.get();
  if(fit_input == "total"   && data_total)   fit_data = data_total.get();

  RooAbsPdf* model_pdf = component_analytic_pdf(component, obs, process, selection);
  const bool use_hist_pdf = (model_pdf == nullptr);
  RooAbsPdf* pdf = choose_pdf_model(pdf_type,
                                    use_hist_pdf,
                                    obs,
                                    *fit_data,
                                    model_pdf,
                                    Form("%s_%i_%s_pdf", process.Data(), selection, component.Data()),
                                    component_title + " PDF",
                                    h);

  if(component == "signal" && should_fit_pdf(pdf_type, use_hist_pdf)) {
    obs.setMax((is_mumem) ? 107. : 95.);
    if(run_component_fit(pdf, *fit_data, pdf_type, use_hist_pdf, false)) return 20;
    obs.setMax(xmax);
  } else {
    if(run_component_fit(pdf, *fit_data, pdf_type, use_hist_pdf)) return 20;
  }

  double norm_val = h->Integral();
  if(component == "cosmic" && h_control) {
    double n_sig_region  = h->Integral() / (xmax - xmin);
    double n_ctrl_region = h_control->Integral() / (xmax - xmin);
    if(n_sig_region == 0.) n_sig_region = n_ctrl_region * ((is_mumem) ? 1./754. : 1./471);
    norm_val = n_sig_region * (xmax - xmin);
  }

  RooRealVar norm(Form("%s_%i_%s_norm", process.Data(), selection, component.Data()),
                  component_title + " norm", norm_val);
  norm.setConstant(true);

  auto frame = obs.frame();
  if(component == "cosmic" && data_control) {
    data_primary.plotOn(frame, RooFit::Name("data_sr"), RooFit::MarkerColor(kRed));
    pdf->plotOn(frame, RooFit::Name("pdf"));
    data_control->plotOn(frame, RooFit::Name("data_cr"));
  } else {
    data_primary.plotOn(frame, RooFit::Name("data"));
    pdf->plotOn(frame, RooFit::Name("pdf"));
  }

  const bool use_fractional_deviation_plot = (component == "signal");
  auto c = plot_fit_frame(frame, obs, "Momentum (MeV/c)", "", "data", "pdf", npot_, livetime_, nmuons_,
                          use_fractional_deviation_plot);
  if(!c) return 10;

  const double ymax = std::max(1.e-9, h->GetMaximum());
  if(component == "rpc_ext" || component == "rpc_int" || component == "pbar") {
    frame->GetYaxis()->SetRangeUser(1.e-2*ymax, 10.*ymax);
  } else if(component == "dio") {
    frame->GetYaxis()->SetRangeUser(1.e-5, 1.e3*npot_/3.6e20);
  } else if(component == "rmc_ext") {
    frame->GetYaxis()->SetRangeUser(0., 1.2*ymax);
  }

  if(component == "signal") {
    TLatex* text = new TLatex();
    text->SetNDC();
    text->SetTextAlign(11);
    text->SetTextSize(0.042);
    text->SetTextFont(42);
    text->SetTextAlign(31);
    const int ntens = -1*std::log10(signal_br_);
    text->DrawLatex(1. - gPad->GetRightMargin() - 0.02,
                    1. - gPad->GetTopMargin() - 0.05,
                    Form("B(#mu^{-}#rightarrowe^{%s}) = %.1f x 10^{-%i}",
                         (is_mumem) ? "-" : "+", signal_br_*std::pow(10, ntens), ntens));
  }

  const TString out_base = Form("%s/%s_fit_%i%s",
                                 figdir.Data(), component.Data(), selection,
                                 (isys > 0) ? Form("_sys_%i", isys) : "");

  if(component == "dio") {
    if(save_fit_plot_pair(c, frame, out_base, 1.e-5, 1.e4*npot_/3.6e20)) return 10;
  } else if(component == "signal") {
    if(save_fit_plot_pair(c, frame, out_base, 1.e-5*ymax, 5.*ymax)) return 10;
  } else if(component == "rmc_ext") {
    if(save_fit_plot_pair(c, frame, out_base, 5.e-3*ymax, 10.*ymax)) return 10;
  } else if(component == "rmc_int") {
    if(save_fit_plot_pair(c, frame, out_base, 1.e-5, 10.*ymax)) return 10;
  } else {
    if(save_fit_plot_pair(c, frame, out_base)) return 10;
  }

  const TString suffix = (isys < 0) ? "" : Form("_sys_%i", isys);
  TH1* h_smoothed_for_ws = did_smooth ? h : nullptr;
  if(component == "signal" || component == "dio") {
    return save_fit_workspace_with_hist(process, selection, tag, component, pdf, obs, norm,
                                        h,
                                        h_raw_for_ws,
                                        h_smoothed_for_ws,
                                        h_t0_raw_for_ws,
                                        h_t0_smoothed_for_ws,
                                        suffix);
  }
  return save_fit_workspace(process, selection, tag, component, component_title, pdf, obs, norm,
                            hist_pdfs_,
                            h_raw_for_ws,
                            h_smoothed_for_ws,
                            h_t0_raw_for_ws,
                            h_t0_smoothed_for_ws,
                            suffix);
}

#endif
