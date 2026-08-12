// Fit the DIO spectrum shape
#ifndef __CONVANA_ANALYSIS_SIGNALFIT__
#define __CONVANA_ANALYSIS_SIGNALFIT__

#include "../tools/utilities.C"
#include "../defaults.C"
#include "fit_workspace_utils.C"
#include "signal_model.C"
#include "../physics.C"

//---------------------------------------------------------------------------------------------------------------------------
int signal_fit(TString process = "mumem", int selection = 20, TString tag = "", const int isys = -1,
               TString pdf_type = "auto", TString tail_model = "exp") {
  if(use_evtana_) set_evtana_defaults();
  init_physics(tag); //initialize normalization info
  const bool is_mumem = process == "mumem";
  const char* figdir = Form("figures/signal_%s%s", process.Data(), (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir, figdir));

  //----------------------------------------------
  // Get the input data

  TH1* h = get_signal_hist(process, selection, "signal", isys);
  if(!h) {
    cout << __func__ << ": Input histogram for selection " << selection << " not found\n";
    return 2;
  }
  const double scale = get_dataset_info(process).norm(npot_*signal_br_);
  printf("Scaling input histogram by %.3e (npot = %.3e, signal_br = %.3e)\n", scale, npot_, signal_br_);
  h->Scale(scale);

  const float xmin(is_mumem ? xmin_em_ : xmin_ep_), xmax(is_mumem ? xmax_em_ : xmax_ep_);
  TH1* h_orig = h;
  h = trim_hist(h, xmin, xmax);
  cout << "Nominal integral: " << h->Integral() << endl;

  // Smooth the right side of the signal spectrum to suppress sparse-bin spikes.
  const double fit_xmin = (is_mumem) ? 105.5 : 93.;
  std::vector<double> tail_params;
  if(tail_model == "exp") tail_params = (is_mumem) ? std::vector<double>{353., -3.4} : std::vector<double>{316., -3.5};
  smooth_right_tail(h,
                    h_orig,
                    Form("%s/signal_smoothing_%i%s.png", figdir, selection, (isys < 0) ? "" : Form("_sys_%i", isys)),
                    tail_model,
                    fit_xmin,
                    xmax,
                    0,
                    0.2,
                    false,
                    tail_params);

  //----------------------------------------------
  // Construct the fit objects

  // Create the observable
  RooRealVar obs(Form("obs_%i", selection), "p", (xmin+xmax)/2., xmin, xmax, "MeV/c");

  // Create the histogram data
  obs.setBins(h->GetNbinsX());
  RooDataHist data_hist("signal_data_hist", "Signal input", obs, h);

  // Create the PDF
  RooAbsPdf* model_pdf = get_signal_model(obs, process, selection, false).pdf_;
  RooAbsPdf* pdf = choose_pdf_model(pdf_type,
                                    hist_pdfs_,
                                    obs,
                                    data_hist,
                                    model_pdf,
                                    Form("%s_%i_signal_pdf", process.Data(), selection),
                                    "Signal PDF",
                                    h);
  if(!pdf) {
    cout << "No PDF returned!\n";
    return 10;
  }

  //----------------------------------------------
  // Perform the fit

  if(should_fit_pdf(pdf_type, hist_pdfs_)) {
    obs.setMax((is_mumem) ? 107. : 95.);
    if(run_component_fit(pdf, data_hist, pdf_type, hist_pdfs_, false)) return 20;
    obs.setMax(xmax);
  }


  //----------------------------------------------
  // Evaluate the predicted rate

  const double n_signal = h->Integral();
  RooRealVar norm(Form("%s_%i_signal_norm", process.Data(), selection), "Signal norm", n_signal);
  norm.setConstant(true);

  //----------------------------------------------
  // Plot the results

  auto frame = obs.frame();
  data_hist.plotOn(frame, RooFit::Name("data"));
  pdf->plotOn(frame, RooFit::Name("pdf"));
  auto c = plot_fit_frame(frame, obs, "momentum (MeV/c)", "", "data", "pdf", npot_, livetime_);
  if(!c) return 10;

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextAlign(11);
  text->SetTextSize(0.042);
  text->SetTextFont(42);
  text->SetTextAlign(31);
  const int ntens = -1*std::log10(signal_br_);
  text->DrawLatex(1. - gPad->GetRightMargin() - 0.02, 1. - gPad->GetTopMargin() - 0.05, Form("B(#mu^{-}#rightarrowe^{%s}) = %.1f x 10^{-%i}", (is_mumem) ? "-" : "+",
                                                                                             signal_br_*std::pow(10, ntens), ntens));
  const float ymax = h->GetMaximum();
  if(save_fit_plot_pair(c, frame, Form("%s/signal_fit_%i%s", figdir, selection,
                                       (isys < 0) ? "" : Form("_sys_%i", isys)),
                        1.e-5*ymax, 5.*ymax)) return 10;

  //----------------------------------------------
  // Save the results

  return save_fit_workspace_with_hist(process, selection, tag, "signal", pdf, norm,
                                      h, (isys < 0) ? "" : Form("_sys_%i", isys));
}

#endif
