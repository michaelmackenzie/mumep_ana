// Fit the DIO spectrum shape
#ifndef __CONVANA_ANALYSIS_DIOFIT__
#define __CONVANA_ANALYSIS_DIOFIT__

#include "../tools/utilities.C"
#include "../tools/functions.C"
#include "../defaults.C"
#include "background_model.C"
#include "fit_workspace_utils.C"
#include "../physics.C"

int dio_fit(TString process = "mumem", int selection = 20, TString tag = "", const int isys = -1,
            TString pdf_type = "auto", TString tail_model = "convolution") {
  if(use_evtana_) set_evtana_defaults();
  init_physics(tag); //initialize normalization info
  const char* figdir = Form("figures/dio%s", (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir, figdir));

  //----------------------------------------------
  // Get the input data

  TH1* h = get_background_hist("dio", selection, "dio_fit", isys);
  if(!h) {
    return 2;
  }

  const bool is_mumem = process == "mumem";
  const float xmin(is_mumem ? xmin_em_ : xmin_ep_), xmax(is_mumem ? xmax_em_ : xmax_ep_);
  TH1* h_orig = h;
  h = trim_hist(h, xmin, xmax); // restrict to the relevant range
  const float scale = get_dataset_info("dio").norm(npot_);
  h->Scale(scale);
  h_orig->Scale(scale);

  // Get the theory (with efficiencies) and resolution function
  const int conv_set = 10; // FIXME: Using a set with wider momentum window
  TH1* theory   = get_background_hist("dio", conv_set, "dio_fit", -1, "MC_GenE");
  TH1* response = get_background_hist("dio", conv_set, "dio_fit", -1, "dP");
  if(!theory || !response) {
    cout << __func__ << ": Theory and/or response not found!\n";
    return 1;
  }

  // Fit the resolution
  TF1* dscb = double_crystal_ball_tf1();
  dscb->SetParameters(1., 0.2, 0., 1., 1., 3., 3.);
  response->Fit(dscb, "R");

  // Convolve theory and resolution
  TH1* convolution = convolve_with_resolution(theory, response);
  convolution = trim_hist(convolution, xmin, xmax); // restrict to the relevant range
  convolution->Scale(h->Integral(h->FindBin(102.), h->GetNbinsX()) /
                     convolution->Integral(convolution->FindBin(102.), convolution->GetNbinsX()));

  // Smooth the right tail with either convolution-derived shape or parametric tail fit.
  const double fit_xmin = 102.0;
  const TString smoothing_plot = Form("%s/dio_smoothing_%i%s.png", figdir, selection,
                                      (isys > 0) ? Form("_sys_%i", isys) : "");
  TString tail_mode = tail_model;
  tail_mode.ToLower();
  if(tail_mode == "convolution") {
    TCanvas c;
    h_orig->Draw("E1");
    h->Draw("E1 same");
    h->SetLineColor(kRed);
    h_orig->SetAxisRange(xmin, xmax, "X");
    convolution->SetLineColor(kOrange);
    convolution->Draw("E1 same");
    c.SetLogy();
    c.SaveAs(smoothing_plot.Data());
    smooth_tail_from_reference(h, convolution, fit_xmin, xmax, 0.2);
  } else if(tail_mode != "none") {
    smooth_right_tail(h,
                      h_orig,
                      smoothing_plot,
                      tail_mode,
                      fit_xmin,
                      xmax,
                      0,
                      0.2,
                      false,
                      {473., -5.});
  }

  //----------------------------------------------
  // Construct the fit objects

  // Create the observable
  RooRealVar obs(Form("obs_%i", selection), "p", (xmin+xmax)/2., xmin, xmax, "MeV/c");
  obs.setBins(h->FindBin(xmax-1.e-6) - h->FindBin(xmin+1.e-6) + 1);

  // Create the histogram data
  RooDataHist data_hist(Form("%s_%i_dio_data_hist", process.Data(), selection), "DIO histogram input", obs, h);

  // Create the PDF
  RooAbsPdf* model_pdf = get_dio_model(obs, process, selection, true).pdf_;
  RooAbsPdf* pdf = choose_pdf_model(pdf_type,
                                    hist_pdfs_,
                                    obs,
                                    data_hist,
                                    model_pdf,
                                    Form("%s_%i_dio_pdf", process.Data(), selection),
                                    "DIO PDF",
                                    h);


  //----------------------------------------------
  // Perform the fit

  if(run_component_fit(pdf, data_hist, pdf_type, hist_pdfs_)) return 20;

  // Plot the underlying PDFs if convolving
  if(pdf->InheritsFrom("RooFFTConvPdf") || pdf->InheritsFrom("RooNumConvPdf")) {
    for(auto comp : *pdf->getComponents()) {
      auto frame = obs.frame();
      ((RooAbsPdf*) comp)->plotOn(frame);
      TCanvas* c = new TCanvas();
      frame->Draw();
      c->SetLogy();
      c->SaveAs(Form("%s/comp_%s_%i%s.png", figdir, comp->GetName(), selection, (isys > 0) ? Form("_sys_%i", isys) : ""));
    }
  }

  //----------------------------------------------
  // Evaluate the predicted rate

  const double n_dio = h->Integral();
  RooRealVar norm(Form("%s_%i_dio_norm", process.Data(), selection), "DIO norm", n_dio);
  norm.setConstant(true);

  //----------------------------------------------
  // Plot the results

  auto frame = obs.frame();
  data_hist.plotOn(frame, RooFit::Name("data"));
  pdf->plotOn(frame, RooFit::Name("pdf"));
  auto c = plot_fit_frame(frame, obs, "momentum (MeV/c)", "", "data", "pdf", npot_, livetime_);
  if(!c) return 10;
  frame->GetYaxis()->SetRangeUser(1.e-5, 1.e3*npot_/3.6e20);
  if(save_fit_plot_pair(c, frame, Form("%s/dio_fit_%i%s", figdir, selection,
                                       (isys > 0) ? Form("_sys_%i", isys) : ""),
                        1.e-5, 1.e4*npot_/3.6e20)) return 10;

  //----------------------------------------------
  // Save the results

  return save_fit_workspace_with_hist(process, selection, tag, "dio", pdf, norm,
                                      h, (isys < 0) ? "" : Form("_sys_%i", isys));
}

#endif
