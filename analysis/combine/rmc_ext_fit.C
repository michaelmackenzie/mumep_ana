// Fit the RMC external spectrum shape
#ifndef __CONVANA_ANALYSIS_RMCEXTFIT__
#define __CONVANA_ANALYSIS_RMCEXTFIT__

#include "../tools/utilities.C"
#include "../defaults.C"
#include "background_model.C"
#include "fit_workspace_utils.C"
#include "../physics.C"

int rmc_ext_fit(TString process = "mumem", int selection = 20, TString tag = "", TString pdf_type = "auto") {
  if(use_evtana_) set_evtana_defaults();
  init_physics(tag); //initialize normalization info
  const char* figdir = Form("figures/rmc_ext%s", (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir, figdir));

  //----------------------------------------------
  // Get the input data

  TH1* h = get_background_hist("rmc_ext_0n", selection, "rmc_fit");

  if(!h) {
    return 2;
  }

  const bool is_mumem = process == "mumem";
  const float xmin(is_mumem ? xmin_em_ : xmin_ep_), xmax(is_mumem ? xmax_em_ : xmax_ep_);
  const float scale = get_dataset_info("rmc_ext_0n").norm(npot_);
  cout << "Scaling inputs by " << scale << endl;
  h->Scale(scale);
  TH1* h_orig = h;
  h = trim_hist(h, xmin, xmax); // restrict to the relevant range

  // Smooth the right tail to stabilize sparse high-momentum bins.
  smooth_right_tail(h,
                    h_orig,
                    Form("%s/rmc_ext_smoothing_%i.png", figdir, selection),
                    "exp",
                    99.5,
                    xmax,
                    1,
                    0.2,
                    true,
                    {473., -4.6});

  //----------------------------------------------
  // Construct the fit objects

  // Create the observable
  RooRealVar obs(Form("obs_%i", selection), "p", (xmin+xmax)/2., xmin, xmax, "MeV/c");
  obs.setBins(h->FindBin(xmax-1.e-6) - h->FindBin(xmin+1.e-6) + 1);

  // Create the histogram data
  RooDataHist data_hist(Form("%s_%i_rmc_ext_data_hist", process.Data(), selection), "RMC (external) input", obs, h);

  // Create the PDF
  RooAbsPdf* shape_pdf = get_rmc_ext_model(obs, process, selection, false).pdf_;
  RooAbsPdf* pdf = choose_pdf_model(pdf_type,
                                    hist_pdfs_,
                                    obs,
                                    data_hist,
                                    shape_pdf,
                                    Form("%s_%i_rmc_ext_pdf", process.Data(), selection),
                                    "RMC (external) PDF",
                                    h);

  //----------------------------------------------
  // Perform the fit

  if(run_component_fit(pdf, data_hist, pdf_type, hist_pdfs_)) return 20;


  //----------------------------------------------
  // Evaluate the predicted rate

  const double n_rmc = h->Integral();
  RooRealVar norm(Form("%s_%i_rmc_ext_norm", process.Data(), selection), "RMC (external) norm", n_rmc);
  norm.setConstant(true);

  //----------------------------------------------
  // Plot the results

  auto frame = obs.frame();
  data_hist.plotOn(frame, RooFit::Name("data"));
  pdf->plotOn(frame, RooFit::Name("pdf"));
  shape_pdf->plotOn(frame, RooFit::Name("shape_pdf"), RooFit::LineColor(kGreen));
  auto c = plot_fit_frame(frame, obs, "momentum (MeV/c)", "", "data", "pdf", npot_, livetime_);
  if(!c) return 10;
  frame->GetYaxis()->SetRangeUser(0., 1.2*h->GetMaximum());
  if(save_fit_plot_pair(c, frame, Form("%s/rmc_ext_fit_%i", figdir, selection),
                        5.e-3*h->GetMaximum(), 10.*h->GetMaximum())) return 10;

  //----------------------------------------------
  // Save the results

  return save_fit_workspace(process, selection, tag, "rmc_ext", "RMC (external)", pdf, obs, norm, hist_pdfs_);
}

#endif
