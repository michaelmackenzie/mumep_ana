// Fit the RMC internal spectrum shape
#ifndef __CONVANA_ANALYSIS_RMCINTFIT__
#define __CONVANA_ANALYSIS_RMCINTFIT__

#include "../tools/utilities.C"
#include "../defaults.C"
#include "background_model.C"
#include "fit_workspace_utils.C"
#include "../physics.C"

int rmc_int_fit(TString process = "mumem", int selection = 20, TString tag = "", TString pdf_type = "auto") {
  if(use_evtana_) set_evtana_defaults();
  init_physics(tag); //initialize normalization info
  const char* figdir = Form("figures/rmc_int%s", (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir, figdir));

  //----------------------------------------------
  // Get the input data

  TH1* h = get_background_hist("rmc_int", selection, "rmc_fit");

  if(!h) {
    return 2;
  }

  const bool is_mumem = process == "mumem";
  const float xmin(is_mumem ? xmin_em_ : xmin_ep_), xmax(is_mumem ? xmax_em_ : xmax_ep_);
  const float scale = get_dataset_info("rmc_int").norm(npot_);
  cout << "Scaling inputs by " << scale << endl;
  h->Scale(scale);
  TH1* h_orig = h;
  h = trim_hist(h, xmin, xmax); // restrict to the relevant range

  smooth_right_tail(h,
                    h_orig,
                    Form("%s/rmc_int_smoothing_%i.png", figdir, selection),
                    "power",
                    89.,
                    xmax,
                    0,
                    0.2,
                    true,
                    {89., 10., -1.8});

  //----------------------------------------------
  // Construct the fit objects

  // Create the observable
  RooRealVar obs(Form("obs_%i", selection), "p", (xmin+xmax)/2., xmin, xmax, "MeV/c");
  obs.setBins(h->FindBin(xmax-1.e-6) - h->FindBin(xmin+1.e-6) + 1);

  // Create the histogram data
  RooDataHist data_hist(Form("%s_%i_rmc_int_data_hist", process.Data(), selection), "RMC (internal) input", obs, h);

  // Create the PDF
  RooAbsPdf* model_pdf = nullptr;
  RooAbsPdf* pdf = choose_pdf_model(pdf_type,
                                    hist_pdfs_,
                                    obs,
                                    data_hist,
                                    model_pdf,
                                    Form("%s_%i_rmc_int_pdf", process.Data(), selection),
                                    "RMC (internal) PDF",
                                    h);

  //----------------------------------------------
  // Perform the fit

  if(run_component_fit(pdf, data_hist, pdf_type, hist_pdfs_)) return 20;


  //----------------------------------------------
  // Evaluate the predicted rate

  const double n_rmc = h->Integral();
  RooRealVar norm(Form("%s_%i_rmc_int_norm", process.Data(), selection), "RMC (internal) norm", n_rmc);
  norm.setConstant(true);

  //----------------------------------------------
  // Plot the results

  auto frame = obs.frame();
  data_hist.plotOn(frame, RooFit::Name("data"));
  pdf->plotOn(frame, RooFit::Name("pdf"));
  auto c = plot_fit_frame(frame, obs, "momentum (MeV/c)", "", "data", "pdf", npot_, livetime_);
  if(!c) return 10;
  if(save_fit_plot_pair(c, frame, Form("%s/rmc_int_fit_%i", figdir, selection),
                        1.e-5, 10.*h->GetMaximum())) return 10;

  //----------------------------------------------
  // Save the results

  return save_fit_workspace(process, selection, tag, "rmc_int", "RMC (internal)", pdf, obs, norm, hist_pdfs_);
}

#endif
