// Fit the antiproton spectrum shape
#ifndef __CONVANA_ANALYSIS_PBARFIT__
#define __CONVANA_ANALYSIS_PBARFIT__

#include "../tools/utilities.C"
#include "../defaults.C"
#include "background_model.C"
#include "fit_workspace_utils.C"
#include "../physics.C"

int pbar_fit(TString process = "mumem", int selection = 20, TString tag = "", TString pdf_type = "auto") {
  if(use_evtana_) set_evtana_defaults();
  init_physics(tag); //initialize normalization info
  const char* figdir = Form("figures/pbar%s", (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir, figdir));

  //----------------------------------------------
  // Get the input data

  TH1* h = get_background_hist("pbar", selection, "pbar");

  if(!h) {
    return 2;
  }

  const bool is_mumem = process == "mumem";
  const float xmin(is_mumem ? xmin_em_ : xmin_ep_), xmax(is_mumem ? xmax_em_ : xmax_ep_);
  h = trim_hist(h, xmin, xmax); // restrict to the relevant range
  const float scale = get_dataset_info("pbar").norm(npot_);
  cout << "Scaling inputs by " << scale << endl;
  h->Scale(scale);

  //----------------------------------------------
  // Construct the fit objects

  // Create the observable
  RooRealVar obs(Form("obs_%i", selection), "p", (xmin+xmax)/2., xmin, xmax, "MeV/c");
  obs.setBins(h->FindBin(xmax-1.e-6) - h->FindBin(xmin+1.e-6) + 1);

  // Create the histogram data
  RooDataHist data_hist("pbar_data_hist", "Antiproton input", obs, h);

  // Create the PDF
  RooAbsPdf* model_pdf = get_pbar_model(obs, process, selection, false).pdf_;
  RooAbsPdf* pdf = choose_pdf_model(pdf_type,
                                    hist_pdfs_,
                                    obs,
                                    data_hist,
                                    model_pdf,
                                    Form("%s_%i_pbar_pdf", process.Data(), selection),
                                    "Antiproton background",
                                    h);

  //----------------------------------------------
  // Perform the fit

  if(run_component_fit(pdf, data_hist, pdf_type, hist_pdfs_)) return 20;


  //----------------------------------------------
  // Evaluate the predicted rate

  const double n_pbar = h->Integral();
  RooRealVar norm(Form("%s_%i_pbar_norm", process.Data(), selection), "Antiproton norm", n_pbar);
  norm.setConstant(true);

  //----------------------------------------------
  // Plot the results

  auto frame = obs.frame();
  data_hist.plotOn(frame, RooFit::Name("data"));
  pdf->plotOn(frame, RooFit::Name("pdf"));
  auto c = plot_fit_frame(frame, obs, "momentum (MeV/c)", "", "data", "pdf", npot_, livetime_);
  if(!c) return 10;
  frame->GetYaxis()->SetRangeUser(1.e-2*h->GetMaximum(), 10.*h->GetMaximum());
  if(save_fit_plot_pair(c, frame, Form("%s/pbar_fit_%i", figdir, selection))) return 10;

  //----------------------------------------------
  // Save the results

  return save_fit_workspace(process, selection, tag, "pbar", "Antiproton", pdf, obs, norm, hist_pdfs_);
}

#endif
