// Fit the RPC internal spectrum shape
#ifndef __CONVANA_ANALYSIS_RPCINTFIT__
#define __CONVANA_ANALYSIS_RPCINTFIT__

#include "../tools/utilities.C"
#include "../defaults.C"
#include "background_model.C"
#include "fit_workspace_utils.C"
#include "../physics.C"

int rpc_int_fit(TString process = "mumem", int selection = 20, TString tag = "") {
  if(use_evtana_) set_evtana_defaults();
  init_physics(tag); //initialize normalization info
  const char* figdir = Form("figures/rpc_int%s", (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir, figdir));

  //----------------------------------------------
  // Get the input data

  TH1* h = get_background_hist("rpc_int", selection, "rpc_fit");

  if(!h) {
    return 2;
  }

  const bool is_mumem = process == "mumem";
  const float xmin(is_mumem ? xmin_em_ : xmin_ep_), xmax(is_mumem ? xmax_em_ : xmax_ep_);
  h = trim_hist(h, xmin, xmax); // restrict to the relevant range
  const float scale = get_dataset_info("rpc_int").norm(npot_);
  cout << "Scaling inputs by " << scale << endl;
  h->Scale(scale);

  //----------------------------------------------
  // Construct the fit objects

  // Create the observable
  RooRealVar obs(Form("obs_%i", selection), "p", (xmin+xmax)/2., xmin, xmax, "MeV/c");
  obs.setBins(h->FindBin(xmax-1.e-6) - h->FindBin(xmin+1.e-6) + 1);

  // Create the histogram data
  RooDataHist data_hist("rpc_int_data_hist", "RPC (internal) input", obs, h);

  // Create the PDF
  RooAbsPdf* pdf = get_rpc_int_model(obs, process, selection, false).pdf_;

  //----------------------------------------------
  // Perform the fit

  const TString name = Form("%s_%i_rpc_int", process.Data(), selection);
  enforce_uniform_if_sparse(h, obs, pdf, Form("%s_pdf", name.Data()), "RPC (internal) background");
  pdf->fitTo(data_hist, RooFit::SumW2Error(true));


  //----------------------------------------------
  // Evaluate the predicted rate

  const double n_rpc = h->Integral();
  RooRealVar norm(Form("%s_%i_rpc_int_norm", process.Data(), selection), "RPC (internal) norm", n_rpc);
  norm.setConstant(true);

  //----------------------------------------------
  // Plot the results

  auto frame = obs.frame();
  data_hist.plotOn(frame, RooFit::Name("data"));
  pdf->plotOn(frame, RooFit::Name("pdf"));
  auto c = plot_fit_frame(frame, obs, "momentum (MeV/c)", "", "data", "pdf", npot_, livetime_);
  if(!c) return 10;
  frame->GetYaxis()->SetRangeUser(1.e-2*h->GetMaximum(), 10.*h->GetMaximum());
  if(save_fit_plot_pair(c, frame, Form("%s/rpc_int_fit_%i", figdir, selection))) return 10;

  //----------------------------------------------
  // Save the results

  return save_fit_workspace(process, selection, tag, "rpc_int", "RPC (internal)", pdf, obs, norm, hist_pdfs_);
}

#endif
