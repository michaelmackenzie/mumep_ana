// Fit the DIO spectrum shape
#ifndef __CONVANA_ANALYSIS_COSMICFIT__
#define __CONVANA_ANALYSIS_COSMICFIT__

#include "../tools/utilities.C"
#include "../defaults.C"
#include "background_model.C"
#include "fit_workspace_utils.C"
#include "../physics.C"

bool use_control_region_ = false; // whether or not to take the distribution from the CRV tagged region
bool fit_total_          = true ; // fit the combined veto + signal regions to get more stable results

int cosmic_fit(TString process = "mumem", int selection = 20, TString tag = "", const int isys = -1, TString pdf_type = "auto") {
  if(use_evtana_) set_evtana_defaults();
  init_physics(tag); //initialize normalization info
  const char* figdir = Form("figures/cosmic%s", (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir, figdir));

  //----------------------------------------------
  // Get the input data

  TH1* h_sig_region  = get_background_hist("cosmic", selection     , "cosmic_fit", isys);
  TH1* h_veto_region = get_background_hist("cosmic", selection+1000, "cosmic_fit", isys);
  if(!h_sig_region || !h_veto_region) {
    return 2;
  }

  const bool is_mumem = process == "mumem";
  const float xmin(is_mumem ? xmin_em_ : xmin_ep_), xmax(is_mumem ? xmax_em_ : xmax_ep_);
  h_sig_region  = trim_hist(h_sig_region , xmin, xmax); // restrict to the relevant range
  h_veto_region = trim_hist(h_veto_region, xmin, xmax);
  const float scale = get_dataset_info("cosmic").norm(livetime_);
  cout << "Scaling input histograms by " << scale << " to account for the livetime\n";
  h_sig_region ->Scale(scale);
  h_veto_region->Scale(scale);
  TH1* h_total = (TH1*) h_sig_region->Clone("h_total");
  h_total->Add(h_veto_region);

  //----------------------------------------------
  // Construct the fit objects

  // Create the observable
  RooRealVar obs(Form("obs_%i", selection), "p", (xmin+xmax)/2., xmin, xmax, "MeV/c");
  obs.setBins(h_sig_region->FindBin(xmax-1.e-6) - h_sig_region->FindBin(xmin+1.e-6) + 1);

  // Create the histogram data
  RooDataHist data_sig_hist ("cosmic_sig_data_hist" , "Cosmic input (signal region)", obs, h_sig_region );
  RooDataHist data_veto_hist("cosmic_veto_data_hist", "Cosmic input (veto region)"  , obs, h_veto_region);
  RooDataHist data_tot_hist ("cosmic_tot_data_hist" , "Cosmic input (total)"        , obs, h_total);

  // Create the PDF
  RooAbsPdf* model_pdf = get_cosmic_model(obs, process, selection, false).pdf_;
  RooAbsPdf* pdf = choose_pdf_model(pdf_type,
                                    hist_pdfs_,
                                    obs,
                                    data_sig_hist,
                                    model_pdf,
                                    Form("%s_%i_cosmic_pdf", process.Data(), selection),
                                    "Cosmic background",
                                    h_sig_region);

  //----------------------------------------------
  // Perform the fit without the veto for the shape

  RooDataHist* fit_data = &data_sig_hist;
  if(fit_total_)            fit_data = &data_tot_hist;
  else if(use_control_region_) fit_data = &data_veto_hist;
  if(run_component_fit(pdf, *fit_data, pdf_type, hist_pdfs_)) return 20;

  //----------------------------------------------
  // Evaluate the predicted rate

  double n_sig_region  = h_sig_region ->Integral() / (xmax - xmin);
  double n_veto_region = h_veto_region->Integral() / (xmax - xmin);
  //FIXME: Set better rate expectations when no events are found
  if(n_sig_region == 0.) {
    n_sig_region = n_veto_region * ((is_mumem) ? 1./754. : 1./471);
  }
  double ratio     = (n_sig_region > 0.) ? n_veto_region / n_sig_region : 0.;
  printf("--- N(signal) = %7.3f / MeV/c; N(veto) = %7.3f / MeV/c; Ratio = %5.2f\n",
         n_sig_region, n_veto_region, ratio);

  RooRealVar norm(Form("%s_%i_cosmic_norm", process.Data(), selection), "Cosmic norm", n_sig_region * (xmax - xmin));
  norm.setConstant(true);

  //----------------------------------------------
  // Plot the results

  auto frame = obs.frame();
  if(use_control_region_) {
    data_veto_hist.plotOn(frame, RooFit::Name("data_cr"));
    pdf->plotOn(frame, RooFit::Name("pdf"));
    data_sig_hist.plotOn(frame, RooFit::Name("data_sr"), RooFit::MarkerColor(kRed));
  } else {
    data_sig_hist.plotOn(frame, RooFit::Name("data_sr"), RooFit::MarkerColor(kRed));
    pdf->plotOn(frame, RooFit::Name("pdf"));
    data_veto_hist.plotOn(frame, RooFit::Name("data_cr"));
  }
  auto c = plot_fit_frame(frame, obs, "momentum (MeV/c)", "", "data", "pdf", npot_, livetime_);
  if(!c) return 10;
  // frame->GetYaxis()->SetRangeUser(1.e-5, 1.e2);
  if(save_fit_plot_pair(c, frame, Form("%s/cosmic_fit_%i%s", figdir, selection,
                                       (isys > 0) ? Form("_sys_%i", isys) : ""))) return 10;

  //----------------------------------------------
  // Save the results

  cout << "\n============================================================\n"
       << "N(signal region cosmics) in 103.4 < p < 104.8 MeV/c: "
       << h_sig_region->Integral(h_sig_region->FindBin(103.401), h_sig_region->FindBin(104.799))
       << "\n============================================================\n";

  return save_fit_workspace(process, selection, tag, "cosmic", "Cosmic ray", pdf, obs, norm,
                            hist_pdfs_, (isys > 0) ? Form("_sys_%i", isys) : "");
}

#endif
