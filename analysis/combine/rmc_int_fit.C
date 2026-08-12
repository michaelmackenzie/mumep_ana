// Fit the RMC internal spectrum shape
#ifndef __CONVANA_ANALYSIS_RMCINTFIT__
#define __CONVANA_ANALYSIS_RMCINTFIT__

#include "../tools/utilities.C"
#include "../defaults.C"
#include "background_model.C"
#include "fit_workspace_utils.C"
#include "../physics.C"

int rmc_int_fit(TString process = "mumem", int selection = 20, TString tag = "") {
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

  // Attempt to smooth the histogram a bit
  const int nsmooth(0);
  h->Smooth(nsmooth);
  {
    TCanvas c;
    h_orig->Draw("E1");
    h->Draw("E1 same");
    h->SetLineColor(kRed);
    h_orig->SetAxisRange(xmin, xmax, "X");

    // Attempt to fit the histogram tail to further smooth it
    const double fit_xmin(89.);
    TF1* tail_func = new TF1("tail_func", "pow((x-[0])/[1],[2])", fit_xmin, xmax);
    tail_func->SetParameters(89, 10, -1.8);
    // TF1* tail_func = new TF1("tail_func", "exp([0] + [1]*x)", fit_xmin, xmax);
    // tail_func->SetParameters(473., -4.6);
    // tail_func->FixParameter(1, -4.6);
    h->Fit(tail_func, "wR");
    tail_func->Draw("same");
    c.SetLogy();
    c.SaveAs(Form("%s/rmc_int_smoothing_%i.png", figdir, selection));

    for(int ibin = h->FindBin(fit_xmin+1.e-6); ibin <= h->FindBin(xmax-1.e-6); ++ibin) {
      const float x = h->GetBinCenter(ibin);
      const float val = tail_func->Eval(x);
      h->SetBinContent(ibin, val);
      h->SetBinError(ibin, 0.2*val); //default to 20% uncertainty on the fit result
    }
  }

  //----------------------------------------------
  // Construct the fit objects

  // Create the observable
  RooRealVar obs(Form("obs_%i", selection), "p", (xmin+xmax)/2., xmin, xmax, "MeV/c");
  obs.setBins(h->FindBin(xmax-1.e-6) - h->FindBin(xmin+1.e-6) + 1);

  // Create the histogram data
  RooDataHist data_hist(Form("%s_%i_rmc_int_data_hist", process.Data(), selection), "RMC (internal) input", obs, h);

  // Create the PDF
  RooAbsPdf* pdf = new RooHistPdf(Form("%s_%i_rmc_int_pdf", process.Data(), selection), "RMC (internal) PDF", obs, data_hist);
    //get_rmc_int_model(obs, process, selection, false).pdf_;

  //----------------------------------------------
  // Perform the fit

  // pdf->fitTo(data_hist, RooFit::SumW2Error(true));


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
