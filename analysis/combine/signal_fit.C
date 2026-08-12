// Fit the DIO spectrum shape
#ifndef __CONVANA_ANALYSIS_SIGNALFIT__
#define __CONVANA_ANALYSIS_SIGNALFIT__

#include "../tools/utilities.C"
#include "../defaults.C"
#include "fit_workspace_utils.C"
#include "signal_model.C"
#include "../physics.C"

//---------------------------------------------------------------------------------------------------------------------------
int signal_fit(TString process = "mumem", int selection = 20, TString tag = "", const int isys = -1) {
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

  // Attempt to smooth the histogram a bit
  {
    TCanvas c;
    h_orig->Draw("E1");
    h_orig->SetAxisRange(xmin, xmax, "X");

    // Attempt to fit the histogram tail to further smooth it
    const double fit_xmin((is_mumem) ? 105.5 : 93.);
    TF1* tail_func = new TF1("tail_func", "exp([0] + [1]*x)", fit_xmin, xmax);
    if(is_mumem) tail_func->SetParameters(353., -3.4);
    else         tail_func->SetParameters(316., -3.5);
    h->Fit(tail_func, "R");
    tail_func->Draw("same");
    c.SetLogy();
    c.SaveAs(Form("%s/signal_smoothing_%i%s.png", figdir, selection, (isys < 0) ? "" : Form("_sys_%i", isys)));

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

  // Create the histogram data
  obs.setBins(h->GetNbinsX());
  RooDataHist data_hist("signal_data_hist", "Signal input", obs, h);

  // Create the PDF
  RooAbsPdf* pdf = (hist_pdfs_) ?
    new RooHistPdf(Form("%s_%i_signal_pdf", process.Data(), selection), "Signal PDF", obs, data_hist) :
    get_signal_model(obs, process, selection, false).pdf_;
  if(!pdf) {
    cout << "No PDF returned!\n";
    return 10;
  }

  //----------------------------------------------
  // Perform the fit

  if(!hist_pdfs_) {
    obs.setMax((is_mumem) ? 107. : 95.);
    pdf->fitTo(data_hist);
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
