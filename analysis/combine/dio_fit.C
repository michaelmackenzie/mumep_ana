// Fit the DIO spectrum shape
#ifndef __CONVANA_ANALYSIS_DIOFIT__
#define __CONVANA_ANALYSIS_DIOFIT__

#include "../tools/utilities.C"
#include "../tools/functions.C"
#include "../defaults.C"
#include "background_model.C"
#include "fit_workspace_utils.C"
#include "../physics.C"

int dio_fit(TString process = "mumem", int selection = 20, TString tag = "", const int isys = -1) {
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

  // Attempt to smooth the histogram a bit
  const int nsmooth(0);
  h->Smooth(nsmooth);
  {
    TCanvas c;
    h_orig->Draw("E1");
    h->Draw("E1 same");
    h->SetLineColor(kRed);
    h_orig->SetAxisRange(xmin, xmax, "X");
    convolution->SetLineColor(kOrange);
    convolution->Draw("E1 same");

    // Attempt to fit the histogram tail to further smooth it
    const double tail_integral = h->Integral(h->FindBin(104.), h->GetNbinsX());
    const double fit_xmin((tail_integral > 0.) ? 102.0 : 102.), fit_xmax(106.);

    // First increase bin errors if it's far below its neighbors
    for(int ibin = h->FindBin(fit_xmin); ibin <= h->FindBin(fit_xmax); ++ibin) {
      if(ibin <= 1 || ibin >= h->GetNbinsX()) continue; // no neighbors
      const double bine       = h->GetBinError(ibin  );
      const double bine_left  = h->GetBinError(ibin-1);
      const double bine_right = h->GetBinError(ibin+1);
      if(bine_left > 0. && bine_right > 0.) {
        if(bine < bine_left && bine < bine_right) h->SetBinError(ibin, (bine_left + bine_right)/2.); //average the bin errors
        else if(bine < 1.e-2*bine_left) h->SetBinError(ibin, (bine_left/5.)); // suspiciously small error
      } else if(bine_left > 0. && bine_left / bine > 100.) h->SetBinError(ibin, bine_left/2.); // if it's very low error, increase it
    }

    TF1* tail_func = new TF1("tail_func", "exp([0] + [1]*x)", fit_xmin, fit_xmax);
    tail_func->SetParameters(473., -5.);
    // tail_func->FixParameter(1, -4.6);
    h->Fit(tail_func, "R");
    tail_func->Draw("same");

    c.SetLogy();
    gStyle->SetOptStat(1001111);
    // h_orig->GetYaxis()->SetRangeUser(1.e-10, 1.e4*npot_/3.6e20);
    c.SaveAs(Form("%s/dio_smoothing_%i%s.png", figdir, selection, (isys > 0) ? Form("_sys_%i", isys) : ""));

    // // Use the fit to smooth
    // for(int ibin = h->FindBin(fit_xmin+1.e-6); ibin <= h->FindBin(xmax-1.e-6); ++ibin) {
    //   const float x = h->GetBinCenter(ibin);
    //   const float val = tail_func->Eval(x);
    //   h->SetBinContent(ibin, val);
    //   h->SetBinError(ibin, 0.2*val); //default to 20% uncertainty on the fit result
    // }

    // Use the convolution to smooth
    for(int ibin = h->FindBin(fit_xmin+1.e-6); ibin <= h->FindBin(xmax-1.e-6); ++ibin) {
      const float x = h->GetBinCenter(ibin);
      const float val = convolution->Interpolate(x) * h->GetBinWidth(ibin) / convolution->GetBinWidth(ibin);
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
  RooDataHist data_hist(Form("%s_%i_dio_data_hist", process.Data(), selection), "DIO histogram input", obs, h);

  // Create the PDF
  RooAbsPdf* pdf = (hist_pdfs_) ?
    new RooHistPdf(Form("%s_%i_dio_pdf", process.Data(), selection), "DIO PDF", obs, data_hist) :
    get_dio_model(obs, process, selection, true).pdf_;


  //----------------------------------------------
  // Perform the fit

  if(!hist_pdfs_) pdf->fitTo(data_hist, RooFit::SumW2Error(true));

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
