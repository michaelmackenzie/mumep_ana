#ifndef __CONVANA_ANALYSIS_SIGNALMODEL__
#define __CONVANA_ANALYSIS_SIGNALMODEL__

#include "../tools/types.C"
#include "../defaults.C"
#include "../physics.C"
#include "../datasets.C"
#include "../tools/utilities.C"
#include "model_io_utils.C"
#include "combine/HiggsAnalysis/CombinedLimit/src/RooLandauCB.cc"
#include "RooTFnBinding.h"

//---------------------------------------------------------------------------------------------------------------------------
TH1* get_signal_hist(const TString process, const int selection, const TString name = "signal", const int isys = -1) {
  return load_component_hist_from_dataset(process, selection, name, isys, var_);
}

//---------------------------------------------------------------------------------------------------------------------------
pdf_info get_signal_model(RooRealVar& obs, const TString process, const int selection, const bool freeze = true, const bool use_hist = false) {
  const char* name = Form("%s_%i_signal", process.Data(), selection);
  const bool is_mumem = process == "mumem";

  double eff = 0.3; //signal efficiency

  TH1* h = get_signal_hist(process, selection, "signal_model");
  if(!h) {
    cout << __func__ << ": Input histogram for selection " << selection << " not found\n";
  } else {
    h = trim_hist(h, obs.getMin(), obs.getMax());
    const double scale = get_dataset_info(process).norm();
    h->Scale(scale);
    eff = h->Integral();
    cout << __func__ << ": Signal efficiency = " << eff << endl;
  }

  pdf_info res;
  RooAbsPdf* pdf(nullptr);
  res.pdf_   = nullptr;
  res.rate_  = eff;
  res.color_ = kBlue;
  res.name_  = "signal";
  res.title_ = (is_mumem) ? "#mu^{-}#rightarrowe^{-}" : "#mu^{-}#rightarrowe^{+}";

  if(!use_hist) {
    delete h; //no longer needed

    const int fit_version = 1; // 0: CB * Landau; 1: Landau CB; 2: CB
    const float signal_peak = (is_mumem) ? 104.0f : 92.3f;

    if(fit_version == 0) { // convolve double-sided crysal ball with energy losses

      // Make a RooLandau PDF, but flip the x-axis to represent energy losses
      auto f = new TF1("flip_signal_input","-x", obs.getMin(), obs.getMax());
      RooAbsReal* sig_flip_obs = RooFit::bindFunction(f,obs);

      RooRealVar* sig_mean     = new RooRealVar(Form("%s_sig_mean"  , name), "mean", -1*signal_peak, -1*signal_peak - 5., -1*signal_peak + 5.);
      RooRealVar* sig_sigma    = new RooRealVar(Form("%s_sig_sigma" , name), "sigma", 0.2, 0., 5.);
      RooAbsPdf*  sig_pdf      = new RooLandau (Form("%s_sig_pdf"  , name), "signal base PDF", *sig_flip_obs, *sig_mean, *sig_sigma);

      // Make a resolution function
      RooRealVar* res_mean     = new RooRealVar(Form("%s_res_mean"  , name), "mean", 0., -1., 1.);
      RooRealVar* res_sigma    = new RooRealVar(Form("%s_res_sigma" , name), "sigma", 0.2, 0., 5.);
      RooRealVar* res_alpha1   = new RooRealVar(Form("%s_res_alpha1", name), "alpha1", 1., 0.1, 10.);
      RooRealVar* res_alpha2   = new RooRealVar(Form("%s_res_alpha2", name), "alpha2", 1., 0.1, 10.);
      RooRealVar* res_n1       = new RooRealVar(Form("%s_res_n1"    , name), "enne1", 5., 0.1, 30.);
      RooRealVar* res_n2       = new RooRealVar(Form("%s_res_n2"    , name), "enne2", 5., 5., 30.);
      // RooFormulaVar* mean_func = new RooFormulaVar(Form("%s_mean_func", name), "mean with offset",
      //                                              "@0*(1 + @1*@2 + @3*@4)", RooArgList(*mean, *elec_ES_shift, *elec_ES_size, *muon_ES_shift, *muon_ES_size));
      RooAbsPdf* res_pdf       = new RooCrystalBall(Form("%s_res_pdf"  , name), "signal resolution PDF", obs, *res_mean, *res_sigma, *res_alpha1, *res_n1, *res_alpha2, *res_n2);

      // Convolve the PDF with the resolution
      obs.setBinning(RooBinning(10000, obs.getMin(), obs.getMax()), "cache"); // for discrete convolution
      pdf           = new RooFFTConvPdf(Form("%s_pdf", name), "Signal PDF", obs, *sig_pdf, *res_pdf, 2);
      ((RooFFTConvPdf*) pdf)->setBufferFraction(1.);
      ((RooFFTConvPdf*) pdf)->setBufferStrategy(RooFFTConvPdf::Extend); //Extend, Flat, or Mirror

      if(is_mumem) {
        sig_mean  ->setVal(-104.167 );
        sig_sigma ->setVal(0.205701 );
        res_alpha1->setVal(0.104305 );
        res_alpha2->setVal(0.353623 );
        res_mean  ->setVal(0.360386 );
        res_n1    ->setVal(3. );
        res_n2    ->setVal(10.  );
        res_sigma ->setVal(0.0147911);
      } else {
        if(selection == 40) {
          sig_mean  ->setVal( -91.4263 ); // +/- 0.0241985
          sig_sigma ->setVal( 0.17123  ); // +/- 0.00130733
          res_alpha1->setVal( 0.315472 ); // +/- 0.0838976
          res_alpha2->setVal( 1.0473   ); // +/- 0.251664
          res_mean  ->setVal( 0.437593 ); // +/- 0.0280891
          res_n1    ->setVal( 30       ); // +/- 0.672383
          res_n2    ->setVal( 5.39243  ); // +/- 0.297799
          res_sigma ->setVal( 0.0772093); // +/- 0.0202807
        }
      }
      sig_mean   ->setConstant(freeze);
      sig_sigma  ->setConstant(freeze);
      res_mean   ->setConstant(freeze);
      res_sigma  ->setConstant(freeze);
      res_alpha1 ->setConstant(freeze);
      res_alpha2 ->setConstant(freeze);
      res_n1     ->setConstant(freeze);
      res_n2     ->setConstant(freeze);
    } else if(fit_version == 1) { // Landau core + power-law tails
      // Make a RooLandauCB PDF

      RooRealVar* sig_mean     = new RooRealVar(Form("%s_sig_mean"  , name), "mean", signal_peak-1., signal_peak - 5., signal_peak + 5.);
      RooRealVar* sig_a        = new RooRealVar(Form("%s_sig_a"     , name), "a", 0.3, 0.1,  1.);
      RooRealVar* sig_b        = new RooRealVar(Form("%s_sig_b"     , name), "b", 4.7, 0.1, 10.);
      RooRealVar* sig_alpha1   = new RooRealVar(Form("%s_sig_alpha1", name), "alpha1", 1.5, 0.1, 2.);
      RooRealVar* sig_alpha2   = new RooRealVar(Form("%s_sig_alpha2", name), "alpha2", 0.5, 0.1, 3.);
      RooRealVar* sig_n1       = new RooRealVar(Form("%s_sig_n1"    , name), "n1", 1.5, 0.1, 10.);
      RooRealVar* sig_n2       = new RooRealVar(Form("%s_sig_n2"    , name), "n2", 5.0, 0.1, 10.);
      pdf                      = new RooLandauCB(Form("%s_pdf"      , name), "signal PDF", obs, *sig_mean, *sig_a, *sig_b, *sig_alpha1, *sig_n1, *sig_alpha2, *sig_n2);

      sig_mean   ->setConstant(freeze);
      sig_a      ->setConstant(freeze);
      sig_b      ->setConstant(freeze);
      sig_mean   ->setConstant(freeze);
      sig_alpha1 ->setConstant(freeze);
      sig_alpha2 ->setConstant(freeze);
      sig_n1     ->setConstant(freeze);
      sig_n2     ->setConstant(freeze);
    } else if(fit_version == 2) { // Double-sided Crystal Ball

      RooRealVar* sig_mean     = new RooRealVar(Form("%s_mean"  , name), "mean", signal_peak, signal_peak - 5., signal_peak + 5.);
      RooRealVar* sig_sigma    = new RooRealVar(Form("%s_sigma" , name), "sigma", 1., 0., 5.);
      RooRealVar* sig_alpha1   = new RooRealVar(Form("%s_alpha1", name), "alpha1", 1., 0.1, 10.);
      RooRealVar* sig_alpha2   = new RooRealVar(Form("%s_alpha2", name), "alpha2", 1., 0.1, 10.);
      RooRealVar* sig_n1       = new RooRealVar(Form("%s_n1"    , name), "enne1", 2., 0.1, 30.);
      RooRealVar* sig_n2       = new RooRealVar(Form("%s_n2"    , name), "enne2", 5., 0.1, 30.);
      pdf       = new RooCrystalBall(Form("%s_pdf"  , name), "Signal PDF", obs, *sig_mean, *sig_sigma, *sig_alpha1, *sig_n1, *sig_alpha2, *sig_n2);

      sig_mean   ->setConstant(freeze);
      sig_sigma  ->setConstant(freeze);
      sig_alpha1 ->setConstant(freeze);
      sig_alpha2 ->setConstant(freeze);
      sig_n1     ->setConstant(freeze);
      sig_n2     ->setConstant(freeze);
    }
  } else { //Use a histogram-based PDF
    auto signal_hist = new RooDataHist(Form("%s_data", name), "Signal Data", obs, h);
    pdf = new RooHistPdf(Form("%s_pdf", name), "Signal PDF", obs, *signal_hist);
  }

  res.pdf_   = pdf;

  return res;
}

#endif
