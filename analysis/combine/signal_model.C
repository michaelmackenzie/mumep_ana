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

    const int fit_version = 1; // 0: CB * Landau; 1: Landau CB; 2: CB; 3: Voigtian; 4: Cruijff
    const float signal_peak = (is_mumem) ? 104.0f : 92.3f;
    RooRealVar*    es_nuis      = new RooRealVar   (Form("%s_%i_es"  , process.Data(), selection), "Energy scale nuisance", 0., -7., 7.); es_nuis->setConstant(true);
    RooRealVar*    sig_mean_nom = new RooRealVar   (Form("%s_sig_mean_nom", name), "mean", signal_peak-1., signal_peak - 5., signal_peak + 5.);
    RooRealVar*    es_size      = new RooRealVar   (Form("%s_es_size"     , name), "Energy scale size", 0.1); es_size->setConstant(true);
    RooFormulaVar* sig_mean     = new RooFormulaVar(Form("%s_sig_mean"    , name), "@0+@1*@2", RooArgList(*sig_mean_nom, *es_size, *es_nuis));

    if(fit_version == 0) { // convolve a Crystal Ball with energy losses
      RooRealVar* l_width  = new RooRealVar(Form("%s_lWidth", name), "Landau Width", 0.3, 0.1, 1.0);
      RooFormulaVar* flipped_obs  = new RooFormulaVar(Form("%s_flippedObs", name), "-1.0 * @0", RooArgList(obs));
      RooFormulaVar* flipped_mean = new RooFormulaVar(Form("%s_flippedMean", name), "-1.0 * @0", RooArgList(*sig_mean));
      RooLandau* landau = new RooLandau(Form("%s_landau", name), "Energy Loss Physics", *flipped_obs, *flipped_mean, *l_width);

      RooRealVar* g_coreL   = new RooRealVar(Form("%s_gCoreL", name), "Resolution Core Left", 0.20, 0.05, 0.4);
      RooRealVar* g_coreR   = new RooRealVar(Form("%s_gCoreR", name), "Resolution Core Right", 0.20, 0.05, 0.4);
      RooRealVar* g_alphaL  = new RooRealVar(Form("%s_gAlphaL", name), "Left Tail Alpha", 1.5, 0.5, 3.5);
      RooRealVar* g_nL      = new RooRealVar(Form("%s_gNL",     name), "Left Tail Power N", 2.5, 1.1, 6.0);
      RooRealVar* g_alphaR  = new RooRealVar(Form("%s_gAlphaR", name), "Right Tail Alpha", 1.63, 0.5, 3.5);
      RooRealVar* g_nR      = new RooRealVar(Form("%s_gNR",     name), "Right Tail Power N", 4.0, 1.1, 6.0);
      RooCrystalBall* resolutionModel = new RooCrystalBall(Form("%s_resModel", name), "Asymmetric DSCB Resolution",
                                                           obs, RooConst(0), *g_coreL, *g_coreR,
                                                           *g_alphaL, *g_nL, *g_alphaR, *g_nR);

      pdf = new RooFFTConvPdf(Form("%s_pdf", name), "Landau (X) Asymmetric DSCB", obs, *landau, *resolutionModel);
      ((RooFFTConvPdf*) pdf)->setBufferFraction(2.0);

      sig_mean_nom ->setConstant(freeze);
      l_width      ->setConstant(freeze);
      g_coreL      ->setConstant(freeze);
      g_coreR      ->setConstant(freeze);
      g_alphaL     ->setConstant(freeze);
      g_nL         ->setConstant(freeze);
      g_alphaR     ->setConstant(freeze);
      g_nR         ->setConstant(true  );
    } else if(fit_version == 1) { // Landau core + power-law tails

      // Make a RooLandauCB PDF
      sig_mean_nom->setVal(signal_peak - 0.7);
      RooRealVar* sig_a        = new RooRealVar(Form("%s_a"       , name), "a", 0.5, 0.1,  1.);
      RooRealVar* sig_b        = new RooRealVar(Form("%s_b"       , name), "b", 2.7217, 0.1, 10.);
      RooRealVar* sig_alpha1   = new RooRealVar(Form("%s_alpha1"  , name), "alpha1", 0.7542, 0.3, 2.);
      RooRealVar* sig_alpha2   = new RooRealVar(Form("%s_alpha2"  , name), "alpha2", 0.5561, 0.3, 2.);
      RooRealVar* sig_n1       = new RooRealVar(Form("%s_n1"      , name), "n1", 2.6869, 0.1, 20.);
      RooRealVar* sig_n2       = new RooRealVar(Form("%s_n2"      , name), "n2", 3.2639, 0.1, 20.);
      pdf                      = new RooLandauCB(Form("%s_pdf"    , name), "signal PDF", obs, *sig_mean, *sig_a, *sig_b, *sig_alpha1, *sig_n1, *sig_alpha2, *sig_n2);

      sig_mean_nom ->setConstant(freeze);
      sig_a        ->setConstant(true  );
      sig_b        ->setConstant(freeze);
      sig_alpha1   ->setConstant(freeze);
      sig_alpha2   ->setConstant(freeze);
      sig_n1       ->setConstant(freeze);
      sig_n2       ->setConstant(freeze);
    } else if(fit_version == 2) { // Double-sided Crystal Ball

      // Asymmetric central Gaussian
      RooRealVar* sig_sigmaL = new RooRealVar(Form("%s_sigmaL", name), "Width Left", 0.6, 0.2, 1.2);
      RooRealVar* sig_sigmaR = new RooRealVar(Form("%s_sigmaR", name), "Width Right", 0.22, 0.1, 0.4);
      RooRealVar* sig_alphaL = new RooRealVar(Form("%s_alphaL", name), "Tail Transition Left", 1.2, 0.4, 2.5);
      RooRealVar* sig_nL     = new RooRealVar(Form("%s_nL",     name), "Power Law Left", 2.0, 1.1, 5.0);
      RooRealVar* sig_alphaR = new RooRealVar(Form("%s_alphaR", name), "Tail Transition Right", 2.2, 1.2, 3.5);
      RooRealVar* sig_nR     = new RooRealVar(Form("%s_nR",     name), "Power Law Right", 2.5, 1.2, 5.0);
      pdf = new RooCrystalBall(Form("%s_pdf", name), "Asymmetric Generalized DSCB", obs,
                               *sig_mean, *sig_sigmaL, *sig_sigmaR,
                               *sig_alphaL, *sig_nL, *sig_alphaR, *sig_nR);

      sig_mean_nom ->setConstant(freeze);
      sig_sigmaL   ->setConstant(freeze);
      sig_sigmaR   ->setConstant(freeze);
      sig_alphaL   ->setConstant(freeze);
      sig_nL       ->setConstant(freeze);
      sig_alphaR   ->setConstant(freeze);
      sig_nR       ->setConstant(freeze);

      // Symmetric central Gaussian version:
      //   sig_mean_nom->setVal(signal_peak - 0.1);
      //   RooRealVar* sig_sigma    = new RooRealVar(Form("%s_sigma" , name), "sigma", 0.6, 0., 5.);
      //   RooRealVar* sig_alpha1   = new RooRealVar(Form("%s_alpha1", name), "alpha1", 0.5, 0.1, 10.);
      //   RooRealVar* sig_alpha2   = new RooRealVar(Form("%s_alpha2", name), "alpha2", 1.4, 0.1, 10.);
      //   RooRealVar* sig_n1       = new RooRealVar(Form("%s_n1"    , name), "enne1", 2.5, 0.1, 10.);
      //   RooRealVar* sig_n2       = new RooRealVar(Form("%s_n2"    , name), "enne2", 2.5, 0.1, 10.);
      //   pdf       = new RooCrystalBall(Form("%s_pdf"  , name), "Signal PDF", obs, *sig_mean, *sig_sigma, *sig_alpha1, *sig_n1, *sig_alpha2, *sig_n2);

      //   sig_mean_nom ->setConstant(freeze);
      //   sig_sigma    ->setConstant(freeze);
      //   sig_alpha1   ->setConstant(freeze);
      //   sig_alpha2   ->setConstant(freeze);
      //   sig_n1       ->setConstant(freeze);
      //   sig_n2       ->setConstant(freeze);
    } else if(fit_version == 3) { // Voigtian
      RooRealVar* sig_a        = new RooRealVar(Form("%s_a"       , name), "Gamma", 1.0, 0.1, 5.0);
      RooRealVar* sig_b        = new RooRealVar(Form("%s_b"       , name), "Sigma", 0.2, 0.05, 1.0);
      pdf                      = new RooVoigtian(Form("%s_pdf"   , name), "Voigtian signal PDF", obs, *sig_mean, *sig_a, *sig_b);

      sig_mean_nom ->setConstant(freeze);
      sig_a        ->setConstant(freeze);
      sig_b        ->setConstant(freeze);
    } else if(fit_version == 4) { // Cruijff PDF
      RooRealVar* sig_a        = new RooRealVar(Form("%s_a", name), "sigmaL", 0.4, 0.05, 1.5);
      RooRealVar* sig_b        = new RooRealVar(Form("%s_b", name), "sigmaR", 0.22, 0.05, 0.5);
      RooRealVar* sig_alpha1   = new RooRealVar(Form("%s_alpha1", name), "alphaL", 0.01, 0.0001, 0.1);
      RooRealVar* sig_alpha2   = new RooRealVar(Form("%s_alpha2", name), "alphaR", 0.005, 0.0001, 0.05);
      sig_mean_nom->setVal(signal_peak - 0.9);

      // Inline mathematical definition of the Cruijff PDF
      TString formula = Form(
                             "exp(-0.5 * (%s - %s)*(%s - %s) / "
                             "((%s < %s) ? (%s*%s + %s*(%s - %s)*(%s - %s)) : "
                             "(%s*%s + %s*(%s - %s)*(%s - %s))))",
                             obs.GetName(), sig_mean->GetName(), obs.GetName(), sig_mean->GetName(),
                             obs.GetName(), sig_mean->GetName(), sig_a->GetName(), sig_a->GetName(), sig_alpha1->GetName(), obs.GetName(), sig_mean->GetName(), obs.GetName(), sig_mean->GetName(),
                             sig_b->GetName(), sig_b->GetName(), sig_alpha2->GetName(), obs.GetName(), sig_mean->GetName(), obs.GetName(), sig_mean->GetName()
                             );

      pdf = new RooGenericPdf(Form("%s_pdf", name), "Cruijff signal PDF", formula,
                              RooArgList(obs, *sig_mean, *sig_a, *sig_b, *sig_alpha1, *sig_alpha2));

      sig_mean_nom ->setConstant(freeze);
      sig_a        ->setConstant(freeze);
      sig_b        ->setConstant(freeze);
      sig_alpha1   ->setConstant(freeze);
      sig_alpha2   ->setConstant(freeze);
    } else { //Use a histogram-based PDF
      auto signal_hist = new RooDataHist(Form("%s_data", name), "Signal Data", obs, h);
      pdf = new RooHistPdf(Form("%s_pdf", name), "Signal PDF", obs, *signal_hist);
    }
  }

  res.pdf_   = pdf;

  return res;
}

#endif
