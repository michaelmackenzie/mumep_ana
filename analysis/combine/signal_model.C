#ifndef __CONVANA_ANALYSIS_SIGNALMODEL__
#define __CONVANA_ANALYSIS_SIGNALMODEL__

#include "../tools/types.C"
#include "../defaults.C"
#include "../physics.C"
#include "../datasets.C"
#include "../tools/utilities.C"
#include "RooTFnBinding.h"

//---------------------------------------------------------------------------------------------------------------------------
TH1* get_signal_hist(const TString process, const int selection, const TString name = "signal", const int isys = -1) {
  auto info = get_dataset_info(process);
  const bool is_mumem = process == "mumem";

  // Retrieve the input file
  TFile* f = TFile::Open(Form("%sConvAna.%s.%s.m%i.%s", hist_path_, hist_func_, info.name_.Data(), hist_mode_, file_type_.Data()), "READ");
  if(!f) return nullptr;

  // Retrieve the input histogram
  TString hist_name = (isys < 0) ? Form("%sHist/trk_%i/%s", dir_path_.Data(), selection, var_.Data()) : Form("%sHist/sys_%i/%s_%i", dir_path_.Data(), selection, var_.Data(), isys);
  TH1* h = (TH1*) f->Get(hist_name.Data());
  if(!h) {
    cout << __func__ << ": Input histogram for selection " << selection << " not found in file " << f->GetName()
         << ": Hist name = " << hist_name.Data() << endl;
    f->Close();
    return nullptr;
  }
  h = (TH1*) h->Clone(Form("%s_%i_%s%s", process.Data(), selection, name.Data(), (isys < 0) ? "" : Form("_sys_%i", isys)));
  h->SetDirectory(0);
  const int rebin = bin_width_ / h->GetBinWidth(1) + 1.e-3;
  if(rebin > 1) h->Rebin(rebin);

  // check the process N(event) counts
  TTree* t_norm = (TTree*) f->Get(Form("%sdata/Norm", dir_path_.Data()));
  if(!t_norm) cout << __func__ << ": Normalization tree for process " << process.Data() << " not found\n";
  else {
    Long64_t nseen(0), ntotal(0);
    t_norm->SetBranchAddress("nseen", &nseen);
    for(Long64_t entry = 0; entry < t_norm->GetEntries(); ++entry) {
      t_norm->GetEntry(entry);
      cout << "N(seen) = " << nseen << " N(total) = " << ntotal << endl;
      ntotal += nseen;
    }
    if(ntotal == 0) {
      cout << __func__ << ": No normalization contained in the normalization tree for process " << process.Data() << endl;
    } else {
      Long64_t nexpect = info.ndigi_;
      if(nexpect == 0) {
        cout << __func__ << ": No estimate for the number of expected events for process " << process.Data() << endl;
      } else {
        if(nexpect != ntotal) {
          const double ratio = nexpect * 1. / ntotal;
          cout << __func__ << ": See " << ntotal << " events but expect " << nexpect << " for process " << process.Data()
               << " --> scaling by " << ratio << endl;
          h->Scale(ratio);
        }
      }
    }
  }

  f->Close();
  return h;
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

    const int fit_version = 0; // 0: CB + Landau; 1: Landau CB
    const float signal_peak = (is_mumem) ? 105.0f : 92.3f;

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
      RooRealVar* res_n2       = new RooRealVar(Form("%s_res_n2"    , name), "enne2", 5., 0.1, 30.);
      // RooFormulaVar* mean_func = new RooFormulaVar(Form("%s_mean_func", name), "mean with offset",
      //                                              "@0*(1 + @1*@2 + @3*@4)", RooArgList(*mean, *elec_ES_shift, *elec_ES_size, *muon_ES_shift, *muon_ES_size));
      RooAbsPdf* res_pdf       = new RooCrystalBall(Form("%s_res_pdf"  , name), "signal resolution PDF", obs, *res_mean, *res_sigma, *res_alpha1, *res_n1, *res_alpha2, *res_n2);

      // Convolve the PDF with the resolution
      obs.setBinning(RooBinning(10000, obs.getMin(), obs.getMax()), "cache"); // for discrete convolution
      pdf           = new RooFFTConvPdf(Form("%s_pdf", name), "Signal PDF", obs, *sig_pdf, *res_pdf, 2);
      ((RooFFTConvPdf*) pdf)->setBufferFraction(1.);
      ((RooFFTConvPdf*) pdf)->setBufferStrategy(RooFFTConvPdf::Extend); //Extend, Flat, or Mirror

      if(is_mumem) {
        if(selection == 20) {
          sig_mean  ->setVal(-104.168 ); // +/- 0.0030234
          sig_sigma ->setVal(0.190398 ); // +/- 0.000811216
          res_alpha1->setVal(0.110681 ); // +/- 0.00466586
          res_alpha2->setVal(0.353623 ); // +/- 0.0101467
          res_mean  ->setVal(0.359745 ); // +/- 0.00294892
          res_n1    ->setVal(14.035   ); // +/- 4.40835
          res_n2    ->setVal(3.75875  ); // +/- 0.0546084
          res_sigma ->setVal(0.0147911); // +/- 0.000412189
        }
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
    // } else if(fit_version == 1) { // Landau core + power-law tails
    //   // Make a RooLandauCB PDF

    //   RooRealVar* sig_mean     = new RooRealVar(Form("%s_sig_mean"  , name), "mean", signal_peak-0.5, signal_peak - 5., signal_peak + 5.);
    //   RooRealVar* sig_a        = new RooRealVar(Form("%s_sig_a"     , name), "a", 0.3, 0.1,  1.);
    //   RooRealVar* sig_b        = new RooRealVar(Form("%s_sig_b"     , name), "b", 4.7, 0.1, 10.);
    //   RooRealVar* sig_alpha1   = new RooRealVar(Form("%s_sig_alpha1", name), "alpha1", 1.5, 0.1, 2.);
    //   RooRealVar* sig_alpha2   = new RooRealVar(Form("%s_sig_alpha2", name), "alpha2", 0.5, 0.1, 2.);
    //   RooRealVar* sig_n1       = new RooRealVar(Form("%s_sig_n1"    , name), "n1", 1.5, 0.1, 10.);
    //   RooRealVar* sig_n2       = new RooRealVar(Form("%s_sig_n2"    , name), "n2", 4.0, 0.1, 10.);
    //   pdf                      = new RooLandauCB(Form("%s_pdf"      , name), "signal PDF", obs, *sig_mean, *sig_a, *sig_b, *sig_alpha1, *sig_n1, *sig_alpha2, *sig_n2);

    //   if(is_mumem) {
    //     if(selection == 20) {
    //       sig_a     ->setVal(  0.306917); // +/- 0.503880
    //       sig_alpha1->setVal(  1.289232); // +/- 0.969600
    //       sig_alpha2->setVal(  0.411932); // +/- 1.076575
    //       sig_b     ->setVal(  4.764794); // +/- 7.068967
    //       sig_mean  ->setVal(104.350141); // +/- 0.365222
    //       sig_n1    ->setVal(  1.533756); // +/- 6.782118
    //       sig_n2    ->setVal(  4.756923); // +/- 7.026505
    //     }
    //   } else {
    //     if(selection == 40) {
    //       sig_a     ->setVal( 0.331757); // +/- 0.583273
    //       sig_alpha1->setVal( 1.513254); // +/- 0.966809
    //       sig_alpha2->setVal( 0.406306); // +/- 1.069646
    //       sig_b     ->setVal( 4.553921); // +/- 6.647130
    //       sig_mean  ->setVal(91.778857); // +/- 0.440913
    //       sig_n1    ->setVal( 1.925420); // +/- 8.981488
    //       sig_n2    ->setVal( 9.995132); // +/- 9.348596
    //     }
    //   }
    //   sig_mean   ->setConstant(freeze);
    //   sig_a      ->setConstant(freeze);
    //   sig_b      ->setConstant(freeze);
    //   sig_mean   ->setConstant(freeze);
    //   sig_alpha1 ->setConstant(freeze);
    //   sig_alpha2 ->setConstant(freeze);
    //   sig_n1     ->setConstant(freeze);
    //   sig_n2     ->setConstant(freeze);
    }
  } else { //Use a histogram-based PDF
    auto signal_hist = new RooDataHist(Form("%s_data", name), "Signal Data", obs, h);
    pdf = new RooHistPdf(Form("%s_pdf", name), "Signal PDF", obs, *signal_hist);
  }

  res.pdf_   = pdf;

  return res;
}

#endif
