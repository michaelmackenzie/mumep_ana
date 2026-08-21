// Build a toy background dataset and construct/plot a RooMultiPdf envelope.

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooFit.h"
#include "RooPlot.h"
#include "RooRealVar.h"

#if __has_include("combine/HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h")
#include "combine/HiggsAnalysis/CombinedLimit/src/RooMultiPdf.cxx"
#include "combine/HiggsAnalysis/CombinedLimit/src/RooBernsteinFast.cc"
#endif

using namespace RooFit;

// create_envelope.C expects these symbols in scope.
bool blind_data_ = false;
int verbose_ = 0;

int count_pdf_params(RooAbsPdf* pdf) {
  if(!pdf) return 0;
  RooArgSet obs;
  RooArgSet* params = pdf->getParameters(obs);
  if(!params) return 0;

  int nparams = 0;
  TIterator* iter = params->createIterator();
  TObject* obj = nullptr;
  while((obj = iter->Next())) {
    RooRealVar* var = dynamic_cast<RooRealVar*>(obj);
    if(var && !var->isConstant()) ++nparams;
  }

  delete iter;
  delete params;
  return nparams;
}

#include "../../analysis/combine/create_envelope.C"

int multi_pdf() {
  const double xmin = 85.;
  const double xmax = 100.;
  const int nbins = static_cast<int>((xmax - xmin) * 4.);

  // Generate a toy background sample.
  TF1* f_bkg = new TF1("f_bkg", "expo", xmin, xmax);
  f_bkg->SetParameters(1., -0.25);

  TH1* h_bkg = new TH1F("h_bkg_multi", "", nbins, xmin, xmax);
  const int nbkg = 10000;
  for(int i = 0; i < nbkg; ++i) h_bkg->Fill(f_bkg->GetRandom());

  // Blind the signal region to emulate sideband-only fitting.
  const double blind_low = 89.5;
  const double blind_high = 92.5;
  for(int ibin = 1; ibin <= h_bkg->GetNbinsX(); ++ibin) {
    const double center = h_bkg->GetBinCenter(ibin);
    if(center >= blind_low && center <= blind_high) {
      h_bkg->SetBinContent(ibin, 0.);
      h_bkg->SetBinError(ibin, 0.);
    }
  }

  h_bkg->SetMarkerStyle(20);
  h_bkg->SetMarkerColor(kBlack);
  h_bkg->SetLineColor(kBlack);

  RooRealVar e_pos("e_pos", "Positron momentum", xmin, xmax, "MeV/c");
  e_pos.setRange("LowSideband", xmin, blind_low);
  e_pos.setRange("HighSideband", blind_high, xmax);
  e_pos.setRange("BlindRegion", blind_low, blind_high);
  e_pos.setRange("full", xmin, xmax);

  RooDataHist bkg_data("bkg_data", "Toy blinded background", RooArgList(e_pos), h_bkg);

  double n_sideband = 0.;
  for(int ibin = 1; ibin <= h_bkg->GetNbinsX(); ++ibin) {
    const double center = h_bkg->GetBinCenter(ibin);
    if(center < blind_low || center > blind_high) n_sideband += h_bkg->GetBinContent(ibin);
  }

  // Configure envelope families.
  use_exp_family_ = true;
  use_power_family_ = true;
  use_laurent_family_ = false;
  use_inv_poly_family_ = false;
  use_poly_family_ = false;
  use_generic_bernstein_ = false;
  use_fast_bernstein_ = false;
  force_fit_order_ = false;
  test_single_function_ = false;
  add_all_fits_ = false;

  RooCategory pdf_index("pdf_index", "pdf index");
  RooMultiPdf* envelope = create_envelope(e_pos, pdf_index, bkg_data, true, "toy_env", 2);
  if(!envelope) {
    printf("%s: Failed to construct envelope.\n", __func__);
    delete h_bkg;
    delete f_bkg;
    return 1;
  }

  TCanvas* c = new TCanvas("c_multi_pdf", "Envelope", 900, 700);
  gStyle->SetOptStat(0);

  RooPlot* frame = e_pos.frame(Title("Toy Data with MultiPdf Envelope"));
  bkg_data.plotOn(frame, Name("data"), MarkerStyle(20));

  const int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kOrange + 7, kMagenta + 1, kCyan + 2, kViolet + 1, kTeal + 3};
  const int ncolors = sizeof(colors) / sizeof(colors[0]);

  const int npdfs = pdf_index.numTypes();
  for(int i = 0; i < npdfs; ++i) {
    pdf_index.setIndex(i);
    TString obj_name = Form("pdf_%i", i);
    envelope->plotOn(frame, Name(obj_name), LineColor(colors[i % ncolors]), LineWidth(2),
                     NormRange("LowSideband,HighSideband"),
                     Normalization(n_sideband, RooAbsReal::NumEvent));
  }

  frame->GetXaxis()->SetTitle("Positron momentum [MeV/c]");
  frame->GetYaxis()->SetTitle("Events / bin");
  frame->Draw();

  TLegend* leg = new TLegend(0.56, 0.56, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(frame->findObject("data"), "Toy data (blinded)", "ep");
  for(int i = 0; i < npdfs; ++i) {
    TString obj_name = Form("pdf_%i", i);
    leg->AddEntry(frame->findObject(obj_name), Form("Envelope PDF %i", i), "l");
  }
  leg->Draw();

  gSystem->Exec("mkdir -p figures");
  c->SaveAs("figures/multi_pdf.png");

  printf("%s: Saved figure to figures/multi_pdf.png with %i PDFs in envelope.\n", __func__, npdfs);
  return 0;
}
