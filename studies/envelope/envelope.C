// Make an example plot with a signal and a background
#include "mumep_ana/analysis/tools/physics.C"
#include "mumep_ana/analysis/tools/functions.C"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"

using namespace RooFit;

void convolve_with_response(TF1* model, TH1* result) {
  if(!model || !result) return;
  result->Reset();
  TF1* res = landau_crystal_ball_tf1();
  res->SetName("res");
  res->SetRange(-15., 15.);
  res->SetParameters(1., -1., 0.5, 3.25, 0.7, 0.5, 2.5, 3.5);
  res->SetParameter(0, 1./res->Integral(-15., 15.));
  cout << "Convolving model " << model->GetName() << endl;

  // Convolve
  const double s_emin = model->GetXmin();
  const double s_emax = model->GetXmax();
  const double r_emin = res->GetXmin();
  const double r_emax = res->GetXmax();
  const double de = 0.01; // 1 keV steps
  double s = s_emin;
  while(s <= s_emax) {
    double r = r_emin;
    while(r <= r_emax) {
      const double p = res->Eval(r)*de*model->Eval(s)*de;
      if(p >= 0.) result->Fill(s + r, p);
      r += de;
    }
    s += de;
  }
  delete res;
  result->Scale(1./result->Integral()/result->GetBinWidth(1));
}

int envelope() {
  const double endpoint = calculate_ep_endpoint(13, 27, 0, 0); // 0 neutron and 0 proton knockout endpoint
  printf("%s: Mu- --> e+ endpoint: %.2f MeV\n", __func__, endpoint);

  // Ground state functions
  TF1* f_gs = new TF1("gs", "1", endpoint-0.001, endpoint+0.001);
  f_gs->SetNpx(1000);

  // Convolve with response
  const double xmin = 85.;
  const double xmax = 100.;
  const int nbins = (xmax - xmin)*4.;
  TH1* h_gs = new TH1D("h_gs", "", nbins, xmin, xmax);
  convolve_with_response(f_gs, h_gs);
  h_gs->Scale(100.);

  // Background PDF
  TF1* f_bkg = new TF1("bkg", "expo", xmin, xmax);
  f_bkg->SetParameters(1.,-0.25);

  // Generate toy data
  const int nbkg = 1e4;
  TH1* h_bkg = new TH1F("h_bkg", "", nbins, xmin, xmax);
  for(int ibkg = 0; ibkg < nbkg; ++ibkg) h_bkg->Fill(f_bkg->GetRandom());

  // ----------------------------------------------------
  // BLINDING PREPARATION (PRE-PROCESS TH1)
  // ----------------------------------------------------
  const double blind_low  = 89.5;
  const double blind_high = 92.5;

  // Zero out the content in the blinded region to enforce sideband-only entries
  for (int i = 1; i <= h_bkg->GetNbinsX(); ++i) {
    double bin_center = h_bkg->GetBinCenter(i);
    if (bin_center >= blind_low && bin_center <= blind_high) {
      h_bkg->SetBinContent(i, 0.0);
      h_bkg->SetBinError(i, 0.0);
    }
  }

  // Set stylistic parameters for the data histogram
  h_bkg->SetMarkerStyle(20);
  h_bkg->SetLineWidth(2);
  h_bkg->SetMarkerColor(kBlack);
  h_bkg->SetLineColor(kBlack);

  // ----------------------------------------------------
  // ROOFIT SETUP & STABLE ALTERNATE HYPOTHESIS FIT
  // ----------------------------------------------------
  RooRealVar e_pos("e_pos", "Positron momentum", xmin, xmax, "MeV/c");

  // Define sub-ranges explicitly for the Fitter calculation
  e_pos.setRange("LowSideband",  xmin, blind_low);
  e_pos.setRange("HighSideband", blind_high, xmax);
  e_pos.setRange("BlindedRegion", blind_low, blind_high);

  // Import the physically blinded TH1 into RooDataHist
  RooDataHist bkg_data("bkg_data", "Blinded Background Dataset", e_pos, h_bkg);

  // --- Model A: Exponential ---
  RooRealVar exp_slope("exp_slope", "Slope of exponential", -0.25, -0.6, -0.01);
  RooExponential bkg_exp("bkg_exp", "Exponential Bkg Model", e_pos, exp_slope);

  // --- Model B: Power Law (Alternate Option replacing Chebychev) ---
  // A power law spectrum (e_pos^(-alpha)) tracks falling shapes safely without sideband artifacts
  RooRealVar alpha("alpha", "Power law exponent", 22.0, 1.0, 50.0);
  RooGenericPdf bkg_pow("bkg_pow", "Power Law Model", "pow(e_pos, -alpha)", RooArgSet(e_pos, alpha));

  // Fit both stable shapes independently across the sidebands
  bkg_exp.fitTo(bkg_data, Range("LowSideband,HighSideband"), PrintLevel(-1));
  bkg_pow.fitTo(bkg_data, Range("LowSideband,HighSideband"), PrintLevel(-1));

  // ----------------------------------------------------
  // PLOTTING & VISUALIZATION
  // ----------------------------------------------------
  TCanvas* c = new TCanvas("c", "Blinded Fit Canvas", 800, 600);
  gStyle->SetOptStat(0);

  RooPlot* frame = e_pos.frame(Title(""));

  // Plot data points (empty bins in signal corridor naturally won't draw)
  bkg_data.plotOn(frame, MarkerStyle(20), Name("data"));

  // --- Draw Model A: Exponential ---
  bkg_exp.plotOn(frame, LineColor(kRed), LineStyle(kDashed), Name("exp_blinded"), Range("BlindedRegion"));
  bkg_exp.plotOn(frame, LineColor(kRed), LineStyle(kSolid), Name("exp_low_sb"), Range("LowSideband"));
  bkg_exp.plotOn(frame, LineColor(kRed), LineStyle(kSolid), Name("exp_high_sb"), Range("HighSideband"));

  // --- Draw Model B: Power Law ---
  bkg_pow.plotOn(frame, LineColor(kGreen+2), LineStyle(kDashed), Name("pow_blinded"), Range("BlindedRegion"));
  bkg_pow.plotOn(frame, LineColor(kGreen+2), LineStyle(kSolid), Name("pow_low_sb"), Range("LowSideband"));
  bkg_pow.plotOn(frame, LineColor(kGreen+2), LineStyle(kSolid), Name("pow_high_sb"), Range("HighSideband"));

  frame->Draw();

  // Draw signal ground state histogram
  h_gs->SetLineColor(kBlue);
  h_gs->SetLineWidth(3);
  h_gs->Draw("hist same");

  // Add structured legend matching the new stable choices
  TLegend* leg = new TLegend(0.40, 0.62, 0.88, 0.88);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(frame->findObject("data"), "Data (Sidebands Only)", "ep");
  leg->AddEntry(frame->findObject("exp_low_sb"), "Model A: Exponential (Fit)", "l");
  leg->AddEntry(frame->findObject("exp_blinded"), "Model A: Exponential (Blinded)", "l");
  leg->AddEntry(frame->findObject("pow_low_sb"), "Model B: Power Law (Fit)", "l");
  leg->AddEntry(frame->findObject("pow_blinded"), "Model B: Power Law (Blinded)", "l");
  leg->AddEntry(h_gs, "Convolved Signal Model (GS)", "l");
  leg->Draw();

  // Save the figure
  gSystem->Exec("mkdir -p figures");
  c->SaveAs("figures/envelope.png");

  return 0;
}
