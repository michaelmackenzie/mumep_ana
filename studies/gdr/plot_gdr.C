// Compare the muon to positron ground state and giant dipole resonance final states
#include "mumep_ana/analysis/tools/physics.C"
#include "mumep_ana/analysis/tools/functions.C"

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

int plot_gdr(int z = 13, int a = 27) {
  const double endpoint = calculate_ep_endpoint(z, a, 0, 0); // 0 neutron and 0 proton knockout endpoint
  printf("%s: Mu- --> e+ endpoint: %.2f MeV\n", __func__, endpoint);

  // GDR function
  const double width =  6.7; // MeV
  const double mode  = 21.1; // MeV
  TF1* f_gdr = new TF1("gdr", "[0]*[1]/((x - [2])^2 + [1]^2)", 0., endpoint);
  f_gdr->SetParameters(1., width/2., endpoint - mode);
  f_gdr->SetNpx(1000);
  f_gdr->SetParameter(0, 1./f_gdr->Integral(0., endpoint));

  // Ground state functions
  TF1* f_gs = new TF1("gs", "1", endpoint-0.001, endpoint+0.001);
  f_gs->SetNpx(1000);

  // Convolve with response
  TH1* h_gdr = new TH1D("h_gdr", "", 1000, 0., 100.);
  convolve_with_response(f_gdr, h_gdr);
  TH1* h_gs = new TH1D("h_gs", "", 1000, 0., 100.);
  convolve_with_response(f_gs, h_gs);

  // Set styles
  TLine* l_gs = new TLine(endpoint, 0., endpoint, 1.);
  h_gdr->SetLineColor( kRed); h_gdr->SetLineWidth(2);
  f_gdr->SetLineColor( kRed); f_gdr->SetLineStyle(kDashed);
  h_gs ->SetLineColor(kBlue); h_gs ->SetLineWidth(2);
  l_gs ->SetLineColor(kBlue); l_gs ->SetLineStyle(kDashed);

  // Draw the results
  TCanvas* c = new TCanvas();
  gStyle->SetOptStat(0);
  TH1* haxis = new TH1F("haxis", ";E_{+} (MeV);p(E_{+})", 1, 50., 100.);
  haxis->Draw();
  f_gdr->Draw("same");
  h_gdr->Draw("hist same");
  l_gs ->Draw("same");
  h_gs ->Draw("hist same");
  haxis->GetYaxis()->SetRangeUser(0., 1.);

  // Add a legend
  TLegend* leg = new TLegend(0.15, 0.7, 0.75, 0.89);
  leg->SetFillStyle(0); leg->SetBorderSize(0);
  leg->SetNColumns(2);
  leg->AddEntry(l_gs , "GS model", "L");
  leg->AddEntry(h_gs , "GS response", "L");
  leg->AddEntry(f_gdr, "GDR model", "L");
  leg->AddEntry(h_gdr, "GDR response", "L");
  leg->Draw();

  // Save the figure
  gSystem->Exec("mkdir -p figures");
  c->SaveAs("figures/spectrum.png");

  return 0;
}
