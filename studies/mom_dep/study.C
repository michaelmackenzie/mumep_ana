// Plot the momentum dependence of variables

int plot(TTree* t, TString var, TString cut = "", TString name = "", TString xvar = "trk_p") {
  if(name == "") name = var;
  TCanvas c;
  t->Draw(var+":"+xvar+" >> h2d", cut, "colz");
  TH2* h2d = (TH2*) gDirectory->Get("h2d");
  if(!h2d) return 1;
  c.SaveAs("figures/" + name + "_2d.png");

  // Normalize the columns
  for(int i = 1; i <= h2d->GetNbinsX(); ++i) {
    double norm = h2d->Integral(i, i, 1, h2d->GetNbinsY());
    if(norm == 0.) continue;
    for(int j = 1; j <= h2d->GetNbinsY(); ++j) {
      h2d->SetBinContent(i,j, h2d->GetBinContent(i,j)/norm);
    }
  }
  h2d->Draw("colz");
  c.SaveAs("figures/" + name + "_2d_norm.png");
  c.SetLogz();
  c.SaveAs("figures/" + name + "_2d_norm_log.png");
  c.SetLogz(false);

  // Draw slices
  const int xbins = h2d->GetNbinsX();
  const int rebin = xbins / 8;
  h2d->RebinX(rebin);
  bool first = true;
  TH1* haxis = nullptr;
  double max_val = 0.;
  TLegend leg(0.11, 0.72, 0.89, 0.89);
  for(int i = 1; i <= h2d->GetNbinsX(); ++i) {
    TH1* h = h2d->ProjectionY(Form("h_%i", i), i, i);
    if(h->Integral() == 0.) continue;
    h->SetLineColor((i > 4) ? i + 1 : i);
    h->Draw((first) ? "hist" : "hist same");
    if(first) haxis = h;
    first = false;
    max_val = max(max_val, h->GetMaximum());
    leg.AddEntry(h, Form("[%.1f,%.1f]", h2d->GetXaxis()->GetBinLowEdge(i), h2d->GetXaxis()->GetBinUpEdge(i)));
  }
  leg.SetFillColor(0); leg.SetLineWidth(0); leg.SetNColumns(3);
  leg.Draw();
  if(haxis) {
    haxis->GetYaxis()->SetRangeUser(0., 1.25*max_val);
    c.SaveAs("figures/" + name + "_slices.png");
    haxis->GetYaxis()->SetRangeUser(max_val/1.e4, 30.*max_val);
    c.SetLogy();
    c.SaveAs("figures/" + name + "_slices_log.png");
    c.SetLogy(false);
    delete haxis;
  }
  delete h2d;
  return 0;
}

int study(const char* file = "/exp/mu2e/data/projects/run1a/mumep_ana/histograms/ConvAna.cnv_ana.fpos0b1s5r0100.m1.root",
          const int set = 40) {

  TFile* f = TFile::Open(file, "READ");
  if(!f) return 1;
  const char* tpath = Form("Ana/Hist/trs_%i/tree", set);
  TTree* t = (TTree*) f->Get(tpath);
  if(!t) {
    cout << "Tree not found! Path = " << tpath << endl;
    f->ls();
    return 1;
  }

  gSystem->Exec("mkdir -p figures");
  gStyle->SetOptStat(0);
  plot(t, "trk_qual", "abs(trk_mc_dp) < 0.2", "trk_qual_good");
  plot(t, "trk_qual", "trk_mc_dp > 0.75", "trk_qual_bad");
  plot(t, "trk_altqual", "abs(trk_mc_dp) < 0.2", "trk_altqual_good");
  plot(t, "trk_altqual", "trk_mc_dp > 0.75", "trk_altqual_bad");
  plot(t, "trk_pid");
  plot(t, "trk_altpid");
  plot(t, "trk_trkonlypid");
  plot(t, "trk_active_ratio");
  plot(t, "trk_null_ratio");
  plot(t, "trk_rmax");
  plot(t, "trk_cos");
  plot(t, "trk_fitcon");
  plot(t, "trk_null_ratio", "", "trk_rmax_vs_null_ratio", "trk_rmax");

  return 0;
}
