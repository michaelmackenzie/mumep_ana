//Plot data from many  Higgs Combine FitDiagnostics toys
//
// file_name     : FitDiagnostics output file (uses tree_fit_sb)
// r_true        : signal strength the toys were generated with
// out_name      : output figure name (".png" appended if missing)
// skip_bad_fits : 0 = use all fits; 1 = also drop fits with asymmetric errors (ratio > 2);
//                 2 = also drop fits with non-zero fit status
// err_mode      : 0 = asymmetric errors; 1 = use rHiErr for both sides; 2 = use symmetric rErr
//
// The pull uses the error on the side facing r_true: (r - r_true)/rLoErr if r > r_true, else /rHiErr.
// Toys whose error on that side is not usable (zero, non-finite, or with the interval
// edge at the rMin/rMax boundary, i.e. the crossing was not found) are excluded from the
// pull but kept in the r distribution, and are counted on the plot.

int plot_combine_fits(const char* file_name, double r_true = 0., TString out_name = "",
                      const int skip_bad_fits = 0, const int err_mode = 0) {

  /////////////////////////////////////////////////////////////////
  // Retrieve the fit data

  TFile* file = TFile::Open(file_name, "READ");
  if(!file || file->IsZombie()) return 1;

  TTree* tree = (TTree*) file->Get("tree_fit_sb");
  if(!tree) {
    cout << "Tree \"tree_fit_sb\" not found in file " << file_name << endl;
    file->Close();
    return 1;
  }

  //Retrieve the r fit range to identify intervals truncated at the boundary
  double r_min = -std::numeric_limits<double>::infinity();
  double r_max =  std::numeric_limits<double>::infinity();
  RooFitResult* fit_res = (RooFitResult*) file->Get("fit_s");
  if(fit_res) {
    RooRealVar* r_var = (RooRealVar*) fit_res->floatParsFinal().find("r");
    if(r_var) {r_min = r_var->getMin(); r_max = r_var->getMax();}
  }
  const double r_tol = (std::isfinite(r_min) && std::isfinite(r_max)) ? 1.e-4*(r_max - r_min) : 0.;

  /////////////////////////////////////////////////////////////////
  // Loop through the fit results

  const Long64_t nentries = tree->GetEntries();
  double r, rErr, rHiErr, rLoErr;
  int fit_status;
  tree->SetBranchAddress("r"         , &r         );
  tree->SetBranchAddress("rErr"      , &rErr      );
  tree->SetBranchAddress("rHiErr"    , &rHiErr    );
  tree->SetBranchAddress("rLoErr"    , &rLoErr    );
  tree->SetBranchAddress("fit_status", &fit_status);

  std::vector<double> r_vals, pulls;
  int n_skipped = 0, n_bad_err = 0, n_r_at_bound = 0;
  for(Long64_t entry = 0; entry < nentries; ++entry) {
    tree->GetEntry(entry);
    if(nentries <= 10) cout << " Entry " << entry << ": (r, r_up, r_down) = (" << r << ", " << rHiErr << ", " << rLoErr << "), status = " << fit_status << endl;

    //over-ride the asymmetric errors if requested
    if(err_mode == 1) rLoErr = rHiErr;
    if(err_mode == 2) rLoErr = rHiErr = rErr;

    //check for fit issues
    if(skip_bad_fits > 0 && rLoErr > 0. && rHiErr > 0.) {
      if(std::fabs(rLoErr/rHiErr) > 2. || std::fabs(rHiErr/rLoErr) > 2.) {
        cout << " Entry " << entry << " has suspicious errors: (r, r_up, r_down) = (" << r << ", " << rHiErr << ", " << rLoErr << "), status = " << fit_status << endl;
        ++n_skipped;
        continue;
      }
    }
    if(fit_status != 0 && skip_bad_fits > 1) {
      cout << " Entry " << entry << " has non-zero fit status: (r, r_up, r_down) = (" << r << ", " << rHiErr << ", " << rLoErr << "), status = " << fit_status << endl;
      ++n_skipped;
      continue;
    }
    if(r - r_min < r_tol || r_max - r < r_tol) ++n_r_at_bound;
    r_vals.push_back(r);

    //Pull using the error on the side facing the true value
    const bool   above   = r > r_true;
    const double err     = (above) ? rLoErr : rHiErr;
    const bool   at_edge = (err_mode != 2) && ((above) ? (r - err - r_min < r_tol) : (r_max - (r + err) < r_tol));
    if(!std::isfinite(err) || err <= 0. || at_edge) {
      cout << " Entry " << entry << " has an unusable " << ((above) ? "low" : "high") << " error: (r, r_up, r_down) = ("
           << r << ", " << rHiErr << ", " << rLoErr << "), status = " << fit_status << endl;
      ++n_bad_err;
      continue;
    }
    pulls.push_back((r - r_true) / err);
  }
  file->Close();

  const int n_r = r_vals.size(), n_pull = pulls.size();
  if(n_r == 0) {
    cout << "No usable fits found in " << file_name << endl;
    return 1;
  }

  /////////////////////////////////////////////////////////////////
  // Summary statistics (unbinned)

  auto quantiles = [](std::vector<double> vals, std::vector<double> probs) {
    std::vector<double> res(probs.size());
    TMath::Quantiles(vals.size(), probs.size(), vals.data(), res.data(), probs.data(), false);
    return res;
  };
  const std::vector<double> r_q = quantiles(r_vals, {0.005, 0.025, 0.16, 0.5, 0.84, 0.975, 0.995});
  const double r_mean  = TMath::Mean(r_vals.begin(), r_vals.end());
  const double r_rms   = (n_r > 1) ? TMath::StdDev(r_vals.begin(), r_vals.end()) : 0.;
  const double r_68    = 0.5*(r_q[4] - r_q[2]);
  const double p_med   = (n_pull > 0) ? quantiles(pulls, {0.5})[0] : 0.;
  const double p_mean  = (n_pull > 0) ? TMath::Mean(pulls.begin(), pulls.end()) : 0.;
  const double p_rms   = (n_pull > 1) ? TMath::StdDev(pulls.begin(), pulls.end()) : 0.;

  /////////////////////////////////////////////////////////////////
  // Fill the histograms

  //Range covers the central 99% of toys, padded and including r_true
  double x_lo = std::min(r_q[0], r_true), x_hi = std::max(r_q[6], r_true);
  const double pad = std::max(0.1*(x_hi - x_lo), 1.e-6);
  x_lo -= pad; x_hi += pad;
  const int nbins = std::max(10, std::min(50, (int) std::sqrt(n_r)));
  TH1* h     = new TH1D("hr"   , "", nbins, x_lo, x_hi);
  TH1* hpull = new TH1D("hpull", "", 40, -4., 4.);
  for(double v : r_vals) h->Fill(v);
  for(double v : pulls ) hpull->Fill(v);
  const int n_out_r    = h    ->GetBinContent(0) + h    ->GetBinContent(nbins+1);
  const int n_out_pull = hpull->GetBinContent(0) + hpull->GetBinContent(41);

  /////////////////////////////////////////////////////////////////
  // Plot the results

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas* c = new TCanvas("c","c",1200,600);
  c->Divide(2,1);

  auto make_text = [](double x1, double y1, double x2, double y2) {
    TPaveText* t = new TPaveText(x1, y1, x2, y2, "NDC");
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    t->SetTextFont(42);
    t->SetTextSize(0.033);
    t->SetTextAlign(12);
    return t;
  };

  //Plot the fit result histogram
  c->cd(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.04);
  h->SetLineWidth(3);
  h->SetLineColor(kBlue);
  // h->SetFillColor(kBlue);
  // h->SetFillStyle(3003);
  h->SetXTitle("#hat{#mu}");
  h->SetYTitle("N(toys)");
  h->GetYaxis()->SetTitleOffset(1.3);
  const double y_top = h->GetMaximum();
  h->GetYaxis()->SetRangeUser(0., 1.5*y_top);
  h->Draw("hist");

  //Quantile bands: central 68% and 95% of toys
  TBox* band_2 = new TBox(r_q[1], 0., r_q[5], y_top);
  TBox* band_1 = new TBox(r_q[2], 0., r_q[4], y_top);
  band_2->SetFillColorAlpha(kOrange  , 0.7);
  band_1->SetFillColorAlpha(kSpring-5, 0.7);
  band_2->Draw();
  band_1->Draw();
  h->Draw("hist same");
  TLine* l_med  = new TLine(r_q[3], 0., r_q[3], y_top);
  TLine* l_true = new TLine(r_true, 0., r_true, y_top);
  l_med ->SetLineWidth(3); l_med ->SetLineStyle(kDashed); l_med->SetLineColor(kBlack);
  l_true->SetLineWidth(3); l_true->SetLineColor(kRed);
  l_med ->Draw();
  l_true->Draw();
  gPad->RedrawAxis();

  TLegend* leg = new TLegend(0.16, 0.66, 0.46, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.033);
  leg->AddEntry(h     , "Toy fits"   , "l");
  leg->AddEntry(l_true, "#mu_{true}" , "l");
  leg->AddEntry(l_med , "Median"     , "l");
  leg->AddEntry(band_1, "68% of toys", "f");
  leg->AddEntry(band_2, "95% of toys", "f");
  leg->Draw();

  TPaveText* t_r = make_text(0.52, 0.68, 0.95, 0.89);
  t_r->AddText(Form("N(toys) = %d", n_r));
  t_r->AddText(Form("Median = %.3g  (#mu_{true} = %.3g)", r_q[3], r_true));
  t_r->AddText(Form("68%%: [%.3g, %.3g]", r_q[2], r_q[4]));
  t_r->AddText(Form("Mean = %.3g #pm %.2g", r_mean, r_rms/std::sqrt(n_r)));
  // t_r->AddText(Form("(Median - #mu_{true}) / #sigma_{68} = %.3f", (r_68 > 0.) ? (r_q[3] - r_true)/r_68 : 0.));
  // if(n_out_r      > 0) t_r->AddText(Form("Outside plot range: %d", n_out_r));
  // if(n_r_at_bound > 0) t_r->AddText(Form("#hat{#mu} at r bound: %d", n_r_at_bound));
  t_r->Draw();

  //Plot the pull histogram
  c->cd(2);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.04);
  hpull->SetLineWidth(2);
  hpull->SetLineColor(kBlack);
  hpull->SetMarkerStyle(20);
  hpull->SetMarkerColor(kBlack);
  hpull->SetXTitle("(#hat{#mu} - #mu_{true}) / #sigma_{#mu}");
  hpull->SetYTitle("N(toys)");
  hpull->GetYaxis()->SetTitleOffset(1.3);
  hpull->GetYaxis()->SetRangeUser(0., 1.6*hpull->GetMaximum());
  hpull->Draw("E1 X0");

  //Reference: unit Gaussian normalized to the number of pulls
  TF1* f_ref = new TF1("f_ref", "[0]*TMath::Gaus(x, 0., 1., true)", -4., 4.);
  f_ref->SetParameter(0, n_pull*hpull->GetBinWidth(1));
  f_ref->SetLineColor(kGray+1);
  f_ref->SetLineStyle(kDashed);
  f_ref->SetLineWidth(2);
  f_ref->SetNpx(400);
  f_ref->Draw("same");

  TF1* f_fit = nullptr;
  if(n_pull > 10) {
    hpull->Fit("gaus", "LQ0");
    f_fit = hpull->GetFunction("gaus");
    if(f_fit) {
      f_fit->SetLineColor(kRed);
      f_fit->SetLineWidth(2);
      f_fit->SetNpx(400);
      f_fit->ResetBit(TF1::kNotDraw);
      f_fit->Draw("same");
    }
  }
  hpull->Draw("E1 X0 same");

  TLegend* leg_p = new TLegend(0.16, 0.74, 0.46, 0.89);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextSize(0.033);
  leg_p->AddEntry(hpull, "Toy pulls", "pe");
  if(f_fit) leg_p->AddEntry(f_fit, "Gaussian fit", "l");
  leg_p->AddEntry(f_ref, "N(0, 1)", "l");
  leg_p->Draw();

  TPaveText* t_p = make_text(0.52, 0.68, 0.95, 0.89);
  // t_p->AddText(Form("N(pulls) = %d", n_pull));
  if(f_fit) {
    t_p->AddText(Form("Fit mean = %.3f #pm %.3f", f_fit->GetParameter(1), f_fit->GetParError(1)));
    t_p->AddText(Form("Fit width = %.3f #pm %.3f", f_fit->GetParameter(2), f_fit->GetParError(2)));
  }
  t_p->AddText(Form("Median = %.3f", p_med));
  t_p->AddText(Form("Mean = %.3f, RMS = %.3f", p_mean, p_rms));
  // if(n_bad_err  > 0) t_p->AddText(Form("Unusable #sigma (excluded): %d", n_bad_err));
  // if(n_out_pull > 0) t_p->AddText(Form("|pull| > 4: %d", n_out_pull));
  // if(n_skipped  > 0) t_p->AddText(Form("Skipped fits: %d", n_skipped));
  t_p->Draw();

  //Print a summary
  printf("plot_combine_fits: %s\n", file_name);
  printf("  toys: %lld, used: %d, skipped: %d, pulls: %d (unusable error: %d), r at bound: %d\n",
         nentries, n_r, n_skipped, n_pull, n_bad_err, n_r_at_bound);
  printf("  r: median = %.4g, 68%% = [%.4g, %.4g], mean = %.4g, (median - r_true)/sigma_68 = %.4f\n",
         r_q[3], r_q[2], r_q[4], r_mean, (r_68 > 0.) ? (r_q[3] - r_true)/r_68 : 0.);
  printf("  pull: median = %.4f, mean = %.4f, rms = %.4f", p_med, p_mean, p_rms);
  if(f_fit) printf(", fit mean = %.4f +- %.4f, fit width = %.4f +- %.4f",
                   f_fit->GetParameter(1), f_fit->GetParError(1), f_fit->GetParameter(2), f_fit->GetParError(2));
  printf("\n");

  if(out_name == "") {out_name = file_name; out_name.ReplaceAll(".root", ".png");}
  else if(!out_name.EndsWith(".png")) out_name += ".png";
  c->SaveAs(out_name.Data());
  return 0;
}
