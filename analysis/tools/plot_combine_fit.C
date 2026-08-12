//Plot fit results from a single FitDiagnostics fit
//create a fit diagnostics root file via:
//$> combine -M FitDiagnostics -d <input card/workspace> --saveShapes --saveWithUncertainties [additional options]

bool unblind_      = false  ;
int  err_mode_     =  0     ; //errors in the pulls: 0: sqrt(data^2 + fit^2); 1: sqrt(data^2 - fit^2)
int  print_stacks_ = -1     ; //print stacked plots: 0 don't, 1 do, -1 only stacks
bool debug_        = false  ; //print debug info
bool do_single_    = false  ; //test printing a single histogram
TString only_mode_ = ""     ; //"fit_s"; //fit version to print ("" to ignore)
bool is_prelim_    = false  ;
TString file_type_ = "pdf"  ;
bool data_minus_b_ = true   ; // Data / bkg vs. Data - bkg subplots
bool is_mumem_     = true   ;

//------------------------------------------------------------------------------------------
// Helper functions
TGraphAsymmErrors* make_data_over_model(const char* name, TGraphAsymmErrors* data, TH1* model) {
  TGraphAsymmErrors* ratio = (TGraphAsymmErrors*) data->Clone(name);
  const int nbins = data->GetN();
  for(int bin = 0; bin < nbins; ++bin) {
    double x, y;
    data->GetPoint(bin, x, y);
    const double err_high = data->GetErrorYhigh(bin);
    const double err_low  = data->GetErrorYlow (bin);
    const double model_v  = model->GetBinContent(bin+1);
    const double model_e  = model->GetBinError  (bin+1);
    if(model_v > 0.) {
      ratio->SetPoint      (bin, x, y / model_v);
      ratio->SetPointEYhigh(bin, err_high / model_v);
      ratio->SetPointEYlow (bin, err_low  / model_v);
    } else {
      ratio->SetPoint      (bin, x, 0.);
      ratio->SetPointEYhigh(bin, 0.);
      ratio->SetPointEYlow (bin, 0.);
    }
  }
  return ratio;
}

TGraphAsymmErrors* make_data_minus_model(const char* name, TGraphAsymmErrors* data, TH1* model, int err_mode = 1) {
  TGraphAsymmErrors* diff = (TGraphAsymmErrors*) data->Clone(name);
  const int nbins = data->GetN();
  for(int bin = 0; bin < nbins; ++bin) {
    double x, y;
    data->GetPoint(bin, x, y);
    const double err_high = data->GetErrorYhigh(bin);
    const double err_low  = data->GetErrorYlow (bin);
    const double model_v  = model->GetBinContent(bin+1);
    const double model_e  = model->GetBinError  (bin+1);
    const double val      = y - model_v;
    const double data_err = (y > model_v) ? err_low : err_high;
    diff->SetPoint      (bin, x, val);
    diff->SetPointEYhigh(bin, err_high);
    diff->SetPointEYlow (bin, err_low);
  }
  return diff;
}

TH1* make_data_pull(const char* name, TGraphAsymmErrors* data, TH1* model, int err_mode = 1) {
  TH1* pull = (TH1*) model->Clone(name);
  pull->Reset();
  const int nbins = data->GetN();
  for(int bin = 0; bin < nbins; ++bin) {
    double x, y;
    data->GetPoint(bin, x, y);
    const double err_high = data->GetErrorYhigh(bin);
    const double err_low  = data->GetErrorYlow (bin);
    const double model_v  = model->GetBinContent(bin+1);
    const double model_e  = model->GetBinError  (bin+1);
    const double data_err = (y > model_v) ? err_low : err_high;
    const double err      = (err_mode == 0) ? sqrt(model_e*model_e + data_err*data_err) :
      (model_e > data_err) ? 0. : sqrt(data_err*data_err - model_e*model_e);
    if(err > 0.) {
      pull->SetBinContent(bin+1, (y - model_v) / err);
    } else {
      pull->SetBinContent(bin+1, 0.);
    }
  }
  return pull;
}

TH1* make_model_uncertainty(const char* name, TH1* model) {
  TH1* unc = (TH1*) model->Clone(name);
  for(int ibin = 1; ibin <= model->GetNbinsX(); ++ibin) {
    if(model->GetBinContent(ibin) <= 0.01) unc->SetBinContent(ibin, 0.);
    else {
      unc->SetBinContent(ibin, 1.);
      if(model->GetBinError(ibin) > 0.) {
        unc->SetBinError(ibin, model->GetBinError(ibin) / model->GetBinContent(ibin));
      } else {
        unc->SetBinError(ibin, 0.01);
      }
    }
  }
  return unc;
}

double hmax(TH1* h) {
  double max_val = h->GetBinContent(1);
  for(int ibin = 2; ibin <= h->GetNbinsX(); ++ibin) max_val = max(max_val, h->GetBinContent(ibin));
  return max_val;
}

double hmin(TH1* h, double cutoff = 0.01) {
  double min_val = h->GetMaximum();
  for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
    if(h->GetBinContent(ibin) < cutoff) continue;
    min_val = min(min_val, h->GetBinContent(ibin));
  }
  return max(cutoff, min_val);
}

double gmax(TGraph* g) {
  const int nbins = g->GetN();
  double max_val = -1.;
  for(int ibin = 0; ibin < nbins; ++ibin) {
    const double val = g->GetY()[ibin];
    max_val = (max_val < 0.) ? val : max(max_val, val);
  }
  return max_val;
}

double gmin(TGraph* g, double cutoff = 0.01) {
  const int nbins = g->GetN();
  double min_val = -1;
  bool first = true;
  for(int ibin = 0; ibin < nbins; ++ibin) {
    const double val = g->GetY()[ibin];
    if(val < cutoff) continue;
    min_val = (first) ? val : min(min_val, val);
    first = false;
  }
  return max(cutoff, min_val);
}

void draw_cms_label(double left_margin = 0.10) {
  //CMS prelim drawing
  TText cmslabel;
  cmslabel.SetNDC();
  cmslabel.SetTextColor(1);
  cmslabel.SetTextSize(0.11);
  cmslabel.SetTextAlign(11);
  cmslabel.SetTextAngle(0);
  cmslabel.SetTextFont(61);
  cmslabel.DrawText(left_margin, 0.915, "CMS");
  if(is_prelim_) {
    cmslabel.SetTextFont(52);
    cmslabel.SetTextSize(0.76*cmslabel.GetTextSize());
    cmslabel.DrawText(left_margin + 0.09, 0.915, "Preliminary");
  }
}

void draw_luminosity() {
  return;
}

void draw_category(TString tag, float left_margin = 0.10) {
  return;
}

void replace_bin_labels(TH1* h, const int neff_bins) {
  if(!h) return;
  // // change bin numbers to text labels
  // const int nbins = h->GetNbinsX();
  // for(int bin = 1; bin <= nbins; ++bin) {
  //   if(neff_bins < 10 || (bin-1) % 2 == 0) h->GetXaxis()->SetBinLabel(bin, Form("%i", bin-1));
  //   else                                   h->GetXaxis()->SetBinLabel(bin, Form(" "));
  // }
}

//------------------------------------------------------------------------------------------
// Get a distribution from the directory list
TH1* get_hist(vector<TDirectoryFile*> dirs, const char* name) {
  TH1* h = nullptr;
  for(auto dir : dirs) {
    TH1* h_tmp = (TH1*) dir->Get(name);
    if(!h_tmp) {
      // cout << __func__ << ": Histogram " << name << " not found in directory " << dir->GetName() << endl;
      return nullptr;
    }
    if(!h) {
      h = (TH1*) h_tmp->Clone(Form("%s_Run2", name));
    } else {
      h->Add(h_tmp);
    }
  }
  // Correct the bin width scaling
  for(int bin = 1; bin <= h->GetNbinsX(); ++bin)
    h->SetBinContent(bin, h->GetBinContent(bin)*h->GetBinWidth(bin));
  return h;
}

//------------------------------------------------------------------------------------------
// Get the data distribution from the directory list
TGraphAsymmErrors* get_data(vector<TDirectoryFile*> dirs) {
  TGraphAsymmErrors* g = nullptr;
  for(auto dir : dirs) {
    TGraphAsymmErrors* g_tmp = (TGraphAsymmErrors*) dir->Get("data");
    if(!g_tmp) {
      cout << __func__ << ": Data not found in directory " << dir->GetName() << endl;
      return nullptr;
    }
    if(!g) {
      g = (TGraphAsymmErrors*) g_tmp->Clone("data_Run2");
    } else {
      //Add the data and errors
      for(int ipoint = 0; ipoint < g->GetN(); ++ipoint) {
        g->SetPointY(ipoint, g->GetPointY(ipoint)+g_tmp->GetPointY(ipoint));
        //use sqrt(N) for the errors
        g->SetPointEYhigh(ipoint, sqrt(g->GetPointY(ipoint)));
        g->SetPointEYlow (ipoint, sqrt(g->GetPointY(ipoint)));
      }
    }
  }
  if(g) {
    // Correct the bin width scaling
    for(int ipoint = 0; ipoint < g->GetN(); ++ipoint) {
      const double width = g->GetErrorXhigh(ipoint) + g->GetErrorXlow(ipoint);
      g->SetPointY     (ipoint, g->GetPointY    (ipoint) * width);
      g->SetPointEYhigh(ipoint, g->GetErrorYhigh(ipoint) * width);
      g->SetPointEYlow (ipoint, g->GetErrorYlow (ipoint) * width);
    }

    //no x errors
    for(int ipoint = 0; ipoint < g->GetN(); ++ipoint) {
      g->SetPointEXhigh(ipoint, 0.);
      g->SetPointEXlow (ipoint, 0.);
    }

  }
  return g;
}

//------------------------------------------------------------------------------------------
// Print an individual stack plot
int print_stack(vector<TDirectoryFile*> dirs, TString tag, TString outdir) {

  //Get the fit results and the data
  TH1* hbackground         = get_hist(dirs, "total_background");
  TH1* hsignal             = get_hist(dirs, "total_signal");
  TH1* htotal              = get_hist(dirs, "total");
  TGraphAsymmErrors* gdata = get_data(dirs);
  if(!hsignal || !hbackground || !htotal || !gdata) {
    cout << "Data not found for tag " << tag.Data() << endl;
    return 1;
  }

  //Build the stack
  THStack* stack = new THStack("hstack", ("stack_" + tag).Data());
  vector<TString> names ;
  vector<TString> titles;
  vector<int>     colors;
  if(is_mumem_) {
    names  = {"pbar"      , "rpc_ext"     , "rpc_int"     , "cosmic"    , "dio"};
    titles = {"Antiproton", "External RPC", "Internal RPC", "Cosmic ray", "DIO"};
    colors = {kGreen-6    , kMagenta-10   , kMagenta+1    , kOrange     , kRed-7};
  } else {
    names  = {"pbar"      , "rpc_ext"     , "rpc_int"     , "cosmic"    , "rmc_ext"       , "rmc_int"       };
    titles = {"Antiproton", "External RPC", "Internal RPC", "Cosmic ray", "RMC (external)", "RMC (internal)"};
    colors = {kGreen-6    , kMagenta-10   , kMagenta+1    , kOrange     , kAtlantic+2     , kAtlantic       };
  }
  for(unsigned i = 0; i < names.size(); ++i) {
    TString name = names[i];
    auto h = get_hist(dirs, name.Data());
    if(!h) { //add zero rate processes
      h = (TH1*) hbackground->Clone((names[i]+"_"+tag).Data());
      h->Reset();
    }
    if(h) {
      const int color = colors[i];
      h->SetLineColor(kBlack);
      h->SetLineWidth(1);
      h->SetFillColor(color);
      h->SetTitle(titles[i].Data());
      stack->Add(h);
    }
  }

  //N(data) points, ensure it matches the background model
  const int nbins = gdata->GetN();
  if(nbins != hbackground->GetNbinsX()) {
    cout << "Data and background have different bin numbers!\n";
    return 2;
  }

  //Determine the x-axis range to use
  const double xmin = htotal->GetBinLowEdge(htotal->FindFirstBinAbove(0.01));
  const double xmax = htotal->GetXaxis()->GetBinUpEdge(htotal->FindLastBinAbove(0.01));


  //Create the canvas to plot on
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TCanvas* c = new TCanvas("c_stack", "c_stack", 1000, 900);
  const float x1(0.23), x2(0.43);
  TPad* pad1 = new TPad("pad1_stack", "pad1_stack", 0., x2, 1., 1.);
  TPad* pad2 = new TPad("pad2_stack", "pad2_stack", 0., x1, 1., x2);
  TPad* pad3 = new TPad("pad3_stack", "pad3_stack", 0., 0., 1., x1);
  pad1->SetRightMargin (0.03); pad2->SetRightMargin (pad1->GetRightMargin()); pad3->SetRightMargin(pad1->GetRightMargin());
  pad1->SetLeftMargin  (0.13); pad2->SetLeftMargin  (pad1->GetLeftMargin ()); pad3->SetLeftMargin (pad1->GetLeftMargin ());
  pad1->SetBottomMargin(0.03); pad2->SetBottomMargin(0.09); pad3->SetBottomMargin(0.33);
  pad1->SetTopMargin   (0.10); pad2->SetTopMargin   (0.08); pad3->SetTopMargin   (0.05);
  pad1->Draw(); pad2->Draw(); pad3->Draw();

  // Draw the data and fit components
  pad1->cd();

  //Configure the data style
  gdata->SetMarkerStyle(20);
  gdata->SetMarkerSize(1.2);
  gdata->SetMarkerColor(kBlack);
  gdata->SetLineColor(kBlack);
  gdata->SetLineWidth(2);

  //Configure the total fit (S+B) style
  htotal->SetLineColor(kBlack);
  htotal->SetMarkerColor(kRed);
  htotal->SetFillColor(kRed);
  htotal->SetFillStyle(3003);
  // htotal->SetMarkerStyle(20);
  htotal->SetMarkerSize(0.);
  htotal->SetLineWidth(1);
  htotal->SetLineColor(kRed);
  htotal->SetTitle("");
  htotal->SetXTitle("");
  htotal->SetYTitle("Events / Bin");
  htotal->GetYaxis()->SetTitleSize(0.08);
  htotal->GetYaxis()->SetTitleOffset(0.70);
  htotal->GetXaxis()->SetLabelSize(0.);
  htotal->GetYaxis()->SetLabelSize(0.065);

  //Configure the background component style
  hbackground->SetLineColor(kRed);
  hbackground->SetMarkerColor(kRed);
  hbackground->SetLineWidth(3);
  hbackground->SetLineStyle(kDashed);

  //Configure the signal component style
  hsignal->SetLineColor(kBlue);
  hsignal->SetMarkerColor(kBlue);
  hsignal->SetLineWidth(3);

  //Draw the results
  htotal->Draw("L");
  htotal->SetLineWidth(0);
  stack->Draw("hist noclear same");
  // htotal->Draw("E1 same");
  // htotal->Draw("E2 same");
  if(unblind_) {
    // hbackground->Draw("hist same");
    hsignal->Draw("hist same");
  }
  gdata->Draw("PZ");
  htotal->GetXaxis()->SetRangeUser(xmin, xmax);
  const int neff_bins = htotal->FindBin(xmax-1.e-3) - htotal->FindBin(xmin+1.e-3) + 1;
  const int ndivisions = (neff_bins > 10) ? 200 + (neff_bins+1) / 2 : neff_bins;
  // htotal->GetXaxis()->SetNdivisions(ndivisions);


  //Add a legend for the summary components and one for the background stack
  TLegend leg_sum(pad1->GetLeftMargin() + 0.03, 0.63, pad1->GetLeftMargin() + 0.29, 0.88);
  leg_sum.AddEntry(gdata, "Data", "PE");
  // leg_sum.AddEntry(htotal, "Total", "F");
  if(unblind_) {
    // leg_sum.AddEntry(hbackground, "Background", "L");
    leg_sum.AddEntry(hsignal, "Signal", "L");
  }
  TLegend leg_bkg(pad1->GetLeftMargin() + 0.30, 0.63, 1. - 0.04 - pad1->GetRightMargin(), 0.88);
  leg_bkg.SetNColumns(2);
  for(auto h : *(stack->GetHists())) leg_bkg.AddEntry(h, h->GetTitle(), "F");

  leg_sum.SetTextSize(0.06);
  leg_sum.SetFillStyle(0);
  leg_sum.SetFillColor(0);
  leg_sum.SetLineColor(0);
  leg_sum.SetLineStyle(0);
  leg_sum.Draw();
  leg_bkg.SetTextSize(0.06);
  leg_bkg.SetFillStyle(0);
  leg_bkg.SetFillColor(0);
  leg_bkg.SetLineColor(0);
  leg_bkg.SetLineStyle(0);
  leg_bkg.Draw();

  pad1->Draw(); //redraw the pad
  gPad->RedrawAxis();

  pad2->cd();

  TGraphAsymmErrors* gRatio_s = data_minus_b_ ? make_data_minus_model("gRatio_s_stack", gdata, htotal, err_mode_) : make_data_over_model("gRatio_s_stack", gdata, htotal);
  TGraphAsymmErrors* gRatio_b = data_minus_b_ ? make_data_minus_model("gRatio_b_stack", gdata, hbackground, err_mode_) : make_data_over_model("gRatio_b_stack", gdata, hbackground);
  gRatio_s->SetMarkerSize(1.5);
  gRatio_b->SetMarkerSize(1.5);
  TH1* hPull_s = make_data_pull("hPull_s_stack", gdata, htotal, err_mode_);
  TH1* hPull_b = make_data_pull("hPull_b_stack", gdata, hbackground, err_mode_);
  //Debug printout if needed
  if(debug_) {
    printf("                           Bin  :     data   sqrt(n)         B +-  sigma_B       S+B +- sigma_S+B  pull_B  pull_S+B\n");

  }
  for(int bin = 0; bin < nbins; ++bin) {
    //Retrieve the data point and corresponding model value
    double x, y;
    gdata->GetPoint(bin, x, y);
    const double err_high = gdata->GetErrorYhigh(bin);
    const double err_low  = gdata->GetErrorYlow (bin);
    const double tot_v    = htotal->GetBinContent(bin+1);
    const double tot_e    = htotal->GetBinError  (bin+1);
    const double bkg_v    = hbackground->GetBinContent(bin+1);
    const double bkg_e    = hbackground->GetBinError  (bin+1);

    if(debug_) {
      printf("%25s bin %2i: %8.0f (%8.2f) %8.1f +- %8.2f  %8.1f +- %8.2f   %5.2f    %5.2f\n",
             tag.Data(), bin+1, y, sqrt(y), bkg_v, bkg_e, tot_v, tot_e, hPull_b->GetBinContent(bin+1), hPull_s->GetBinContent(bin+1));
    }
    if(std::fabs(hPull_s->GetBinContent(bin+1)) > 3.f) {
      printf(">>> High pull!\n");
      printf("                           Bin  :     data   sqrt(n)       S+B +- sigma_S+B  Delta +- err     pull_S+B\n");
      const double err_high_calc = (y > tot_v) ? err_low : err_high;
      const double err_s_calc = (err_mode_ == 0) ? sqrt(tot_e*tot_e + err_high_calc*err_high_calc) : (tot_e > err_high_calc) ? 0. : sqrt(err_high_calc*err_high_calc - tot_e*tot_e);
      printf("%25s bin %2i: %8.0f (%8.2f) %8.1f +- %8.2f  %6.1f +- %6.2f   %5.2f\n",
             tag.Data(), bin+1, y, sqrt(y), tot_v, tot_e, y-tot_v, err_s_calc, hPull_s->GetBinContent(bin+1));
    }
  }

  auto gRatio = (unblind_) ? gRatio_s : gRatio_b;

  TH1* hBkg_unc = make_model_uncertainty("hBkg_unc_stack", hbackground);

  TH1* hRatio_s = (TH1*) htotal->Clone("hRatio_s_stack");
  for(int ibin = 1; ibin <= hsignal->GetNbinsX(); ++ibin) {
    const double b = hbackground->GetBinContent(ibin);
    const double t = htotal     ->GetBinContent(ibin);
    if(!data_minus_b_ && b <= 1.e-5) hRatio_s->SetBinContent(ibin, 0.);
    else                             hRatio_s->SetBinContent(ibin, (data_minus_b_) ? t - b : t / b);
  }
  hRatio_s->SetFillColor(0);
  hRatio_s->SetLineColor(kBlue);
  hRatio_s->SetLineStyle(kDashed);
  hRatio_s->SetLineWidth(3);

  const float max_model_diff = (unblind_) ? gmax(gRatio_b) : gmax(gRatio_s);
  const float min_model_diff = (unblind_) ? gmin(gRatio_b,-1.e10) : gmin(gRatio_s,-1.e10);
  const float range = max_model_diff - min_model_diff;
  const float ratio_min = data_minus_b_ ? min_model_diff - 0.2*range : 0.f;
  const float ratio_max = data_minus_b_ ? max_model_diff + 0.2*range : 2.f;

  hBkg_unc->Draw("hist"); //"E2");
  hBkg_unc->SetLineWidth(0);
  if(unblind_) {
    hRatio_s->Draw("hist same");
    hRatio_s->GetYaxis()->SetRangeUser(ratio_min, ratio_max);
  }
  gRatio->SetLineColor(gdata->GetLineColor());
  gRatio->SetLineWidth(gdata->GetLineWidth());
  gRatio->SetMarkerColor(gdata->GetMarkerColor());
  gRatio->SetMarkerStyle(gdata->GetMarkerStyle());
  gRatio->SetMarkerSize(gdata->GetMarkerSize());
  gRatio->Draw("PEZ");
  hBkg_unc->GetXaxis()->SetRangeUser(xmin, xmax);

  TLine* line = new TLine(xmin, data_minus_b_ ? 0. : 1., xmax, data_minus_b_ ? 0. : 1.);
  line->SetLineColor(kBlack);
  line->SetLineWidth(2);
  line->SetLineStyle(kDashed);
  line->Draw("same");

  float txt_scale = (1.-x2)/(x2-x1);
  hBkg_unc->GetYaxis()->SetRangeUser(ratio_min, ratio_max);
  hBkg_unc->SetTitle("");
  hBkg_unc->SetXTitle("");
  hBkg_unc->GetXaxis()->SetLabelSize(0.);
  if(!data_minus_b_) hBkg_unc->GetYaxis()->SetNdivisions(205);
  else               hBkg_unc->GetYaxis()->SetNdivisions(505);
  hBkg_unc->GetYaxis()->SetLabelSize(txt_scale*htotal->GetYaxis()->GetLabelSize());
  hBkg_unc->GetYaxis()->SetTitleSize(txt_scale*htotal->GetYaxis()->GetTitleSize());
  hBkg_unc->GetYaxis()->SetTitleOffset(0.23);
  if(unblind_) hBkg_unc->SetYTitle(data_minus_b_ ? "Data-Bkg" : "Data/Bkg");
  else         hBkg_unc->SetYTitle(data_minus_b_ ? "Data-Fit" : "Data/Fit");
  hBkg_unc->GetYaxis()->CenterTitle(true);


  //Make a pull plot
  pad3->cd();

  pad3->cd();

  auto hPull = (TH1*) hPull_s;
  replace_bin_labels(hPull, neff_bins);
  const double max_pull = hPull->GetMaximum();
  const double min_pull = hPull->GetMinimum();
  hPull->SetLineColor(kAtlantic);
  hPull->SetFillColor(kAtlantic);
  hPull->SetFillStyle(1000);
  hPull->Draw("hist");
  hPull->GetYaxis()->SetRangeUser(min(-3., 1.1*min_pull),max(3., 1.1*max_pull));
  hPull->SetTitle("");
  hPull->SetXTitle("Momentum (MeV/c)");
  hPull->SetYTitle("Pull");

  txt_scale = (1.-x2)/(x1);
  hPull->GetXaxis()->SetLabelOffset(0.012);
  hPull->GetXaxis()->SetLabelSize(1.*txt_scale*htotal->GetYaxis()->GetLabelSize());
  hPull->GetYaxis()->SetNdivisions(505);
  hPull->GetYaxis()->SetLabelSize(txt_scale*htotal->GetYaxis()->GetLabelSize());
  hPull->GetXaxis()->SetTitleSize(txt_scale*htotal->GetYaxis()->GetTitleSize());
  hPull->GetYaxis()->SetTitleSize(hPull->GetXaxis()->GetTitleSize());
  hPull->GetXaxis()->SetTitleOffset(0.76);
  hPull->GetYaxis()->SetTitleOffset(0.27);
  hPull->GetYaxis()->CenterTitle(true);

  //Add a reference line for perfect agreement
  TLine* line_2 = new TLine(xmin, 0., xmax, 0.);
  line_2->SetLineColor(kBlack);
  line_2->SetLineWidth(2);
  line_2->SetLineStyle(kDashed);
  line_2->Draw("same");

  //Add the CMS label
  pad1->cd();
  // draw_cms_label(pad1->GetLeftMargin());
  // draw_luminosity();
  // draw_category(outdir + tag, pad1->GetLeftMargin());

  //Print a linear and a log version of the distribution
  double min_val = std::max(0.1, std::min(gmin(gdata), hmin(htotal)));
  double max_val = std::max(gdata->GetMaximum(), hmax(htotal));
  htotal->GetYaxis()->SetRangeUser(0., 1.75*max_val);
  c->SaveAs(Form("%s%s_stack.%s", outdir.Data(), tag.Data(), file_type_.Data()));
  c->SaveAs(Form("%s%s_stack.root", outdir.Data(), tag.Data()));
  double plot_min = 1.e-4*min_val;
  double plot_max = plot_min*pow(10, 1.75*log10(max_val/plot_min));
  htotal->GetYaxis()->SetRangeUser(plot_min, plot_max);
  pad1->SetLogy();
  c->SaveAs(Form("%s%s_stack_logy.%s", outdir.Data(), tag.Data(), file_type_.Data()));
  c->SaveAs(Form("%s%s_stack_logy.root", outdir.Data(), tag.Data()));

  double ndata = 0.;
  for (int i = 0; i < gdata->GetN(); ++i) ndata += gdata->GetPointY(i);
  printf("%s: N(data) = %.1f\n", __func__, ndata);

  //Clean up after printing
  delete c;
  delete gRatio_s;
  delete gRatio_b;
  delete hPull_s;
  delete hPull_b;
  delete hBkg_unc;
  delete hRatio_s;
  delete stack;
  delete line;
  delete line_2;

  return 0;
}

//------------------------------------------------------------------------------------------
// Print an individual histogram
int print_hist(vector<TDirectoryFile*> dirs, TString tag, TString outdir) {

  //Get the fit results and the data
  TH1* hbackground         = get_hist(dirs, "total_background");
  TH1* hsignal             = get_hist(dirs, "total_signal");
  TH1* htotal              = get_hist(dirs, "total");
  TGraphAsymmErrors* gdata = get_data(dirs);
  if(!hsignal || !hbackground || !htotal || !gdata) {
    cout << "Data not found for tag " << tag.Data() << endl;
    return 1;
  }

  //N(data) points, ensure it matches the background model
  const int nbins = gdata->GetN();
  if(nbins != hbackground->GetNbinsX()) {
    cout << "Data and background have different bin numbers!\n";
    return 2;
  }

  //Determine the x-axis range to use
  const double xmin = htotal->GetBinLowEdge(htotal->FindFirstBinAbove(0.01));
  const double xmax = htotal->GetXaxis()->GetBinUpEdge(htotal->FindLastBinAbove(0.01));

  //Create the canvas to plot on
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TCanvas* c = new TCanvas("c", "c", 1000, 900);
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0.40, 1., 1.00);
  TPad* pad2 = new TPad("pad2", "pad2", 0., 0.20, 1., 0.40);
  TPad* pad3 = new TPad("pad3", "pad3", 0., 0.00, 1., 0.20);
  pad1->SetRightMargin(0.03); pad2->SetRightMargin(0.03); pad3->SetRightMargin(0.03);
  pad1->SetBottomMargin(0.02); pad2->SetBottomMargin(0.05); pad3->SetBottomMargin(0.25);
  pad2->SetTopMargin(0.03); pad3->SetTopMargin(0.05);
  pad1->Draw(); pad2->Draw(); pad3->Draw();

  // Draw the data and fit components
  pad1->cd();

  //Configure the data style
  gdata->SetMarkerStyle(20);
  gdata->SetMarkerSize(0.8);
  gdata->SetLineWidth(3);

  //Configure the total fit (S+B) style
  htotal->SetLineColor(kRed);
  htotal->SetMarkerColor(kRed);
  htotal->SetFillColor(kRed);
  htotal->SetFillStyle(3003);
  // htotal->SetMarkerStyle(20);
  htotal->SetMarkerSize(0.);
  htotal->SetLineWidth(3);
  htotal->SetLineColor(kRed);
  htotal->SetTitle("");
  htotal->SetXTitle("");
  htotal->SetYTitle("Events / Bin");
  htotal->GetYaxis()->SetTitleSize(0.05);
  htotal->GetYaxis()->SetTitleOffset(0.92);
  htotal->GetXaxis()->SetLabelSize(0.);

  //Configure the background component style
  hbackground->SetLineColor(kRed);
  hbackground->SetMarkerColor(kRed);
  hbackground->SetLineWidth(3);
  hbackground->SetLineStyle(kDashed);

  //Configure the signal component style
  hsignal->SetLineColor(kBlue);
  hsignal->SetMarkerColor(kBlue);
  hsignal->SetLineWidth(3);

  //Draw the results
  htotal->Draw("E1");
  htotal->Draw("E2 same");
  if(unblind_) {
    hbackground->Draw("hist same");
    hsignal->Draw("hist same");
  }
  gdata->Draw("P");
  htotal->GetXaxis()->SetRangeUser(xmin, xmax);


  //Add a legend
  TLegend leg(0.6, 0.5, 0.85, 0.85);
  leg.AddEntry(gdata, "Data", "PLE");
  leg.AddEntry(htotal, "Background+signal fit", "LF");
  if(unblind_) {
    leg.AddEntry(hbackground, "Background component", "L");
    leg.AddEntry(hsignal, "Signal component", "L");
  }
  leg.SetTextSize(0.05);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetLineColor(0);
  leg.SetLineStyle(0);
  leg.Draw();

  //Make the ratio plots
  pad2->cd();

  //Make the data / total fit and data / background component distributions
  TGraphAsymmErrors* gRatio_s = data_minus_b_ ? make_data_minus_model("gRatio_s", gdata, htotal, err_mode_) : make_data_over_model("gRatio_s", gdata, htotal);
  TGraphAsymmErrors* gRatio_b = data_minus_b_ ? make_data_minus_model("gRatio_b", gdata, hbackground, err_mode_) : make_data_over_model("gRatio_b", gdata, hbackground);
  TH1* hPull_s = make_data_pull("hPull_s", gdata, htotal, err_mode_);
  TH1* hPull_b = make_data_pull("hPull_b", gdata, hbackground, err_mode_);
  //Debug printout if needed
  if(debug_) {
    printf("Bin:     data        B +-  sigma_B       S+B +- sigma_S+B  pull_B  pull_S+B\n");

  }
  for(int bin = 0; bin < nbins; ++bin) {
    //Retrieve the data point and corresponding model value
    double x, y;
    gdata->GetPoint(bin, x, y);
    const double err_high = gdata->GetErrorYhigh(bin);
    const double err_low  = gdata->GetErrorYlow (bin);
    const double tot_v    = htotal->GetBinContent(bin+1);
    const double tot_e    = htotal->GetBinError  (bin+1);
    const double bkg_v    = hbackground->GetBinContent(bin+1);
    const double bkg_e    = hbackground->GetBinError  (bin+1);

    //Debug printout if needed
    if(debug_) {
      printf(" %2i: %8.0f %8.1f +- %8.2f  %8.1f +- %8.2f   %5.2f    %5.2f\n",
             bin+1, y, bkg_v, bkg_e, tot_v, tot_e, hPull_b->GetBinContent(bin+1), hPull_s->GetBinContent(bin+1));
    }
  }

  auto gRatio = (unblind_) ? gRatio_b : gRatio_s;

  TH1* hBkg_unc = make_model_uncertainty("hBkg_unc", hbackground);
  TH1* hRatio_s = (TH1*) htotal->Clone("hRatio_s");
  hRatio_s->Divide(hbackground);
  for(int ibin = 1; ibin <= hsignal->GetNbinsX(); ++ibin) {
    if(hbackground->GetBinContent(ibin) <= 0.1) hRatio_s->SetBinContent(ibin, 0.);
  }
  hRatio_s->SetFillColor(0);
  hRatio_s->SetLineColor(kBlue);
  hRatio_s->SetLineStyle(kDashed);

  const float ratio_min = data_minus_b_ ? -3.f : 0.9f;
  const float ratio_max = data_minus_b_ ?  3.f : 1.1f;

  hBkg_unc->Draw("E2");
  if(unblind_) {
    hRatio_s->Draw("hist same");
    hRatio_s->GetYaxis()->SetRangeUser(ratio_min, ratio_max);
  }
  gRatio->Draw("P");
  hBkg_unc->GetXaxis()->SetRangeUser(xmin, xmax);

  TLine* line = new TLine(xmin, data_minus_b_ ? 0. : 1., xmax, data_minus_b_ ? 0. : 1.);
  line->SetLineColor(kBlack);
  line->SetLineWidth(2);
  line->SetLineStyle(kDashed);
  line->Draw("same");

  hBkg_unc->GetYaxis()->SetRangeUser(ratio_min, ratio_max);
  hBkg_unc->SetTitle("");
  hBkg_unc->SetXTitle("");
  hBkg_unc->GetXaxis()->SetLabelSize(0.);
  hBkg_unc->GetYaxis()->SetLabelSize(0.10);
  hBkg_unc->GetYaxis()->SetTitleSize(0.15);
  hBkg_unc->GetYaxis()->SetTitleOffset(0.30);
  if(unblind_) hBkg_unc->SetYTitle(data_minus_b_ ? "Data-Bkg" : "Data/Bkg");
  else         hBkg_unc->SetYTitle(data_minus_b_ ? "Data-Fit" : "Data/Fit");


  pad3->cd();

  auto hPull = (TH1*) hPull_s;
  const double max_pull = hPull->GetMaximum();
  const double min_pull = hPull->GetMinimum();
  hPull->SetLineColor(kAtlantic);
  hPull->SetFillColor(kAtlantic);
  hPull->SetFillStyle(1000);
  hPull->Draw("hist");
  hPull->GetYaxis()->SetRangeUser(min(-3., 1.1*min_pull),max(3., 1.1*max_pull));
  hPull->SetTitle("");
  hPull->SetXTitle("Momentum (MeV/c)");
  hPull->SetYTitle("Pull");
  hPull->GetXaxis()->SetLabelSize(0.10);
  hPull->GetYaxis()->SetLabelSize(0.10);
  hPull->GetXaxis()->SetTitleSize(0.15);
  hPull->GetYaxis()->SetTitleSize(0.15);
  hPull->GetXaxis()->SetTitleOffset(0.75);
  hPull->GetYaxis()->SetTitleOffset(0.29);

  //Add a reference line for perfect agreement
  TLine* line_2 = new TLine(xmin, 0., xmax, 0.);
  line_2->SetLineColor(kBlack);
  line_2->SetLineWidth(2);
  line_2->SetLineStyle(kDashed);
  line_2->Draw("same");

  //Add the CMS label
  pad1->cd();
  // draw_cms_label(pad1->GetLeftMargin());
  // draw_luminosity(year);

  //Print a linear and a log version of the distribution
  double min_val = std::max(0.1, std::min(gmin(gdata), hmin(htotal)));
  double max_val = std::max(gdata->GetMaximum(), hmax(htotal));
  htotal->GetYaxis()->SetRangeUser(0., 1.2*max_val);
  c->SaveAs(Form("%s%s.%s", outdir.Data(), tag.Data(), file_type_.Data()));
  c->SaveAs(Form("%s%s.root", outdir.Data(), tag.Data()));
  double plot_min = std::min(std::max(0.2, 0.2*hmax(hsignal)), 0.2*min_val);
  double plot_max = plot_min*pow(10, 1.7*log10(max_val/plot_min));
  htotal->GetYaxis()->SetRangeUser(plot_min, plot_max);
  pad1->SetLogy();
  c->SaveAs(Form("%s%s_logy.%s", outdir.Data(), tag.Data(), file_type_.Data()));
  c->SaveAs(Form("%s%s_logy.root", outdir.Data(), tag.Data()));

  double ndata = 0.;
  for (int i = 0; i < gdata->GetN(); ++i) ndata += gdata->GetPointY(i);
  printf("%s: N(data) = %.1f\n", __func__, ndata);

  //Clean up after printing
  delete c;
  delete gRatio_s;
  delete gRatio_b;
  delete hPull_s;
  delete hPull_b;
  delete hBkg_unc;
  delete hRatio_s;
  delete line;
  delete line_2;

  return 0;
}

//------------------------------------------------------------------------------------------
// Print the fit results for each category in a fit configuration directory
int print_dir(TDirectoryFile* dir, TString tag, TString outdir) {
  int status(0);
  if(!dir) return 1;

  //List of categories
  TList* list = dir->GetListOfKeys();
  if(!list) return 10;
  bool subdir(false); //whether there are sub-directories or not
  for(TObject* o : *list) {
    TObject* obj = dir->Get(o->GetName());
    if(!obj) continue;

    //Check if this object is a directory with categories, if so recursively process it
    bool isdir = obj->InheritsFrom(TDirectoryFile::Class());
    if(isdir) {
      auto next_dir = (TDirectoryFile*) obj;
      status += print_dir(next_dir, (tag + "_") + obj->GetName(), outdir);
      if(do_single_) return status;
      subdir = true;
    }
  }

  //If this directory doesn't contain a sub-directory, print the histograms within the category
  if(!subdir && tag.Contains(only_mode_)) { //histogram directory
    if(is_prelim_) tag += "_prelim";
    if(print_stacks_ >= 0) status += print_hist({dir}, tag, outdir);
    if(print_stacks_ != 0) status += print_stack({dir}, tag, outdir); //stacked histogram
  }
  return status;
}

//------------------------------------------------------------------------------------------
// Print all fit figures
int plot_combine_fit(TString fname, TString outdir = "figures", bool unblind = false) {
  unblind_ = unblind;
  is_mumem_ = !(fname.Contains("mumep") || outdir.Contains("mumep"));

  //Get the fit file
  TFile* file = TFile::Open(fname.Data(), "READ");
  if(!file) return 1;

  TDirectoryFile *prefit = (TDirectoryFile*) file->Get("shapes_prefit");
  TDirectoryFile *fit_b  = (TDirectoryFile*) file->Get("shapes_fit_b");
  TDirectoryFile *fit_s  = (TDirectoryFile*) file->Get("shapes_fit_s");

  if(!prefit || (unblind && (!fit_b || !fit_s))) {
    cout << "Fit directories not found!\n";
    return 2;
  }

  //Create the figure directory
  if(!outdir.EndsWith("/")) outdir += "/";
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", outdir.Data(), outdir.Data()));

  // Reduce number of entries listed in linear, move exponent
  TGaxis::SetMaxDigits(3);
  TGaxis::SetExponentOffset(-0.06, 0.01, "Y");

  //Print the fit configurations: Pre-fit, background-only fit, and background+signal fit
  int status(0);
  status += print_dir(prefit, "prefit", outdir);
  if(!do_single_) {
    if(fit_b) status += print_dir(fit_b , "fit_b" , outdir);
    if(fit_s) status += print_dir(fit_s , "fit_s" , outdir);
  }

  cout << "Plotting status = " << status << endl;
  return status;
}
