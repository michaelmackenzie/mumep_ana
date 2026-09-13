// Build the statistical model

#include "../tools/types.C"
#include "signal_model.C"
#include "background_model.C"
#include "../physics.C"
#include "systematics.C"
#include "create_envelope.C"
#include "../tools/write_datacard.C"
// #include "combine/HiggsAnalysis/CombinedLimit/src/RooLandauCB.cc"

bool print_      = true;
bool write_card_ = true;

//---------------------------------------------------------------------------------------------------------------------------
struct RateUnc_t {
  TString name;
  double value;
  bool isBeam;
  TString process;
  RateUnc_t(TString name, double value, bool isBeam = false, TString process = "") :
    name(name), value(value), isBeam(isBeam), process(process) {}
};

std::vector<RateUnc_t> rate_uncertainties(TString process) {
  std::vector<RateUnc_t> sys;
  sys.push_back(RateUnc_t("lumi"  , 0.1  , true               ));
  sys.push_back(RateUnc_t("cosmic", 0.2  , false, "cosmic"    ));
  sys.push_back(RateUnc_t("dio"   , 0.025, false, "dio"       ));
  sys.push_back(RateUnc_t("rpc"   , 0.27 , false, "rpc"       ));
  sys.push_back(RateUnc_t("pbar"  , 1.   , false, "pbar"      ));
  sys.push_back(RateUnc_t("rmc"   , 0.079, false, "rmc"       )); // 1.40 +- 0.11 (TRIUMF, 1999)
  sys.push_back(RateUnc_t("rmc"   , 0.045, false, "rmc_int"   )); // 0.0069 ± 0.00031
  return sys;
}

const TH1* preferred_momentum_hist(const pdf_info& info) {
  if(info.smoothed_hist_)   return info.smoothed_hist_;
  if(info.normalized_hist_) return info.normalized_hist_;
  return info.hist_;
}

const TH1* preferred_time_hist(const pdf_info& info) {
  if(info.t0_smoothed_hist_)   return info.t0_smoothed_hist_;
  if(info.t0_normalized_hist_) return info.t0_normalized_hist_;
  return info.t0_raw_hist_;
}

TH2* make_independent_2d_hist(const TH1* momentum_hist,
                              const TH1* time_hist,
                              const TString& name,
                              const TString& title,
                              const double target_rate) {
  if(!momentum_hist || !time_hist) return nullptr;

  const double momentum_integral = momentum_hist->Integral();
  const double time_integral = time_hist->Integral();
  if(momentum_integral <= 0. || time_integral <= 0.) return nullptr;

  const int momentum_bins = momentum_hist->GetNbinsX();
  const int time_bins = time_hist->GetNbinsX();
  auto h2 = new TH2D(name, title,
                     momentum_bins, momentum_hist->GetXaxis()->GetBinLowEdge(1), momentum_hist->GetXaxis()->GetBinUpEdge(momentum_bins),
                     time_bins, time_hist->GetXaxis()->GetBinLowEdge(1), time_hist->GetXaxis()->GetBinUpEdge(time_bins));
  h2->GetXaxis()->SetTitle(momentum_hist->GetXaxis()->GetTitle());
  h2->GetYaxis()->SetTitle(time_hist->GetXaxis()->GetTitle());

  for(int ix = 1; ix <= momentum_bins; ++ix) {
    const double px = momentum_hist->GetBinContent(ix) / momentum_integral;
    for(int iy = 1; iy <= time_bins; ++iy) {
      const double py = time_hist->GetBinContent(iy) / time_integral;
      const double value = target_rate * px * py;
      h2->SetBinContent(ix, iy, value);
      h2->SetBinError(ix, iy, std::sqrt(std::max(0., value)));
    }
  }
  return h2;
}

bool add_independent_2d_inputs(RooWorkspace& ws,
                               const pdf_info& info,
                               const TString& process,
                               const int selection,
                               const TString& component_name,
                               RooRealVar& momentum_obs,
                               RooRealVar& time_obs,
                               const TString& figdir,
                               TFile* out_file = nullptr) {
  const TH1* momentum_hist = preferred_momentum_hist(info);
  const TH1* time_hist = preferred_time_hist(info);
  if(!momentum_hist || !time_hist) {
    cout << __func__ << ": Missing 1D inputs for 2D model component " << component_name.Data() << endl;
    return false;
  }

  const TString base_name = Form("%s_%i_%s", process.Data(), selection, component_name.Data());
  TH2* h2 = make_independent_2d_hist(momentum_hist,
                                     time_hist,
                                     Form("%s_hist", base_name.Data()),
                                     Form("%s 2D histogram", info.title_.Data()),
                                     info.rate_);
  if(!h2) {
    cout << __func__ << ": Unable to build 2D histogram for component " << component_name.Data() << endl;
    return false;
  }

  RooArgList vars;
  vars.add(momentum_obs);
  vars.add(time_obs);
  RooDataHist data(Form("%s_2d_data_hist", base_name.Data()), Form("%s 2D data histogram", info.title_.Data()), vars, h2);
  RooHistPdf pdf(Form("%s_pdf", base_name.Data()), Form("%s 2D PDF", info.title_.Data()), vars, data);

  if(figdir != "") {
    const TString plot_dir = Form("%s/2d", figdir.Data());
    gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", plot_dir.Data(), plot_dir.Data()));

    auto save_2d_plot = [&](const TString& suffix, const bool logz) {
      TCanvas c(Form("c_%s%s", base_name.Data(), suffix.Data()), Form("c_%s%s", base_name.Data(), suffix.Data()), 1200, 1000);
      c.SetRightMargin(0.16);
      c.SetLeftMargin(0.12);
      c.SetBottomMargin(0.12);
      c.SetTopMargin(0.06);
      if(logz) {
        const double min_positive = std::max(1.e-12, h2->GetMinimum(1.e-12));
        h2->SetMinimum(min_positive);
        c.SetLogz();
      } else {
        h2->SetMinimum(0.);
      }
      h2->Draw("COLZ");
      c.SaveAs(Form("%s/%s%s.png", plot_dir.Data(), base_name.Data(), suffix.Data()));
    };

    save_2d_plot("", false);
    save_2d_plot("_logz", true);
  }

  ws.import(data);
  ws.import(pdf);
  if(out_file) {
    out_file->cd();
    h2->Write();
  }
  delete h2;
  return true;
}

//---------------------------------------------------------------------------------------------------------------------------
void print_model(TString figdir, const int selection, RooRealVar& obs, RooAbsData* data,
                 pdf_info& signal_model,
                 std::vector<pdf_info>& background_model, const bool is_mumem) {
  const double signal_scale = (is_mumem) ? 20. : 1.7e3;

  TCanvas* c = new TCanvas("c_model", "c_model", 1200, 1000);
  auto frame = obs.frame();
  frame->SetTitle(Form("%s model", signal_model.title_.Data()));
  frame->SetXTitle("Momentum (MeV/c)");

  // data the data
  gStyle->SetEndErrorSize(0);
  data->plotOn(frame, RooFit::Name("data"), RooFit::XErrorSize(0), RooFit::LineWidth(2));

  // draw the signal
  auto sig_pdf = signal_model.pdf_;
  sig_pdf->plotOn(frame, RooFit::Name(signal_model.name_),
                  RooFit::LineColor(signal_model.color_), RooFit::FillColor(signal_model.color_),
                  RooFit::FillStyle(3005),
                  RooFit::Normalization(signal_scale*signal_model.rate_, RooAbsReal::NumEvent));

  // draw the backgrounds
  RooArgList bkg_pdfs;
  RooArgList bkg_rates;
  for(auto& bkg : background_model) {
    if(!bkg.pdf_) continue;
    bkg_pdfs.add(*bkg.pdf_);
    bkg_rates.add(*(new RooRealVar(Form("%s_rate", bkg.name_.Data()), "", bkg.rate_)));
  }
  RooAddPdf tot_bkg("tot_bkg", "Total background", bkg_pdfs, bkg_rates);
  tot_bkg.plotOn(frame, RooFit::Invisible(), RooFit::Name("bkg"));

  for(auto& bkg : background_model) {
    if(!bkg.pdf_) {
      cout << __func__ << ": Background " << bkg.name_.Data() << " has an undefined PDF!\n";
      continue;
    }
    bkg.pdf_->plotOn(frame, RooFit::Name(bkg.name_),
                     RooFit::LineColor(bkg.color_), RooFit::Normalization(bkg.rate_, RooAbsReal::NumEvent));
    cout << "bkg " << bkg.name_.Data() << " norm " << bkg.rate_ << endl;
  }
  frame->SetYTitle("");
  frame->Draw();

  c = plot_fit_frame(frame, obs, "Momentum (MeV/c)", Form("Events / %.1f MeV/c", bin_width_), "data", "bkg", npot_, livetime_, nmuons_);
  auto pad1 = (TPad*) c->GetPrimitive("pad1");
  gPad->SetTickx(1);
  gPad->SetTicky(1);

  // add a legend
  TLegend* leg = new TLegend((pad1) ? 0.03 + pad1->GetLeftMargin() : 0.13, 0.75, (pad1) ? 0.98 - pad1->GetRightMargin() : 0.88, (pad1) ? 0.97 - pad1->GetTopMargin() : 0.9);
  leg->SetNColumns(3); leg->SetLineWidth(0); leg->SetFillColor(0); leg->SetLineColor(0); leg->SetFillStyle(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.055);
  leg->AddEntry(signal_model.name_, "Signal", "L");
  for(auto& bkg : background_model) leg->AddEntry(bkg.name_, bkg.title_, "L");
  leg->Draw();

  c->SaveAs(Form("%s/input_pdfs_%i.png", figdir.Data(), selection));
  if(pad1) {
    frame->GetYaxis()->SetRangeUser(1.e-6, 100.*max(frame->GetMaximum(), 1.e2));
    pad1->SetLogy();
  }
  c->SaveAs(Form("%s/input_pdfs_%i_log.png", figdir.Data(), selection));
  delete frame;
  delete c;

  // Plot a stack model
  c = standard_canvas("c");
  gStyle->SetOptStat(0);

  THStack* stack = new THStack("stack","");
  std::sort(background_model.begin(), background_model.end(), [](const pdf_info& a, const pdf_info& b) {
        return a.rate_ < b.rate_;
  });

  leg = new TLegend(0.03 + c->GetLeftMargin(), 0.75, 0.98 - c->GetRightMargin(), 0.97 - c->GetTopMargin());
  leg->SetNColumns(3); leg->SetLineWidth(0); leg->SetFillColor(0); leg->SetLineColor(0); leg->SetFillStyle(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.043);

  auto h_sig = signal_model.pdf_->createHistogram("h_sig", obs);
  h_sig->Scale(signal_model.rate_);
  h_sig->SetLineColor(signal_model.color_);
  h_sig->SetLineWidth(2);
  h_sig->SetFillStyle(3004);
  h_sig->SetFillColor(signal_model.color_);
  leg->AddEntry(h_sig, signal_model.title_, "L");

  for(auto& bkg : background_model) {
    TH1* h = bkg.pdf_->createHistogram(bkg.name_, obs);
    h->Scale(bkg.rate_);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->SetFillColor(bkg.color_);
    h->SetFillStyle(kSolid);
    stack->Add(h);
    leg->AddEntry(h, bkg.title_, "F");
  }


  h_sig->GetXaxis()->SetLabelFont(132);
  h_sig->GetXaxis()->SetTitleFont(132);
  h_sig->GetYaxis()->SetLabelFont(132);
  h_sig->GetYaxis()->SetTitleFont(132);
  h_sig->Draw("hist");
  stack->Draw("hist noclear same");
  h_sig->Draw("hist same");
  h_sig->SetTitle(Form(";Momentum (MeV/c);Events / (%.2g MeV/c)", h_sig->GetBinWidth(1)));

  leg->Draw();

  c->SaveAs(Form("%s/input_stack_%i.png", figdir.Data(), selection));
  const double ymin = 1.e-4;
  const double max_val = max(((TH1*) stack->GetStack()->Last())->GetMaximum(), h_sig->GetMaximum());
  const double ymax = ymin*std::pow(max_val/ymin, 1./0.73);
  h_sig->GetYaxis()->SetRangeUser(ymin, ymax);
  c->SetLogy();
  c->SaveAs(Form("%s/input_stack_%i_log.png", figdir.Data(), selection));
  delete c;
  delete h_sig;
  delete stack;
}

//---------------------------------------------------------------------------------------------------------------------------
int build_model(TString process = "mumem", int selection = 20, TString tag = "") {
  if(use_evtana_) set_evtana_defaults();
  init_physics(tag);
  process.ToLower();

  // Create the observable
  const bool is_mumem = process == "mumem";
  const float xmin(is_mumem ? xmin_em_ : xmin_ep_), xmax(is_mumem ? xmax_em_ : xmax_ep_);
  RooRealVar obs(Form("obs_%i", selection), "p", (xmin+xmax)/2., xmin, xmax, "MeV/c");
  const int nbins = (bin_width_ > 0.) ? (xmax - xmin)/bin_width_ + 1.e-6 : 100;
  // cout << "----- nbins = " << nbins << endl;
  // return 1;
  obs.setBins(nbins);

  // Number of signal events expected to be generated
  const double n_signal_exp = nmuons_ * muon_capture_fraction_ * signal_br_;
  TString figdir = Form("figures/%s%s", process.Data(), (tag != "") ? ("_"+tag).Data() : "");

  // Retrieve the signal data
  auto signal_model     = read_model          ("signal", process, selection, tag);
  auto background_model = get_background_model(obs     , process, selection, tag);
  RooAbsData* data      = get_data            (obs     , process, selection, tag);

  if(do_2d_fit_ && !hist_pdfs_) {
    cout << __func__ << ": 2D fits require hist_pdfs_ = true for now." << endl;
    return 1;
  }

  auto sig_pdf = signal_model.pdf_;

  if(!sig_pdf || background_model.empty()) {
    cout << "Model PDFs not found!\n";
    return 1;
  }
  for(auto& bkg : background_model) {
    if(!bkg.pdf_) {
      cout << "Background " << bkg.name_.Data() << " PDF not found!\n";
      return 1;
    }
  }

  // Generate toy data
  if(!data) {
    if(unbinned_) {
      auto generated_data = static_cast<RooDataSet*>(nullptr);
      for(auto& bkg : background_model) {
        auto gen_data = bkg.pdf_->generate(obs, bkg.rate_);
        if(!gen_data) {
          cout << __func__ << ": Gen data for process " << bkg.name_.Data() << " is null!\n";
        } else if(!generated_data) {
          generated_data = gen_data;
          generated_data->SetName("data_obs");
        } else {
          generated_data->append(*gen_data);
        }
      }
      data = generated_data;
    } else {
      auto generated_data = static_cast<RooDataHist*>(nullptr);
      for(auto& bkg : background_model) {
        auto gen_data = bkg.pdf_->generateBinned(obs, bkg.rate_);
        if(!gen_data) {
          cout << __func__ << ": Gen data for process " << bkg.name_.Data() << " is null!\n";
        } else if(!generated_data) {
          generated_data = gen_data;
          generated_data->SetName("data_obs");
        } else {
          generated_data->add(*gen_data);
        }
      }
      data = generated_data;
    }
  }

  // Build the envelope for mu- --> e+, if requested
  // std::vector<TString> explicit_processes = {"cosmic"}; // Processes to model outside of the envelope
  std::vector<TString> explicit_processes = {}; // Processes to model outside of the envelope
  double p_blind_min = 90.;
  double p_blind_max = 93.;
  envelope_t envelope;
  if(use_env_ && process == "mumep") {
    auto data_hist = dynamic_cast<RooDataHist*>(data);
    if(!data_hist) {
      cout << __func__ << ": The background envelope requires binned data (unbinned_ = false)!\n";
      return 1;
    }

    // Split the background model into the processes modeled explicitly and the ones the
    // envelope absorbs into its data-driven shape
    auto is_explicit = [&](const TString& name) {
      for(auto& proc : explicit_processes) if(proc == name) return true;
      return false;
    };
    std::vector<pdf_info> explicit_model;
    for(auto& bkg : background_model) {
      if(is_explicit(bkg.name_)) explicit_model.push_back(bkg);
    }
    if(explicit_model.size() != explicit_processes.size()) {
      cout << __func__ << ": Only found " << explicit_model.size() << " / " << explicit_processes.size()
           << " explicitly modeled processes in the background model!\n";
      return 1;
    }

    envelope = build_background_envelope(obs, *data_hist, explicit_model,
                                         Form("%s_%i_env", process.Data(), selection),
                                         p_blind_min, p_blind_max, verbose_, figdir);
    if(!envelope.pdf_ || !envelope.norm_) {
      cout << __func__ << ": Failed to build the background envelope!\n";
      return 1;
    }

    // The envelope stands in for every background it absorbed, alongside the explicit ones
    pdf_info env_info;
    env_info.pdf_   = envelope.pdf_;
    env_info.norm_  = envelope.norm_;
    env_info.rate_  = envelope.rate_;
    env_info.name_  = "env";
    env_info.title_ = envelope.title_;
    env_info.color_ = kGray;
    env_info.hist_  = envelope.pdf_->createHistogram(Form("%s_%i_env_hist", process.Data(), selection), obs);
    background_model = explicit_model;
    background_model.push_back(env_info);
  }

  // Draw the inputs
  if(print_) {
    gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir.Data(), figdir.Data()));
    print_model(figdir, selection, obs, data, signal_model, background_model, process == "mumem");
  }

  // Open the output file
  gSystem->Exec("[ ! -d workspaces ] && mkdir workspaces");
  TString ws_file = Form("workspaces/workspace_%s_%i%s.root", process.Data(), selection, (tag != "") ? ("_"+tag).Data() : "");
  TString comp_file = ws_file; // for comparisons with other analyses
  comp_file.ReplaceAll("workspace_", "comp_");

  // Construct the comparison file
  TFile* fcomp = new TFile(comp_file, "RECREATE");
  auto comp_dir = fcomp->mkdir("category");
  comp_dir->cd();
  TH1* hsig_comp = (TH1*) signal_model.hist_->Clone("signal");
  hsig_comp->SetTitle(signal_model.title_);
  hsig_comp->Scale(signal_model.rate_ / hsig_comp->Integral());
  hsig_comp->Write();
  for(auto& bkg : background_model) {
    TH1* hbkg_comp = (TH1*) bkg.hist_->Clone(bkg.name_);
    hbkg_comp->Scale(bkg.rate_ / hbkg_comp->Integral());
    hbkg_comp->SetTitle(bkg.title_);
    hbkg_comp->Write();
  }

  auto workflow_dir = fcomp->mkdir("workflow");
  auto raw_dir = workflow_dir->mkdir("raw");
  auto normalized_dir = workflow_dir->mkdir("normalized");
  auto smoothed_dir = workflow_dir->mkdir("smoothed");

  auto write_workflow_hist = [&](TH1* h_in,
                                 const TString& out_name,
                                 const TString& out_title,
                                 const double target_rate,
                                 TDirectory* out_dir,
                                 const bool normalize_to_rate = true) {
    if(!h_in || !out_dir) return;
    TH1* h_out = (TH1*) h_in->Clone(out_name.Data());
    h_out->SetTitle(out_title.Data());
    if(normalize_to_rate) {
      const double integral = h_out->Integral();
      if(integral > 0.) h_out->Scale(target_rate / integral);
    }
    out_dir->cd();
    h_out->Write();
    fcomp->cd();
  };

  write_workflow_hist(signal_model.raw_hist_, "signal", signal_model.title_, signal_model.rate_, raw_dir, false);
  write_workflow_hist(signal_model.normalized_hist_, "signal", signal_model.title_, signal_model.rate_, normalized_dir);
  write_workflow_hist(signal_model.smoothed_hist_, "signal", signal_model.title_, signal_model.rate_, smoothed_dir);
  if(include_t0_) {
    write_workflow_hist(signal_model.t0_raw_hist_,
                        "signal_t0",
                        Form("%s t0 raw", signal_model.title_.Data()),
                        signal_model.rate_,
                        raw_dir,
                        false);
    write_workflow_hist(signal_model.t0_normalized_hist_,
                        "signal_t0",
                        Form("%s t0 normalized", signal_model.title_.Data()),
                        signal_model.rate_,
                        normalized_dir);
    write_workflow_hist(signal_model.t0_smoothed_hist_,
                        "signal_t0_exp",
                        Form("%s t0 exponential fit", signal_model.title_.Data()),
                        signal_model.rate_,
                        smoothed_dir);
  }
  for(auto& bkg : background_model) {
    write_workflow_hist(bkg.raw_hist_, bkg.name_, bkg.title_, bkg.rate_, raw_dir, false);
    write_workflow_hist(bkg.normalized_hist_, bkg.name_, bkg.title_, bkg.rate_, normalized_dir);
    write_workflow_hist(bkg.smoothed_hist_, bkg.name_, bkg.title_, bkg.rate_, smoothed_dir);
    if(include_t0_) {
      write_workflow_hist(bkg.t0_raw_hist_,
                          Form("%s_t0", bkg.name_.Data()),
                          Form("%s t0 raw", bkg.title_.Data()),
                          bkg.rate_,
                          raw_dir,
                          false);
      write_workflow_hist(bkg.t0_normalized_hist_,
                          Form("%s_t0", bkg.name_.Data()),
                          Form("%s t0 normalized", bkg.title_.Data()),
                          bkg.rate_,
                          normalized_dir);
      write_workflow_hist(bkg.t0_smoothed_hist_,
                          Form("%s_t0_exp", bkg.name_.Data()),
                          Form("%s t0 exponential fit", bkg.title_.Data()),
                          bkg.rate_,
                          smoothed_dir);
    }
  }

  fcomp->cd();
  TH1* hrmue_comp = new TH1F("rmue", "R_{#mue}", 1, 0., 1.);
  hrmue_comp->Fill(0.5, signal_br_);
  hrmue_comp->Write();
  TH1* hpot_comp = new TH1F("pot", "POT", 1, 0., 1.);
  hpot_comp->Fill(0.5, npot_);
  hpot_comp->Write();
  TH1* hlivetime_comp = new TH1F("livetime", "Livetime", 1, 0., 1.);
  hlivetime_comp->Fill(0.5, livetime_);
  hlivetime_comp->Write();
  TH1* hnmuons_comp = new TH1F("nmuons", "N(muon stops)", 1, 0., 1.);
  hnmuons_comp->Fill(0.5, nmuons_);
  hnmuons_comp->Write();
  fcomp->Close();

  // Construct the output workspace
  TFile* fout  = new TFile(ws_file  , "RECREATE");
  fout->cd();

  RooWorkspace ws("workspace", "workspace");
  ws.import(obs);
  ws.import(*data);
  if(do_2d_fit_) {
    ws.import(*signal_model.norm_);
    for(auto& bkg : background_model) ws.import(*bkg.norm_);
  } else {
    ws.import(*sig_pdf); ws.import(*signal_model.norm_);
    signal_model.hist_->Write();
    for(auto& bkg : background_model) {
      ws.import(*bkg.pdf_); ws.import(*bkg.norm_);
      bkg.hist_->Write();
    }
    // Combine profiles over the envelope function choice using this index; importing the
    // RooMultiPdf already pulls it in, so only add it if it somehow did not come along
    if(envelope.cat_ && !ws.cat(envelope.cat_->GetName())) ws.import(*envelope.cat_);
  }

  if(do_2d_fit_) {
    if(!include_t0_) {
      cout << __func__ << ": 2D fits require include_t0_ = true." << endl;
      return 1;
    }

    const TH1* signal_time_hist = preferred_time_hist(signal_model);
    if(!signal_time_hist) {
      cout << __func__ << ": Signal time histogram is missing, cannot build 2D inputs." << endl;
      return 1;
    }

    RooRealVar obs_t(Form("obs_t_%i", selection), "t0",
                     0.5 * (signal_time_hist->GetXaxis()->GetBinLowEdge(1) + signal_time_hist->GetXaxis()->GetBinUpEdge(signal_time_hist->GetNbinsX())),
                     signal_time_hist->GetXaxis()->GetBinLowEdge(1),
                     signal_time_hist->GetXaxis()->GetBinUpEdge(signal_time_hist->GetNbinsX()));
    obs_t.SetTitle("t0");
    ws.import(obs_t);

    if(!add_independent_2d_inputs(ws, signal_model, process, selection, signal_model.name_, obs, obs_t, figdir, fout)) {
      return 1;
    }
    for(auto& bkg : background_model) {
      if(!add_independent_2d_inputs(ws, bkg, process, selection, bkg.name_, obs, obs_t, figdir, fout)) {
        return 1;
      }
    }
  }

  // Add systematic uncertainties
  map<TString,map<TString, bool>> sys_map;
  if(include_sys_) {
    // Shape-based uncertainties
    for(int isys = 1; isys < mumep_ana::kMaxSystematics; ++isys) {
      TString sys_name = fSystematics.GetName(isys);
      if(sys_name == "") continue;
      const bool is_up = fSystematics.IsUp(isys);
      auto total_infos = background_model; total_infos.push_back(signal_model);
      for(auto& info : total_infos) {
        auto sys_pdf = read_model(info.name_, process, selection, tag, isys);
        if(!sys_pdf.pdf_ ) continue;
        if(!sys_pdf.norm_) continue;
        if(!sys_pdf.hist_) continue;
        fout->cd();
        sys_map[sys_name][info.name_] = true;
        sys_pdf.hist_->SetName(Form("%s_%s%s", sys_pdf.hist_->GetName(), sys_name.Data(), (is_up) ? "Up" : "Down"));
        sys_pdf.hist_->Write();
        sys_pdf.pdf_->SetName(Form("%s_%s%s", sys_pdf.pdf_->GetName(), sys_name.Data(), (is_up) ? "Up" : "Down"));
        ws.import(*sys_pdf.pdf_);
        sys_pdf.norm_->SetName(Form("%s_%s%s_norm", sys_pdf.pdf_->GetName(), sys_name.Data(), (is_up) ? "Up" : "Down"));
        ws.import(*sys_pdf.norm_);
      }
    }

    // Rate-based uncertainties
    auto rate_sys = rate_uncertainties(process);
    for(auto& sys : rate_sys) {
      // Add the signal model
      if(sys.isBeam || sys.process.Contains("signal")) {
        RooRealVar sig_impact(Form("%s_%i_signal_RateSys_%s", process.Data(), selection, sys.name.Data()),
                              Form("signal uncertainty from %s", sys.name.Data()), sys.value);
        ws.import(sig_impact);
      }
      // Add the background model
      for(auto& bkg : background_model) {
        bool include = false;
        include |= sys.isBeam && !bkg.name_.Contains("cosmic");
        include |= sys.process != "" && sys.process.Contains(bkg.name_);
        if(include) {
          cout << "Including rate uncertainty " << sys.name << " for process " << bkg.name_ << endl;
          RooRealVar bkg_impact(Form("%s_%i_%s_RateSys_%s", process.Data(), selection, bkg.name_.Data(), sys.name.Data()),
                                Form("%s uncertainty from %s", bkg.name_.Data(), sys.name.Data()), sys.value);
          ws.import(bkg_impact);
        }
      }
    }
  }

  // add a reference to the signal branching fraction and N(POT)/livetime used
  RooRealVar ref_signal_br("ref_signal_br", "BR(Signal) reference", signal_br_);
  ws.import(ref_signal_br);
  RooRealVar npot("npot", "N(POT)", npot_);
  ws.import(npot);
  RooRealVar livetime("livetime", "Livetime", livetime_);
  ws.import(livetime);
  RooRealVar nmuons("nmuons", "nmuons", nmuons_);
  ws.import(nmuons);
  RooRealVar sig_eff("signal_eff", "signal efficiency", signal_model.rate_ / n_signal_exp);
  ws.import(sig_eff);
  fout->cd();
  ws.Write();
  fout->Close();

  // Write the data card if requested
  if(write_card_) {
    // construct the card info list
    std::vector<card_info_t> card_info;
    card_info.push_back(card_info_t(signal_model.name_, signal_model.rate_, selection));
    for(auto& bkg : background_model) {
      const bool floating = envelope.pdf_ && bkg.pdf_ == envelope.pdf_;
      card_info.push_back(card_info_t(bkg.name_, bkg.rate_, selection, floating));
    }
    std::vector<TString> extra_lines;
    if(envelope.cat_) extra_lines.push_back(Form("%-10s discrete", envelope.cat_->GetName()));
    if(write_datacard(process, card_info, ws_file, sys_map, "", extra_lines)) {
      cout << __func__ << ": Data card writing failed!\n";
      return 1;
    }
  }


  // Create a cut-and-count selection to go along with this selection
  {
    TH1* sig_hist = sig_pdf->createHistogram("sig_hist", obs);
    sig_hist->Scale(signal_model.rate_ / sig_hist->Integral());
    TH1* bkg_hist = (TH1*) sig_hist->Clone("bkg_hist");
    bkg_hist->Reset();
    vector<TH1*> bkg_hists;
    for(auto& bkg : background_model) {
      TH1* h_tmp = bkg.pdf_->createHistogram(Form("bkg_hist_%s", bkg.name_.Data()), obs);
      h_tmp->Scale(bkg.rate_/h_tmp->Integral());
      bkg_hist->Add(h_tmp);
      bkg_hists.push_back(h_tmp);
    }
    // Determine the "best" observable region
    const int nbins = sig_hist->GetNbinsX();
    int best_low(1), best_high(1);
    double best_value = -1.;
    for(int bin_low = 1; bin_low <= nbins; ++bin_low) {
      for(int bin_high = bin_low; bin_high <= nbins; ++bin_high) {
        const double nsig = sig_hist->Integral(bin_low, bin_high);
        const double nbkg = bkg_hist->Integral(bin_low, bin_high);
        if(nbkg <= 0. || nsig <= 0.) continue;
        const double value = nsig / sqrt(nbkg + nsig); // S/sqrt(S+B)
        if(value > best_value) {
          best_value = value;
          best_low = bin_low;
          best_high = bin_high;
        }
      }
    }
    const double nsig  = sig_hist->Integral(best_low, best_high);
    const double nbkg  = bkg_hist->Integral(best_low, best_high);
    const double xlow  = sig_hist->GetXaxis()->GetBinLowEdge(best_low);
    const double xhigh = sig_hist->GetXaxis()->GetBinUpEdge(best_high);
    printf("----------------------------------------------\n");
    printf("Cut-and-count values:\n");
    printf("N(signal) = %.4f\nN(background) = %.4f\n", nsig, nbkg);
    printf("Region: %.2f - %.2f\n", xlow, xhigh);
    printf("----------------------------------------------\n");
    delete sig_hist;
    delete bkg_hist;

    // Make the card
    TString filler = std::string((background_model.size()+2)*10 + 15, '-');
    TString outname = ws_file;
    outname.ReplaceAll(".root", "_cc.txt");
    if(outname.Contains("/")) outname = outname(outname.Last('/')+1, outname.Sizeof());
    outname.ReplaceAll("workspace", "combine");
    outname = "datacards/" + outname;

    // construct the card info list
    std::vector<card_info_t> card_info;
    card_info.push_back(card_info_t(signal_model.name_, nsig, selection));
    for(size_t index = 0; index < background_model.size(); ++index) {
      const auto& bkg = background_model[index];
      const double nbkg_i  = bkg_hists[index]->Integral(best_low, best_high);
      card_info.push_back(card_info_t(bkg.name_, nbkg_i, selection));
    }
    const double nexp = nsig + nbkg;
    const double sig_eff_cc = nsig / n_signal_exp;
    if(write_counting_datacard(process, card_info, outname, (int) nexp, npot_, livetime_, nmuons_, signal_br_, sig_eff_cc, xlow, xhigh)) {
      cout << __func__ << ": Counting data card writing failed!\n";
      return 1;
    }
  }

  // Print out summary info:
  printf("----------------------------------------------\n");
  printf("N(POT)   = %.1e\n", npot_);
  printf("Livetime = %.1e s\n", livetime_);
  printf("N(muons) = %.1e\n", nmuons_);
  printf("%-10s (%-20s): Rate = %.4f\n", signal_model.name_.Data(), signal_model.title_.Data(), signal_model.rate_);
  for(auto& bkg : background_model)
    printf("%-10s (%-20s): Rate = %.4f\n", bkg.name_.Data(), bkg.title_.Data(), bkg.rate_);
  printf("----------------------------------------------\n");
  return 0;
}
