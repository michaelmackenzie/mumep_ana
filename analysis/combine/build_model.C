// Build the statistical model

#include "../tools/types.C"
#include "signal_model.C"
#include "background_model.C"
#include "../physics.C"
#include "systematics.C"
#include "../tools/write_datacard.C"

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
  sys.push_back(RateUnc_t("lumi"  , 0.1  , true           ));
  sys.push_back(RateUnc_t("cosmic", 0.2  , false, "cosmic"));
  sys.push_back(RateUnc_t("dio"   , 0.025, false, "dio"   ));
  sys.push_back(RateUnc_t("rpc"   , 0.27 , false, "rpc"   ));
  sys.push_back(RateUnc_t("pbar"  , 1.   , false, "pbar"  ));
  return sys;
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

  // Retrieve the signal data
  auto signal_model     = read_model          ("signal", process, selection, tag);
  auto background_model = get_background_model(obs     , process, selection, tag);
  auto data             = get_data            (obs     , process, selection, tag);

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

  RooArgList bkg_pdfs;
  RooArgList bkg_rates;
  for(auto& bkg : background_model) {
    if(!bkg.pdf_) continue;
    bkg_pdfs.add(*bkg.pdf_);
    bkg_rates.add(*(new RooRealVar(Form("%s_rate", bkg.name_.Data()), "", bkg.rate_)));
  }
  RooAddPdf tot_bkg("tot_bkg", "Total background", bkg_pdfs, bkg_rates);


  // Generate toy data
  if(!data) {
    for(auto& bkg : background_model) {
      auto gen_data = bkg.pdf_->generateBinned(obs, bkg.rate_);
      if(!data) {
        data = gen_data;
        data->SetName("data_obs");
      } else if(gen_data) {
        data->add(*gen_data);
      } else {
        cout << __func__ << ": Gen data for process " << bkg.name_.Data() << " is null!\n";
      }
    }
  }

  // Draw the inputs
  if(print_) {
    const double signal_scale = (process == "mumem") ? 20. : 1.7e3;
    TString figdir = Form("figures/%s%s", process.Data(), (tag != "") ? ("_"+tag).Data() : "");
    gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir.Data(), figdir.Data()));

    TCanvas* c = new TCanvas("c_model", "c_model", 1200, 1000);
    auto frame = obs.frame();
    frame->SetTitle(Form("%s model", signal_model.title_.Data()));
    frame->SetXTitle("Momentum (MeV/c)");

    // data the data
    data->plotOn(frame, RooFit::Name("data"));

    // draw the signal
    sig_pdf->plotOn(frame, RooFit::Name(signal_model.name_),
                    RooFit::LineColor(signal_model.color_), RooFit::FillColor(signal_model.color_),
                    RooFit::FillStyle(3005),
                    RooFit::Normalization(signal_scale*signal_model.rate_, RooAbsReal::NumEvent));

    // draw the backgrounds
    tot_bkg.plotOn(frame, RooFit::Invisible(), RooFit::Name("bkg"));
    for(auto& bkg : background_model) {
      if(!bkg.pdf_) {
        cout << __func__ << ": Background " << bkg.name_.Data() << " has an undefined PDF!\n";
        return 1;
      }
      bkg.pdf_->plotOn(frame, RooFit::Name(bkg.name_),
                       RooFit::LineColor(bkg.color_), RooFit::Normalization(bkg.rate_, RooAbsReal::NumEvent));
      cout << "bkg " << bkg.name_.Data() << " norm " << bkg.rate_ << endl;
    }
    frame->SetYTitle("");
    frame->Draw();

    c = plot_fit_frame(frame, obs, "Momentum (MeV/c)", Form("Events / %.1f MeV/c", bin_width_), "data", "bkg", npot_, livetime_, nmuons_);
    auto pad1 = (TPad*) c->GetPrimitive("pad1");

    // add a legend
    TLegend* leg = new TLegend((pad1) ? pad1->GetLeftMargin() : 0.1, 0.75, (pad1) ? 1. - pad1->GetRightMargin() : 0.9, (pad1) ? 1. - pad1->GetTopMargin() : 0.9);
    leg->SetNColumns(3); leg->SetLineWidth(0); leg->SetFillColor(0); leg->SetLineColor(0); leg->SetFillStyle(0);
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
  ws.import(*sig_pdf); ws.import(*signal_model.norm_);
  signal_model.hist_->Write();
  for(auto& bkg : background_model) {
    ws.import(*bkg.pdf_); ws.import(*bkg.norm_);
    bkg.hist_->Write();
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
    for(auto& bkg : background_model) card_info.push_back(card_info_t(bkg.name_, bkg.rate_, selection));
    if(write_datacard(process, card_info, ws_file, sys_map)) {
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
  printf("%-10s (%-15s): Rate = %.4f\n", signal_model.name_.Data(), signal_model.title_.Data(), signal_model.rate_);
  for(auto& bkg : background_model)
    printf("%-10s (%-15s): Rate = %.4f\n", bkg.name_.Data(), bkg.title_.Data(), bkg.rate_);
  printf("----------------------------------------------\n");
  return 0;
}
