#ifndef __MUMEP_ANA_ANALYSIS_MAKE_PLOTS__
#define __MUMEP_ANA_ANALYSIS_MAKE_PLOTS__

// Make standard plots
#include "plotter/Plotter.C"
static Plotter* plotter_ = nullptr;

//---------------------------------------------------------------------------------------------------
// Print the relevant process codes in the model
int print_proc_info(const int selection) {
  if(!plotter_) return 1;
  TString hist = "primary_code";
  TString type = "evt";
  auto backgrounds = plotter_->get_histograms(hist, type, selection,  1);
  auto signals     = plotter_->get_histograms(hist, type, selection, -1);
  auto datas       = plotter_->get_histograms(hist, type, selection,  0);

  printf("Process information for set %4i:\n", selection);
  for(auto h : backgrounds) {
    cout << "Background " << h->GetTitle() << ":\n";
    for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
      if(h->GetBinContent(ibin) > 0.) printf("  Process %3i: rate = %10g\n", (int) (h->GetBinCenter(ibin) + 0.5), h->GetBinContent(ibin));
    }
  }

  for(auto h : signals) {
    cout << "Signal " << h->GetTitle() << ":\n";
    for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
      if(h->GetBinContent(ibin) > 0.) printf("  Process %3i: rate = %10g\n", (int) (h->GetBinCenter(ibin) + 0.5), h->GetBinContent(ibin));
    }
  }

  for(auto h : datas) {
    cout << "Data " << h->GetTitle() << ":\n";
    for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
      if(h->GetBinContent(ibin) > 0.) printf("  Process %3i: rate = %10g\n", (int) (h->GetBinCenter(ibin) + 0.5), h->GetBinContent(ibin));
    }
  }

  return 0;
}

//---------------------------------------------------------------------------------------------------
// Print the dataset counts for a set
int print_dataset_info(const int selection) {
  if(!plotter_) return 1;
  TString hist = "p";
  TString type = "trk";

  printf("Dataset information for set %4i:\n", selection);
  for(auto& input : plotter_->data_) {
    int set_offset = input.set_offset_;
    if(set_offset > 0 && selection > set_offset && (selection / set_offset) % 2 == 1) set_offset = 0; // in the control region
    const int hist_set = selection; //(plotter_->use_offsets_) ? selection + set_offset : selection;
    TString hist_path = Form("%sHist/%s_%i/%s", plotter_->dir_path_.Data(), type.Data(), hist_set, hist.Data());
    TH1* h = (TH1*) input.f_->Get(hist_path);
    if(!h) continue;
    double counts = h->Integral(0, h->GetNbinsX()+1);
    double rate = counts*input.norm_*input.scale_;
    if(plotter_->uses_livetime_scale(input.name_)) rate *= livetime_;
    else                                           rate *= npot_    ;
    printf("%-25s (%10s): %10.4g (%10.4g)\n", input.label_.Data(), input.name_.Data(), rate, counts);
  }

  return 0;
}

//---------------------------------------------------------------------------------------------------
int make_plots(const bool mumem = true, vector<int> sets = {7, 10, 20}, TString dataset = "mds3c", TString tag = "r0102") {
  if(plotter_) {
    delete plotter_;
    plotter_ = nullptr;
  }
  plotter_ = new Plotter();
  TString figdir = "figures/plots/";
  if(use_evtana_) figdir += "evtana_";
  figdir += (mumem) ? "mumem" : "mumep";
  if(dataset != "") figdir += "_" + dataset;
  if(tag != "") figdir += "_" + tag;
  plotter_->figdir_ = figdir;
  plotter_->signal_ = (mumem) ? "mumem" : "mumep";
  plotter_->configure_style(true, 3, true, 2, {"cosmic"});
  if(mumem) plotter_->bkgs_ = {"rpc_ext", "rpc_int", "pbar", "rmc_ext", "cosmic", "dio"};
  else      plotter_->bkgs_ = {"rpc_ext", "rpc_int", "pbar", "cosmic", "rmc_ext", "rmc_int"};
  if(use_evtana_) {
    plotter_->configure_for_evtana();
    plotter_->bkgs_ = {"rpc_int", "rpc_ext", "rmc_ext_0n", "rmc_ext_1n", "rmc_int_0n", "rmc_int_1n", "cosmic", "dio"}; // only some are available
    hist_mode_ = 1;
  }
  if(plotter_->init(dataset, tag)) {
    delete plotter_;
    plotter_ = nullptr;
    return 1;
  }
  if(dataset.BeginsWith("mds1")) plotter_->stack_signal_ = 1; // include the signal in the stacks
  printf("------------------------------------------------------\n");
  printf("Normalization: N(POT) = %.2e, livetime = %.2e\n", npot_, livetime_);
  printf("------------------------------------------------------\n");

  int status(0);
  const bool mds = dataset.Contains("mds");
  TCanvas* c;
  if(!mumem) signal_br_ = 1.7e-13;
  const double base_br(signal_br_);
  plotter_->update_signal_br(signal_br_);
  plotter_->use_offsets_ = false; //don't use control regions for initial counts
  print_proc_info(0);
  print_proc_info(6);
  if(use_evtana_) print_dataset_info(60);
  for(int set : sets) {
    if(set < 0) continue;
    signal_br_ = (set%20 == 1) ? base_br/10. : base_br*100.;
    if(mds)  signal_br_ = base_br;
    plotter_->update_signal_br(signal_br_);
    double p_min((mumem) ? (mds && set < 20) ? 95. : 100. : 87.), p_max((mumem) ? 110. : 97.);
    if(use_evtana_ && set == 60) {p_min = 97.; p_max = 110.;}
    if(set < 10)       plotter_->use_offsets_ = false; //if using control region offsets
    else if(set == 30) plotter_->use_offsets_ = false;
    else if(set == 60) plotter_->use_offsets_ = false;
    else if(set >= 150 && set <= 180) plotter_->use_offsets_ = false;// PID sets
    else if(tag.Contains("mdc2025")) plotter_->use_offsets_ = false; // Don't use in Run 1a config
    // else               plotter_->use_offsets_ = true ; // default to using them
    for(int logy = 0; logy < 2; ++logy) {
      if(mumem) {c = plotter_->print_stack(plot_t("p"   , "trk", set, 1,-120., -95.,  1., -1., logy, false, "q*p", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);}
      else      {c = plotter_->print_stack(plot_t("p"   , "trk", set, 1,  75., 120.,  1., -1., logy, false, "p", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);}
      c = plotter_->print_stack(plot_t("p_2"            , "trk", set, 5,p_min, p_max, 1., -1., logy, false, "p", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
      // c = plotter_->print_stack(plot_t("p_corr"         , "trk", set, 5,p_min, p_max, 1., -1., logy, false, "p", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("obs"            , "trk", set, 1,p_min, p_max, 1., -1., logy, false, "p", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("pt"             , "trk", set, 1,  50.,  110., 1., -1., logy, false, "p_{T}", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("t0"             , "trk", set, 5, 500., 1700., 1., -1., logy, false, "t_{0}", "ns")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("d0"             , "trk", set, 5, -200., 200., 1., -1., logy, false, "D_{0}", "mm")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("rMax"           , "trk", set, 5, 350., 800. , 1., -1., logy, false, "R_{max}", "mm")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("tanDip"         , "trk", set, 5, 1.  , -1.  , 1., -1., logy, false, "tan(dip)", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("cosTheta"       , "trk", set, 2, 0.  ,  1.  , 1., -1., logy, false, "cos(#theta)", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("trkQual"        , "trk", set, 2, 0.  ,  1.  , 1., -1., logy, false, "track quality", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("trkQual_1"      , "trk", set, 2, 0.  ,  1.  , 1., -1., logy, false, "track quality", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("pid"            , "trk", set, 2, 0.  ,  1.  , 1., -1., logy, false, "PID", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("trkpid"         , "trk", set, 2, 0.  ,  1.  , 1., -1., logy, false, "tracker PID", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("cosmic_id"      , "trk", set, 2, 0.  ,  1.  , 1., -1., logy, false, "Cosmic ID", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("fitCons"        , "trk", set, 5, 0.  ,  1.  , 1., -1., logy, false, "p(#chi^2)", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("fitCons_log"    , "trk", set, 5, -5. ,  1.  , 1., -1., logy, false, "log(p(#chi^2))", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("clusterE"       , "trk", set, 4, 0.  , 110. , 1., -1., logy, false, "Cluster energy", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
      // c = plotter_->print_stack(plot_t("clusterDisk"    , "trk", set, 0, 1.  ,  -1. , 1., -1., logy, false, "Cluster disk", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("dt"             , "trk", set, 2, -10.,  5.  , 1., -1., logy, false, "#Deltat", "ns")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("nActive"        , "trk", set, 1,   1.,  -1. , 1., -1., logy, false, "N(active hits)", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("tzslopesig"     , "trk", set, 2, -10.,  10. , 1., -1., logy, false, "TZ slope / #sigma", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("tzsloperatio"   , "trk", set, 2,  -5.,   5. , 1., -1., logy, false, "TZ slope / expected", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("pST_diff"       , "trk", set, 2, 0., 5.     , 1., -1., logy, false, "#Deltap(ST)", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("us_dt"          , "trk", set, 2,-50.,250.   , 1., -1., logy, false, "#Deltat(upstream track)", "ns")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("us_dp"          , "trk", set, 2,-10., 5.    , 1., -1., logy, false, "#Deltap(upstream track)", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("crv_min_deltat" , "trk", set, 2,-100., 100. , 1., -1., logy, false, "#Deltat(CRV)", "ns")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("st_inters"      , "trk", set, 0, 0., 15.    , 1., -1., logy, false, "N(ST intersections)", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("track_id"       , "trk", set, 1, 0., 14.    , 1., -1., logy, false, "", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("track_exl_id"   , "trk", set, 1, 0., 14.    , 1., -1., logy, false, "", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("inst_lumi"      , "evt", set, 5, 0., 1.e8   , 1., -1., logy, false, "N(POT)", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("trig_bits"      , "evt", set, 1, 1, -1.     , 1., -1., logy, false, "Trigger bits", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("event_weight_log", "evt", set, 5, 1, -1.    , 1., -1., logy, false, "log10(event weight)", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("MC_PDG_0"       , "trk", set, 1, 1., -1.    , 1., -1., logy, false, "PDG ID", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("MC_trajectory"  , "trk", set, 1, 1., -1.    , 1., -1., logy, false, "MC Trajectory", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("MC_PFront"      , "trk", set, 2, 80., 120.  , 1., -1., logy, false, "MC p(Front)", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("MC_PSTOut"      , "trk", set, 2, 1., -1.    , 1., -1., logy, false, "MC p(ST)", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("MC_GenE"        , "trk", set, 2, 80., 120.  , 1., -1., logy, false, "MC E(gen)", "MeV")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("primary_code"   , "evt", set, 0, 1., -1.    , 1., -1., logy, false, "Process code", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("primary_type"   , "evt", set, 0, 1., -1.    , 1., -1., logy, false, "", "")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("dP"             , "trk", set, 2,-5.,  5.    , 1., -1., logy, false, "p - p(MC front)", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);

      if(set == 67 && !logy) {
        for(int eff = 0; eff < 2; ++eff) {
          c = plotter_->print_roc(plot_t("tanDip"   , "trk", set), false, eff == 1); if(!c) ++status; //else Empty_Canvas(c);
          c = plotter_->print_roc(plot_t("cosmic_id", "trk", set),  true, eff == 1); if(!c) ++status; //else Empty_Canvas(c);
        }
      }
      if((set >= 150 && set <= 180) || set == 30) { // PID/CRV sets
        c = plotter_->print_stack(plot_t("crv_deltat_crv"    , "trk", set, 2,-100., 100. , 1., -1., logy, false, "#Deltat(CRV)", "ns")); if(!c) ++status; else Empty_Canvas(c);
        c = plotter_->print_stack(plot_t("crv_deltat_st"     , "trk", set, 2,-100., 100. , 1., -1., logy, false, "#Deltat(CRV through ST)", "ns")); if(!c) ++status; else Empty_Canvas(c);
        c = plotter_->print_stack(plot_t("crv_deltat_calo"   , "trk", set, 2,-100., 100. , 1., -1., logy, false, "#Deltat(CRV through Calo)", "ns")); if(!c) ++status; else Empty_Canvas(c);
        c = plotter_->print_stack(plot_t("crv_deltat_calo_mu", "trk", set, 2,-100., 100. , 1., -1., logy, false, "#Deltat(CRV through Calo (muon))", "ns")); if(!c) ++status; else Empty_Canvas(c);
        c = plotter_->print_stack(plot_t("crv_deltat_extrap" , "trk", set, 2,-100., 100. , 1., -1., logy, false, "#Deltat(CRV through extrapolation)", "ns")); if(!c) ++status; else Empty_Canvas(c);
        c = plotter_->print_component(plot_t("crv_deltat_crv"     , "trk", set, 2,-100., 100. , 1., -1., logy, false, "#Deltat(CRV)", "ns"), "cosmic"); handle_canvas(c, status);
        c = plotter_->print_component(plot_t("crv_deltat_st"      , "trk", set, 2,-100., 100. , 1., -1., logy, false, "#Deltat(CRV through ST)", "ns"), "cosmic"); handle_canvas(c, status);
        c = plotter_->print_component(plot_t("crv_deltat_calo"    , "trk", set, 2,-100., 100. , 1., -1., logy, false, "#Deltat(CRV through Calo)", "ns"), "cosmic"); handle_canvas(c, status);
        c = plotter_->print_component(plot_t("crv_deltat_calo_mu" , "trk", set, 2,-100., 100. , 1., -1., logy, false, "#Deltat(CRV through Calo (muon))", "ns"), "cosmic"); handle_canvas(c, status);
      }

      // // Individual components
      // c = plotter_->print_component(plot_t("t0", "trk", set, 5, 500., 1200., 1., -1., logy, false, "t_{0}", "ns"), "rpc_ext"); if(!c) ++status; else Empty_Canvas(c);
      // c = plotter_->print_component(plot_t("t0", "trk", set, 5, 500., 1200., 1., -1., logy, false, "t_{0}", "ns"), "rpc_int"); if(!c) ++status; else Empty_Canvas(c);

      // Systematics
      if(set == 20 || set == 34 || set == 35) {
        // plotter_->stack_signal_ = false;
        c = plotter_->print_systematic(plot_t("obs", "trk", set, 5, 100., 110., 1., -1., logy, false, "p", "MeV/c"), 1, 2 ); if(!c) ++status; else Empty_Canvas(c);
      }

    }

    // print summary info about the process codes:
    if(use_evtana_) print_dataset_info(set);
    else            print_proc_info   (set);
  }
  return status;
}

#endif
