#ifndef __MUMEP_ANA_ANALYSIS_MAKE_PAPER_PLOTS__
#define __MUMEP_ANA_ANALYSIS_MAKE_PAPER_PLOTS__

// Make plots targeted for the paper
#include "make_plots.C"

//---------------------------------------------------------------------------------------------------
int make_paper_plots(const bool mumem = false, TString tag = "paper_r0101") {
  combine_rmc_ = mumem;
  if(plotter_) {
    delete plotter_;
    plotter_ = nullptr;
  }
  plotter_ = new Plotter();
  TString figdir = "figures/plots/";
  if(use_evtana_) figdir += "evtana_";
  figdir += (mumem) ? "mumem" : "mumep";
  if(tag != "") figdir += "_" + tag;
  plotter_->figdir_ = figdir;
  plotter_->signal_ = (mumem) ? "mumem" : "mumep";
  plotter_->configure_style(true, 3, true, 2, {"cosmic"});
  if(mumem) plotter_->bkgs_ = {"rpc_ext", "rpc_int", "pbar", "rmc_ext", "cosmic", "dio"};
  else      plotter_->bkgs_ = {"rpc_ext", "rpc_int", "pbar", "cosmic", "rmc_ext", "rmc_int"};
  if(use_evtana_) {
    plotter_->configure_for_evtana();
    if(mumem)
      plotter_->bkgs_ = {"pbar", "rpc_int", "rpc_ext", "rmc_ext_0n", "rmc_ext_1n", "rmc_int_0n", "rmc_int_1n", "cosmic", "dio"}; // only some are available
    else
      plotter_->bkgs_ = {"pbar", "rpc_int", "rpc_ext", "cosmic", "rmc_ext_0n", "rmc_ext_1n", "rmc_int_0n", "rmc_int_1n"}; // only some are available
    hist_mode_ = 1;
  }
  if(plotter_->init("", tag)) {
    delete plotter_;
    plotter_ = nullptr;
    return 1;
  }

  plotter_->stack_signal_ = 0;
  printf("------------------------------------------------------\n");
  printf("Normalization: N(POT) = %.2e, livetime = %.2e\n", npot_, livetime_);
  printf("------------------------------------------------------\n");

  int status(0);
  TCanvas* c;
  if(!mumem) signal_br_ = 1.7e-12;
  const double base_br(signal_br_);
  plotter_->update_signal_br(signal_br_);
  plotter_->use_offsets_ = false; //don't use control regions for initial counts
  plotter_->ratio_plot_ = false;

  for(int logy = 0; logy < 2; ++logy) {
    if(mumem) {
      c = plotter_->print_stack(plot_t("p_2", "trk", 75, 5, 100., 110., 1., -1., logy, false, "p", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
    } else {
      c = plotter_->print_stack(plot_t("p_2", "trk",  5, 5, 87., 97., 1., -1., logy, false, "p", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
      c = plotter_->print_stack(plot_t("p_2", "trk", 40, 5, 87., 97., 1., -1., logy, false, "p", "MeV/c")); if(!c) ++status; else Empty_Canvas(c);
    }
  }
  return status;
}

#endif
