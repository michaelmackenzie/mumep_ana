// Exposure of the analysis macros. The values are set by combine/samples.py (export_to_root, from samples.Exposure and
// Mu2eEvtAna's plotter/physics/Mu2ePhysics.C) before a macro runs: run the macros with run_macro.py
#ifndef __CONVANA_ANALYSIS_PHYSICS__
#define __CONVANA_ANALYSIS_PHYSICS__
#include "defaults.C"

double  livetime_           = 0.    ; // livetime (s)
double  npot_               = 0.    ; // N(POT)
double  nmuons_             = 0.    ; // N(stopped muons)
double  nmuons_per_pot_     = 0.    ; // stopped muons per POT (the stopping target's with exposure_target_)
double  exposure_reference_ = 1.e-15; // signal branching fraction the exposure is normalized to
TString exposure_name_      = ""    ; // samples.Exposure name: run1a, run2, mds3c, or mds3d
bool    exposure_two_bb_    = false ; // the Run 1A/2 exposure at the 2BB duty cycle and intensity
TString exposure_target_    = ""    ; // stopping-target design of samples.TARGETS ("" for the nominal target)
double  signal_br_          = 1.e-15; // example signal branching fraction for normalization
double  muon_capture_fraction_ = 0. ; // on aluminum (used by legacy/build_model.C)

// Exposure named by a dataset/tag string ("" if it names an exposure that isn't defined)
TString exposure_from_tag(TString tag) {
  tag.ToLower();
  TString name = "run1a";
  if(tag.Contains("run2") || tag.Contains("run_2")) name = "run2";
  if(tag.Contains("mds")) {
    if     (tag.Contains("mds3c")) name = "mds3c";
    else if(tag.Contains("mds3d")) name = "mds3d";
    else                           name = ""; // MDS1/2 ensembles: Stntuple inputs, no longer defined
  }
  if(tag.Contains("mdc2025") || tag.Contains("run1a")) name = "run1a";
  return name;
}

// Check that samples.py set the exposure the tag names (Run 1A unless it names run2, mds3c, or mds3d) with the 2BB
// setting of use_2bb_, and set the signal branching fraction to the exposure's reference
bool init_physics(TString tag) {
  if(exposure_name_ == "") {
    printf("init_physics: The exposure is not set. It is set by combine/samples.py (export_to_root):\n"
           "  run the macro with run_macro.py, e.g. ./run_macro.py 'make_plots.C(true, {75})'\n");
    return false;
  }
  const TString name = exposure_from_tag(tag);
  if(name == "") {
    printf("init_physics: Tag \"%s\" names an exposure that isn't defined (see samples.Exposure)\n", tag.Data());
    return false;
  }
  if(name != exposure_name_) {
    printf("init_physics: Tag \"%s\" uses the %s exposure, but the %s exposure is set (run_macro.py --exposure %s)\n",
           tag.Data(), name.Data(), exposure_name_.Data(), name.Data());
    return false;
  }
  if(use_2bb_ != exposure_two_bb_) {
    printf("init_physics: use_2bb_ = %i, but the exposure was set with 2BB = %i (run_macro.py --2bb)\n",
           use_2bb_, exposure_two_bb_);
    return false;
  }
  signal_br_ = exposure_reference_;
  tag.ToLower();
  if(tag.Contains("mdc2025") || tag.Contains("run1a")) signal_br_ = 1.e-13; // plotting convention of these tags
  return true;
}

#endif
