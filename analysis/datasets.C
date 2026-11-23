#ifndef __CONVANA_ANALYSIS_DATASETINFO__
#define __CONVANA_ANALYSIS_DATASETINFO__
#include "defaults.C"
#include "tools/types.C"
#include "physics.C"
#include <map>
#include <vector>

// The dataset list lives in combine/samples.py (DATASETS, TARGETS, DATA). Its export_to_root fills datasets_, and the
// exposure (physics.C) and histogram directories (defaults.C), before a macro runs; run the macros with run_macro.py:
//   ./run_macro.py 'make_plots.C(true, {75})'
//   ./run_macro.py --target st49hole1to5 'optimize_selection_fast.C(true, 81, false, "st49hole1to5", 1e3)'
// The rates (DatasetInfo_t::theory_) are per POT (per unit branching fraction for the signal, 1 for the cosmic rays,
// whose N(gen) is the generated livetime), at the exposure's stops per POT (nmuons_per_pot_).
std::map<TString, DatasetInfo_t> datasets_;
const char* hist_func_    = "cnv_ana"; // for histogram file naming
bool        combine_rpc_  = true  ; // Combine RPC samples into one label
bool        combine_rmc_  = false ; // Combine RMC samples into one label
bool        do_knockouts_ = false ; // Separate distict RMC knockout channels

// Check that samples.py filled the dataset list
bool require_dataset_info() {
  if(!datasets_.empty()) return true;
  printf("require_dataset_info: The dataset list is empty. It is filled by combine/samples.py (export_to_root):\n"
         "  run the macro with run_macro.py, e.g. ./run_macro.py 'make_plots.C(true, {75})'\n");
  return false;
}

// Check that the dataset list was made for the stopping-target design ("" for the nominal samples). init_physics checks
// the rest of the exposure.
bool check_dataset_exposure(const TString target = "") {
  if(!require_dataset_info()) return false;
  if(target != exposure_target_) {
    printf("check_dataset_exposure: The dataset list is for target \"%s\" but target \"%s\" is requested"
           " (run_macro.py --target)\n", exposure_target_.Data(), target.Data());
    return false;
  }
  return true;
}

DatasetInfo_t get_dataset_info(TString name) {
  if(!require_dataset_info()) return DatasetInfo_t(1., 0, 0., 1., 1., "");
  if(datasets_.count(name) != 0) return datasets_[name];
  cout << __func__ << ": No dataset with name " << name.Data() << " found!\n";
  return DatasetInfo_t(1., 0, 0., 1., 1., "");
}

//---------------------------------------------------------------------------------------------------------------------------
// RMC neutron-knockout handling
// The mu- --> e+ model splits each RMC source into its 0 and 1 neutron knockout contributions,
// while mu- --> e- keeps a single component per source (using the 0n sample).
bool split_rmc_knockout(const TString process) {
  return process == "mumep";
}

// Component names to include for an RMC source ("rmc_ext" or "rmc_int")
std::vector<TString> rmc_components(const TString process, const TString base) {
  if(split_rmc_knockout(process)) return {base + "_0n", base + "_1n"};
  return {base};
}

// Knockout label ("0n"/"1n") of a component name, "" if it isn't knockout-split
TString rmc_knockout(const TString component) {
  if(component.EndsWith("_0n")) return "0n";
  if(component.EndsWith("_1n")) return "1n";
  return "";
}

// Component name for an RMC source and knockout label
TString rmc_component_name(const TString base, const TString knockout) {
  return (knockout == "") ? base : base + "_" + knockout;
}

//---------------------------------------------------------------------------------------------------------------------------
void set_style(const TString name, TString& title, int& color) {
  title = name;
  color = kRed;
  if(name == "signal") {
    color = kBlue;
    title = "Signal";
  } else if(name == "mumem") {
    color = kBlue;
    title = "#mu^{-}#rightarrowe^{-}";
  } else if(name == "mumep") {
    color = kBlue;
    title = "#mu^{-}#rightarrowe^{+}";
  } else if(name == "rmc") { // signal version
    title = "RMC";
    color = kRed-7;
  } else if(name == "dio") {
    title = "DIO";
    color = kMagenta-10;
  } else if(name == "ipa_dio") {
    title = "IPA DIO";
    color = kMagenta-8;
  } else if(name.BeginsWith("cosmic")) {
    title = "Cosmic ray";
    color = kAzure-4;
  } else if(name == "pbar") {
    title = "Antiproton";
    color = kMagenta-3;
  } else if(name == "rpc_ext") {
    title = (combine_rpc_) ? "RPC" : "RPC (external)";
    color = (combine_rpc_) ? kGreen-6 : kGreen-6;
  } else if(name == "rpc_int") {
    title = (combine_rpc_) ? "RPC" : "RPC (internal)";
    color = (combine_rpc_) ? kGreen-6 : kGreen-9;
  } else if(name.BeginsWith("rmc_ext")) {
    const TString knockout = rmc_knockout(name);
    title = (combine_rmc_) ? "RMC" : ((!do_knockouts_ || knockout == "") ? "RMC (external)" : Form("RMC (%s-ext)", knockout.Data()));
    color = (do_knockouts_ && knockout == "1n") ? kRed-5 : kRed-7;
  } else if(name.BeginsWith("rmc_int")) {
    const TString knockout = rmc_knockout(name);
    title = (combine_rmc_) ? "RMC" : ((!do_knockouts_ || knockout == "") ? "RMC (internal)" : Form("RMC (%s-int)", knockout.Data()));
    color = (do_knockouts_ && knockout == "1n") ? kRed-8 : kRed-9;
  } else if(name.BeginsWith("mnbs")) {
    title = "Pileup";
    color = kRed;
  }
}

#endif
