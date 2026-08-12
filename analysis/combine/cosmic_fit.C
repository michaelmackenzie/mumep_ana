// Fit the DIO spectrum shape
#ifndef __CONVANA_ANALYSIS_COSMICFIT__
#define __CONVANA_ANALYSIS_COSMICFIT__

#include "component_fit.C"

bool use_control_region_ = false; // whether or not to take the distribution from the CRV tagged region
bool fit_total_          = true ; // fit the combined veto + signal regions to get more stable results

int cosmic_fit(TString process = "mumem", int selection = 20, TString tag = "", const int isys = -1,
               TString pdf_type = "default", std::vector<int> control_region_sets = {}) {
  TString fit_input = "primary";
  if(fit_total_) fit_input = "total";
  else if(use_control_region_) fit_input = "control";

  return fit_component_model(process, selection, tag,
                             "cosmic", "Cosmic ray",
                             pdf_type, isys, "none",
                             {}, control_region_sets, fit_input);
}

#endif
