// Fit the RMC internal spectrum shape
#ifndef __CONVANA_ANALYSIS_RMCINTFIT__
#define __CONVANA_ANALYSIS_RMCINTFIT__

#include "component_fit.C"

int rmc_int_fit(TString process = "mumem", int selection = 20, TString tag = "", TString pdf_type = "default",
                TString tail_model = "default", std::vector<int> shape_sets = {},
                std::vector<int> control_region_sets = {}) {
  return fit_component_model(process, selection, tag,
                             "rmc_int", "RMC (internal)",
                             pdf_type, -1, tail_model,
                             shape_sets, control_region_sets, "primary");
}

#endif
