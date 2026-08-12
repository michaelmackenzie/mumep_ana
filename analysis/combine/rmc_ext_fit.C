// Fit the RMC external spectrum shape
#ifndef __CONVANA_ANALYSIS_RMCEXTFIT__
#define __CONVANA_ANALYSIS_RMCEXTFIT__

#include "component_fit.C"

int rmc_ext_fit(TString process = "mumem", int selection = 20, TString tag = "", TString pdf_type = "default",
                TString tail_model = "default", std::vector<int> shape_sets = {},
                std::vector<int> control_region_sets = {}) {
  return fit_component_model(process, selection, tag,
                             "rmc_ext", "RMC (external)",
                             pdf_type, -1, tail_model,
                             shape_sets, control_region_sets, "primary");
}

#endif
