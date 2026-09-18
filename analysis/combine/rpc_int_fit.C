// Fit the RPC internal spectrum shape
#ifndef __CONVANA_ANALYSIS_RPCINTFIT__
#define __CONVANA_ANALYSIS_RPCINTFIT__

#include "component_fit.C"

// hist_file: optional alternate histogram file to fit, instead of the default one in hist_path_
int rpc_int_fit(TString process = "mumem", int selection = 20, TString tag = "", TString pdf_type = "default",
                std::vector<int> shape_sets = {}, std::vector<int> control_region_sets = {},
                TString hist_file = "") {
  return fit_component_model(process, selection, tag,
                             "rpc_int", "RPC (internal)",
                             pdf_type, -1, "none",
                             shape_sets, control_region_sets, "primary", hist_file);
}

#endif
