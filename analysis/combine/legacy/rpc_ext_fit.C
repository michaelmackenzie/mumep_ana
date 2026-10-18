// Fit the RPC external spectrum shape
#ifndef __CONVANA_ANALYSIS_RPCEXTFIT__
#define __CONVANA_ANALYSIS_RPCEXTFIT__

#include "component_fit.C"

// hist_file: optional alternate histogram file to fit, instead of the default one in hist_path_
int rpc_ext_fit(TString process = "mumem", int selection = 20, TString tag = "", TString pdf_type = "default",
                std::vector<int> shape_sets = {}, std::vector<int> control_region_sets = {},
                TString hist_file = "") {
  return fit_component_model(process, selection, tag,
                             "rpc_ext", "RPC (external)",
                             pdf_type, -1, "none",
                             shape_sets, control_region_sets, "primary", hist_file);
}

#endif
