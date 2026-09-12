// Fit the RMC internal spectrum shape
#ifndef __CONVANA_ANALYSIS_RMCINTFIT__
#define __CONVANA_ANALYSIS_RMCINTFIT__

#include "component_fit.C"

// knockout: "" for the single (0n) component, or "0n"/"1n" for the knockout-split components
int rmc_int_fit(TString process = "mumem", int selection = 20, TString tag = "", TString pdf_type = "default",
                TString tail_model = "default", std::vector<int> shape_sets = {},
                std::vector<int> control_region_sets = {}, TString knockout = "") {
  const TString component = rmc_component_name("rmc_int", knockout);
  TString title; int color;
  set_style(component, title, color);
  return fit_component_model(process, selection, tag,
                             component, title,
                             pdf_type, -1, tail_model,
                             shape_sets, control_region_sets, "primary");
}

#endif
