// Fit the antiproton spectrum shape
#ifndef __CONVANA_ANALYSIS_PBARFIT__
#define __CONVANA_ANALYSIS_PBARFIT__

#include "component_fit.C"

int pbar_fit(TString process = "mumem", int selection = 20, TString tag = "", TString pdf_type = "default",
             std::vector<int> shape_sets = {}, std::vector<int> control_region_sets = {}) {
    if(use_evtana_) {
        cout << __func__ << ": skipping pbar fit when use_evtana_ is true" << endl;
        return 0;
    }
  return fit_component_model(process, selection, tag,
                             "pbar", "Antiproton",
                             pdf_type, -1, "none",
                             shape_sets, control_region_sets, "primary");
}

#endif
