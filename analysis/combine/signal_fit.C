// Fit the DIO spectrum shape
#ifndef __CONVANA_ANALYSIS_SIGNALFIT__
#define __CONVANA_ANALYSIS_SIGNALFIT__

#include "component_fit.C"

//---------------------------------------------------------------------------------------------------------------------------
int signal_fit(TString process = "mumem", int selection = 20, TString tag = "", const int isys = -1,
               TString pdf_type = "default", TString tail_model = "default", std::vector<int> shape_sets = {}) {
  return fit_component_model(process, selection, tag,
                             "signal", "Signal",
                             pdf_type, isys, tail_model,
                             shape_sets, {}, "primary");
}

#endif
