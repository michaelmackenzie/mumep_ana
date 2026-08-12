#ifndef __CONVANA_ANALYSIS_SYSTEMATICS__
#define __CONVANA_ANALYSIS_SYSTEMATICS__
#include "../../stn/ana/ana/Systematics.hh"
#include "dio_fit.C"
#include "cosmic_fit.C"
#include "signal_fit.C"

// Systematic information
ConvAna::Systematics fSystematics;

int systematics(TString process = "mumem", int selection = 20, TString tag = "", TString sys = "") {
  if(use_evtana_) set_evtana_defaults();
  std::vector<TString> names = { // systematics to consider
    "Scale"
  };

  int status(0);
  for(auto name : names) {
    if(sys != "" && sys != name) continue;
    const int up = fSystematics.GetNum(name);
    cout << __func__ << ": Performing DIO fit for systematic " << name.Data() << " up (" << up << ")\n";
    status += dio_fit(process, selection, tag, up);
    if(fSystematics.GetName(up+1) != name) continue;
    cout << __func__ << ": Performing DIO fit for systematic " << name.Data() << " down (" << up+1 << ")\n";
    status += dio_fit(process, selection, tag, up+1);
    cout << __func__ << ": Performing Cosmic fit for systematic " << name.Data() << " up (" << up << ")\n";
    status += cosmic_fit(process, selection, tag, up);
    if(fSystematics.GetName(up+1) != name) continue;
    cout << __func__ << ": Performing Cosmic fit for systematic " << name.Data() << " down (" << up+1 << ")\n";
    status += cosmic_fit(process, selection, tag, up+1);
    cout << __func__ << ": Performing Signal fit for systematic " << name.Data() << " up (" << up << ")\n";
    status += signal_fit(process, selection, tag, up);
    if(fSystematics.GetName(up+1) != name) continue;
    cout << __func__ << ": Performing Signal fit for systematic " << name.Data() << " down (" << up+1 << ")\n";
    status += signal_fit(process, selection, tag, up+1);
  }

  return status;
}

#endif
