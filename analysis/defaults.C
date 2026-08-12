#ifndef __CONVANA_ANALYSIS_DEFAULTS__
#define __CONVANA_ANALYSIS_DEFAULTS__

const char*  hist_path_   = "/exp/mu2e/data/projects/run1a/mumep_ana/histograms/";
bool         hist_pdfs_   = true ; // Use functions or histograms in the model
bool         include_sys_ = false; // Evaluate systematics
bool         use_evtana_  = true ; // Use Mu2eEvtAna inputs (EventNtuple)
bool         run1a_range_ = false; // Use Run 1A paper range (only in EvtAna mode)

// Histogram file info
int          hist_mode_   = 2;
TString      dir_path_    = "Ana/ConvAna_ConvAna/"; // path in histogram files
TString      file_type_   = "hist"; // histogram file extension

// Fit info
TString      var_         = "obs"; // Observable name
double       bin_width_   = 0.25;  // Expected bin width, rebin to achieve if possible
double       xmin_em_     = 100.;  // Momentum range for mu- --> e- fit
double       xmax_em_     = 110.;
double       xmin_ep_     =  87.;  // Momentum range for mu- --> e+ fit
double       xmax_ep_     =  97.;

void set_evtana_defaults() {
  dir_path_  = "Ana/";
  file_type_ = "root";
  hist_mode_ = 1;
  include_sys_ = false; // systematics aren't fully implemented
  if(run1a_range_) {
    // xmin_em_ = 97.;
    bin_width_ = 0.25; // to better match the paper fit
    // bin_width_ = 0.1; // to better match the cut-and-count range
  }
}
#endif
