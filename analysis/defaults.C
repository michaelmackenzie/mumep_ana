#ifndef __CONVANA_ANALYSIS_DEFAULTS__
#define __CONVANA_ANALYSIS_DEFAULTS__

// Histogram directory, overridden by the CONVANA_HIST_PATH environment variable (e.g. for make_selection_hists.C outputs)
TString      hist_path_str_ = (gSystem->Getenv("CONVANA_HIST_PATH")) ? TString(gSystem->Getenv("CONVANA_HIST_PATH")) + "/"
                                                                       : TString("/exp/mu2e/data/projects/run1a/mumep_ana/histograms/");
const char*  hist_path_   = hist_path_str_.Data();

// Boolean flag from the environment (1/true/yes or 0/false/no), or the given default if it is not set
bool env_flag(const char* name, const bool default_value) {
  const char* value = gSystem->Getenv(name);
  if(!value) return default_value;
  TString v(value); v.ToLower();
  if(v == "1" || v == "true"  || v == "yes") return true;
  if(v == "0" || v == "false" || v == "no" ) return false;
  printf("env_flag: Unrecognized value %s=%s, using the default (%i)\n", name, value, default_value);
  return default_value;
}

bool         hist_pdfs_   = env_flag("CONVANA_HIST_PDFS", true ); // Use functions or histograms in the model
bool         include_sys_ = false; // Evaluate systematics
bool         use_evtana_  = true ; // Use Mu2eEvtAna inputs (EventNtuple)
bool         use_2bb_     = env_flag("CONVANA_USE_2BB", false); // Use the Mu2eEvtAna inputs re-weighted to 2BB (mode 2), with the 2BB exposure (only in EvtAna mode)
bool         run1a_range_ = false; // Use Run 1A paper range (only in EvtAna mode)
bool         include_t0_  = true ; // Add t0 histograms to the output
bool         do_2d_fit_   = false; // Perform a (p, t) fit instead of a 1D p fit
bool         unbinned_    = false; // Perform unbinned fits
bool         use_env_     = env_flag("CONVANA_USE_ENV", false); // Model the mu- --> e+ background with an envelope
bool         blind_data_  = true ; // Blind the signal region when evaluating envelope fits
int          verbose_     = 1    ; // Printout level for the envelope construction

// Histogram file info
int          hist_mode_   = 2;
TString      dir_path_    = "Ana/ConvAna_ConvAna/"; // path in histogram files
TString      file_type_   = "hist"; // histogram file extension

// Fit info
TString      var_         = "obs"; // Observable name
double       bin_width_   =  0.2 ; // Expected bin width, rebin to achieve if possible
double       t_bin_width_ =   25.; // Time distribution bin width
double       xmin_em_     =  100.; // Momentum range for mu- --> e- fit
double       xmax_em_     =  110.;
double       xmin_ep_     =   87.; // Momentum range for mu- --> e+ fit
double       xmax_ep_     =   97.;
double       xmin_ci_em_  =   95.; // Momentum range for the charge-inverted (e+ track) control region of mu- --> e-
double       xmax_ci_em_  =  105.;
double       xmin_ci_ep_  =   87.; // Momentum range for the charge-inverted (e- track) control region of mu- --> e+
double       xmax_ci_ep_  =   97.;
double       t_min_       =  475.; // Time range for 2D fits
double       t_max_       = 1650.;

void set_evtana_defaults() {
  dir_path_  = "Ana/";
  file_type_ = "root";
  hist_mode_ = (use_2bb_) ? 2 : 1; // mode 2: beam intensity re-weighted to 2BB
  include_sys_ = false; // systematics aren't fully implemented
  if(run1a_range_) {
    // xmin_em_ = 97.;
    bin_width_ = 0.2; // to better match the paper fit
  }
}
#endif
