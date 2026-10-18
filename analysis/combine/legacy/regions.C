// Fit regions: a selection set is either a signal region or a control region of a base set (< 1000).
// Control regions add flag offsets to the base set, which combine, e.g. 75 + 1000 + 2000 = 3075.
// The offsets follow the histogram filling in Mu2eEvtAna (src/ConvAna.cc, inc/GlobalConstants.h):
//   + 1000 (kCRVVetoOffset): inverted CRV veto, the track has a coincident CRV cluster
//   + 2000 (kTimeCutOffset): inverted time cut. For set 75 this is the early time region,
//                            t0 <= 540 ns; for sets 20/40 it is t0 outside 475-1650 ns
//   + 4000 (kChargeOffset) : charge-inverted tracks
#ifndef __CONVANA_ANALYSIS_REGIONS__
#define __CONVANA_ANALYSIS_REGIONS__

#include "../../defaults.C"

// Control region flag offsets. perform_fits.sh and full_loop.sh use the same convention.
const int crv_inverted_offset_    = 1000;
const int time_inverted_offset_   = 2000;
const int charge_inverted_offset_ = 4000;

// Whether a selection set has the given control region flag
bool has_region_flag(const int selection, const int offset) {
  return ((selection / 1000) & (offset / 1000)) != 0;
}

struct region_t {
  int    base_set_        = 0;
  bool   crv_inverted_    = false; // tracks are tagged by the CRV
  bool   time_inverted_   = false; // tracks fail the t0 cut (early time for set 75)
  bool   charge_inverted_ = false; // tracks have the opposite charge of the signal
  bool   negative_tracks_ = true ; // the region selects e- tracks
  bool   has_signal_      = true ; // the process's signal is modeled in the region
  double xmin_            = 0.   ; // momentum range of the fit
  double xmax_            = 0.   ;
  bool is_control() const { return crv_inverted_ || time_inverted_ || charge_inverted_; }
};

region_t get_region(TString process, const int selection) {
  process.ToLower();
  const bool is_mumem = process == "mumem";
  region_t region;
  region.base_set_        = selection % 1000;
  region.crv_inverted_    = has_region_flag(selection, crv_inverted_offset_);
  region.time_inverted_   = has_region_flag(selection, time_inverted_offset_);
  region.charge_inverted_ = has_region_flag(selection, charge_inverted_offset_);
  region.negative_tracks_ = is_mumem != region.charge_inverted_;
  // The signal has ~no tracks of the opposite charge, but does leak into the CRV-tagged and
  // time-inverted regions (~0.1% and ~17% of the signal region yield for mu- --> e- set 75)
  region.has_signal_      = !region.charge_inverted_;
  if(region.charge_inverted_) {
    region.xmin_ = (is_mumem) ? xmin_ci_em_ : xmin_ci_ep_;
    region.xmax_ = (is_mumem) ? xmax_ci_em_ : xmax_ci_ep_;
  } else {
    region.xmin_ = (is_mumem) ? xmin_em_ : xmin_ep_;
    region.xmax_ = (is_mumem) ? xmax_em_ : xmax_ep_;
  }
  return region;
}

#endif
