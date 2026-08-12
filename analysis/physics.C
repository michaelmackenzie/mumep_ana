// Relevant physics constants
#ifndef __CONVANA_ANALYSIS_PHYSICS__
#define __CONVANA_ANALYSIS_PHYSICS__

//--------------------------------------------------------
// Livetime/muon stop normalization

// General configurations

double duty_cycle_1bb_         = 0.322; // taken from https://github.com/Mu2e/Production/blob/main/JobConfig/ensemble/python/normalizations.py
double duty_cycle_2bb_         = 0.246;
double npot_event_1bb_         = 1.58e7; // average values
double npot_event_2bb_         = 3.93e7;
double npot_rate_1bb_          = npot_event_1bb_ / 1.695e-6; // instantaneous N(POT) per second
double npot_rate_2bb_          = npot_event_2bb_ / 1.695e-6;
double walltime_month_         = 2.6166667e6; // appox 1/12th of a year * duty factor: 3.14e7/12 (1BB)
double walltime_week_          = 603420.0; // 7x24x60x60
double livetime_month_1bb_     = duty_cycle_1bb_ * walltime_month_;
double livetime_month_2bb_     = duty_cycle_2bb_ * walltime_month_;
double livetime_week_1bb_      = duty_cycle_1bb_ * walltime_week_ ;
double livetime_week_2bb_      = duty_cycle_2bb_ * walltime_week_ ;
double npot_month_1bb_         = npot_rate_1bb_ * livetime_month_1bb_;
double npot_month_2bb_         = npot_rate_2bb_ * livetime_month_2bb_;

// SU2020 info
double livetime_su2020_        = 1.11e7;
double npot_su2020_            = 3.8e19;
double nmuons_per_pot_su2020_  = 0.00159;
double nmuons_su2020_          = npot_su2020_*nmuons_per_pot_su2020_;

// MDSd info
double livetime_mdsd_        = 9.96697e+06;
double npot_mdsd_            = 2.997584962406015e+19;
double nmuons_per_pot_mdsd_  = 0.0015576895;
double nmuons_mdsd_          = npot_mdsd_*nmuons_per_pot_mdsd_;

// MDSf info
double livetime_mdsf_        = 623819.;
double npot_mdsf_            = 1.8761473684210527e+18;
double nmuons_per_pot_mdsf_  = 0.0015576895;
double nmuons_mdsf_          = npot_mdsd_*nmuons_per_pot_mdsd_;

// MDSg info
double livetime_mdsg_        = 623819.;
double npot_mdsg_            = 1.8761473684210527e+18;
double nmuons_per_pot_mdsg_  = 0.0015576895;
double nmuons_mdsg_          = npot_mdsd_*nmuons_per_pot_mdsd_;

// MDS2a info
double livetime_mds2a_       = 1.06939e+07;
double npot_mds2a_           = 3.216210526315789e+19;
double nmuons_per_pot_mds2a_ = 0.0015576895;
double nmuons_mds2a_         = npot_mds2a_*nmuons_per_pot_mds2a_;

// MDS2b info
double livetime_mds2b_       = 1.06939e+07;
double npot_mds2b_           = 3.216210526315789e+19;
double nmuons_per_pot_mds2b_ = 0.0015576895;
double nmuons_mds2b_         = npot_mds2a_*nmuons_per_pot_mds2a_;

// MDS2c info
double livetime_mds2c_       = 1.06939e+07 * (812./822.);
double npot_mds2c_           = 3.216210526315789e+19 * (812./822.);
double nmuons_per_pot_mds2c_ = 0.0015576895;
double nmuons_mds2c_         = npot_mds2a_*nmuons_per_pot_mds2a_;

// Run 1a info
double livetime_run1a_       = 28.*24.*60.*60.*duty_cycle_1bb_; // 28 days of wall time, 100% up
double npot_run1a_           = livetime_run1a_ * npot_rate_1bb_; // livetime * (POT/event)/1695.e-9
double nmuons_per_pot_run1a_ = 0.000767114;
double nmuons_run1a_         = npot_run1a_*nmuons_per_pot_run1a_;

// Normalization used, defaulting to Run 1A
double livetime_       = livetime_run1a_       ; // default to SU2020
double npot_           = npot_run1a_           ;
double nmuons_         = nmuons_run1a_         ;
double nmuons_per_pot_ = nmuons_per_pot_run1a_ ;
double signal_br_      = 1.e-15                ; // example signal branching fraction for normalization

//--------------------------------------------------------
// Branching fractions

// MDS1g:
// pistops/POT 0.00223279
// pi time eff 0.06459918608102208
// pi surv prob 0.0008919234313257805
// BRRPC= 0.0215

double muon_capture_fraction_ = 0.609    ; // on aluminum
double muon_capture_carbon_   = 0.0701   ;
double ipa_nmuons_per_pot_    = 2.062e-08; // stopped muon / POT in IPA
double ipa_dio_frac_70_       = 2.538e-06; //fraction of carbon DIO spectrum above 70 MeV
double pion_stop_rate_        = 0.0018801*0.51656; //Pion stops in the target (infinite lifetime)
double pion_survive_frac_     = 2393.60487 / 1e10; // sum of sampled weights / N(sampled pions)
double rpc_frac_50_           = 0.9888   ; //fraction of RPC > 50 MeV, cutoff used in simulation
double rpc_br_                = 0.0215   ;
double rpc_int_br_            = 0.0069   ; //using BR(internal RPC) / BR(RPC) on hydrogen
double pbar_stops_per_pot_    = 4.7e-18  ; //N(pbar at ST) / POT
double dio_frac_95_           = 3.637e-11; //DIO spectrum fraction above 95 MeV

// RMC data
double br_rmc_                = 1.40e-5  ; //For E > 57 MeV, relative to muon capture
double rmc_frac_57_           = 0.14800  ; //fraction of RMC > 57 (kmax = 90.1)
double rmc_frac_85_           = 0.0010641; //fraction of RMC > 85 (kmax = 90.1)
double rmc_int_br_            = 0.0069   ; //using BR(internal RPC) / BR(RPC) on hydrogen
double rmc_conv_              = 380711./1.e8; // fraction of RMC photons that convert

// RMC phase-space model results
double rmc_ps_0n_57_r_        = 0.099    ; // R(0 knockout | E > 57) / R(RMC | E > 57)
double rmc_ps_1n_57_r_        = 0.901    ; // R(1 knockout | E > 57) / R(RMC | E > 57)
double rmc_ps_0n_frac_57_     = 0.22887  ; // Fraction above 57 MeV
double rmc_ps_1n_frac_57_     = 0.061620 ; // Fraction above 57 MeV
double rmc_ps_0n_frac_80_     = 0.03319  ; // Fraction above 80 MeV
double rmc_ps_1n_frac_80_     = 0.0013175; // Fraction above 80 MeV


void init_physics(TString tag) {
  tag.ToLower();

  // livetime_       = livetime_su2020_      ; // default to SU2020
  // npot_           = npot_su2020_          ;
  // nmuons_         = nmuons_su2020_        ;
  // nmuons_per_pot_ = nmuons_per_pot_su2020_;

  livetime_       = livetime_run1a_        ; // default to Run 1A
  npot_           = npot_run1a_            ;
  nmuons_per_pot_ = nmuons_per_pot_run1a_  ;
  nmuons_         = nmuons_per_pot_ * npot_;

  if(tag.Contains("mds1d")) {
    const float ad_hoc(1.f);
    livetime_       = duty_cycle_1bb_*livetime_mdsd_ ; // MDS1d
    npot_           = ad_hoc*npot_mdsd_          ;
    nmuons_         = ad_hoc*nmuons_mdsd_        ;
    nmuons_per_pot_ = nmuons_per_pot_mdsd_       ;
    signal_br_      = 1.e-13                     ;
  }

  if(tag.Contains("mds1f")) {
    const float ad_hoc(0.9f); // about 10% of files were lost when producing MDS1f
    livetime_       = duty_cycle_1bb_*livetime_mdsf_ ; // MDS1f
    npot_           = ad_hoc*npot_mdsf_          ;
    nmuons_         = ad_hoc*nmuons_mdsf_        ;
    nmuons_per_pot_ = nmuons_per_pot_mdsf_       ;
    signal_br_      = 1.e-13                     ;
  }

  if(tag.Contains("mds1g")) {
    const float ad_hoc(1.f);
    livetime_       = duty_cycle_1bb_*livetime_mdsg_ ; // MDS1g
    npot_           = ad_hoc*npot_mdsg_          ;
    nmuons_         = ad_hoc*nmuons_mdsg_        ;
    nmuons_per_pot_ = nmuons_per_pot_mdsg_       ;
    signal_br_      = 1.e-13                     ;
  }

  if(tag.Contains("mds2a")) {
    const float ad_hoc(1.f);
    livetime_       = livetime_mds2a_             ; // MDS2a
    npot_           = ad_hoc*npot_mds2a_          ;
    nmuons_         = ad_hoc*nmuons_mds2a_        ;
    nmuons_per_pot_ = nmuons_per_pot_mds2a_       ;
    signal_br_      = 1.e-13                      ;
  }

  if(tag.Contains("mds2b")) {
    const float ad_hoc(1.f);
    livetime_       = livetime_mds2b_             ; // MDS2b
    npot_           = ad_hoc*npot_mds2b_          ;
    nmuons_         = ad_hoc*nmuons_mds2b_        ;
    nmuons_per_pot_ = nmuons_per_pot_mds2b_       ;
    signal_br_      = 1.e-13                      ;
  }

  if(tag.Contains("mds2c")) {
    const float ad_hoc(1.f);
    livetime_       = livetime_mds2c_             ; // MDS2c
    npot_           = ad_hoc*npot_mds2c_          ;
    nmuons_         = ad_hoc*nmuons_mds2c_        ;
    nmuons_per_pot_ = nmuons_per_pot_mds2c_       ;
    signal_br_      = 1.e-13                      ;
  }

  if(tag.Contains("mds3c") && !use_evtana_) { // MDS3c am had fewer events than MDS3c ar
    const float ad_hoc(1.f); // 0.55f
    livetime_       = 4401645.8                      ; // MDS3c
    npot_           = ad_hoc*livetime_*npot_rate_1bb_;
    nmuons_         = npot_*nmuons_per_pot_run1a_    ;
    nmuons_per_pot_ = nmuons_per_pot_run1a_          ;
    signal_br_      = 1.e-13                         ;
  }

  if(tag.Contains("mds3c") && use_evtana_) { // MDS3c am had fewer events than MDS3c ar
    const float ad_hoc(1.f);
    livetime_       = 4401645.8*4921433./4863968.    ; // MDS3c
    npot_           = ad_hoc*livetime_*npot_rate_1bb_;
    nmuons_         = npot_*nmuons_per_pot_run1a_    ;
    nmuons_per_pot_ = nmuons_per_pot_run1a_          ;
    signal_br_      = 1.e-13                         ;
  }

  if(tag.Contains("mdc2025") || tag.Contains("run1a")) {
    const float ad_hoc(1.f);
    livetime_       = livetime_run1a_             ; // Run 1a
    npot_           = ad_hoc*npot_run1a_          ;
    nmuons_         = ad_hoc*nmuons_run1a_        ;
    nmuons_per_pot_ = nmuons_per_pot_run1a_       ;
    signal_br_      = 1.e-13                      ;
  }
}

#endif
