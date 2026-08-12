#ifndef __CONVANA_ANALYSIS_DATASETINFO__
#define __CONVANA_ANALYSIS_DATASETINFO__
#include "defaults.C"
#include "tools/types.C"
#include "physics.C"
#include <map>

std::map<TString, DatasetInfo_t> datasets_;
const char* hist_func_    = "cnv_ana"; // for histogram file naming
TString     rmc_spectrum_ = "s0v0"; // RMC spectrum to use
bool        combine_rpc_  = true  ; // Combine RPC samples into one label
bool        combine_rmc_  = false ; // Combine RMC samples into one label

// Background model options
bool        physical_dio_     = false ; // Use flat electron or physical sample
bool        physical_rpc_ext_ = true  ; // Use infinite lifetime or physical sample
bool        physical_rpc_int_ = true  ; // Use infinite lifetime or physical sample
bool        physical_rmc_ext_ = false ; // Use flat photon or physical sample
bool        physical_rmc_int_ = false ; // Use flat photon or physical sample

void init_dataset_info() {
  if(datasets_.size() > 0) return;
  const double rate_dio((1.-muon_capture_fraction_)*nmuons_per_pot_), rate_dio_95(rate_dio*dio_frac_95_),
    rate_rpc(pion_stop_rate_*rpc_br_), rate_rpc_int(rate_rpc*rpc_int_br_),
    rate_phys_rpc(rate_rpc*pion_survive_frac_), rate_phys_rpc_int(rate_phys_rpc*rpc_int_br_),
    rate_rmc(nmuons_per_pot_*muon_capture_fraction_*br_rmc_), rate_rmc_conv(rate_rmc*rmc_conv_),
    rate_rmc_int(rate_rmc*rmc_int_br_); // FIXME: rmc_int_br in convolution sometimes
  const double rate_rmc_c(rate_rmc * (0.999 - 0.992)/2.); // calo aimed RMC photons
  const double rate_ipa_dio((1.-muon_capture_carbon_)*ipa_nmuons_per_pot_*ipa_dio_frac_70_);
  const double rate_rmc_85(rate_rmc*rmc_frac_85_), rate_rmc_85_int(rate_rmc_85*rmc_int_br_);
  const double rate_pbar(pbar_stops_per_pot_);
  const double rate_sig(muon_capture_fraction_*nmuons_per_pot_); // R_mue is relative to BR(muon capture)
  const double rate_b1(1./1.6e7); // 1 / N(POT per event);

  // RMC phase-space rates
  const double rate_rmc_ps_0n_80 = rate_rmc * (rmc_ps_0n_57_r_ / rmc_ps_0n_frac_57_) * rmc_ps_0n_frac_80_;
  const double rate_rmc_ps_1n_80 = rate_rmc * (rmc_ps_1n_57_r_ / rmc_ps_1n_frac_57_) * rmc_ps_1n_frac_80_;

  // Retrieve N(gen) using scripts/samCountGenEvents.sh and N(events) using scripts/samCountEvents.sh
  //                                  N(gen)      N(digi) emin emax    rate          stn dataset           digi dataset

  // const double livetime_digi = (49535.e5 / 2585823777.) * 556000.; // N(gen digi) / N(gen sim) * livetime (sim) (CosmicAll)
  const double t_cry4ab1 = 4437500.0; // 1.37e7 // livetime
  // mcs evaluation: 4437713

  if(!use_evtana_) { // Stntuple inputs
    // MDC2025
    datasets_["dio-flat"]    = DatasetInfo_t(7.5e6     , 2495714,70., 120., rate_dio         , "fele0b1s5r0102", "dig.mu2e.FlateMinusMix1BBTriggerable.MDC2025af_best_v1_1.art"             );
    datasets_["dio-phys"]    = DatasetInfo_t(2.5e7     , 9366309, 0.,   1., rate_dio_95      , "dio00b0s5r0102", "dig.mu2e.DIOtail95OnSpill.MDC2025ap_best_v1_1.art"                        );
    datasets_["cosmic"]      = DatasetInfo_t(t_cry4ab1 , 4454595, 0.,   1.,       1.         , "cry4ab1s5r0102", "dig.mu2e.CosmicSignalMix1BBTriggerable.MDC2025af_best_v1_1.art"           );
    datasets_["cosmic_all"]  = DatasetInfo_t(1065094.  ,136334666,0.,   1.,       1.         , "cry4bb1s5r0101", "dig.mu2e.CosmicAllMix1BBTriggerable.MDC2025af_best_v1_1.art"              );
    // datasets_["rmc_int-flat"]= DatasetInfo_t(1e6       ,  440372,90., 110., rate_rmc_int     , "fpos0b1s5r0102", "dig.mu2e.FlatePlusMix1BBTriggerable.MDC2025af_best_v1_1.art"              );
    datasets_["rmc_int-flat"]= DatasetInfo_t(7500000   , 2282799,70., 120., rate_rmc_int     , "fpos0b0s5r0102", "dig.mu2e.FlatePlusOnSpill.MDC2025ap_best_v1_1.art"                        );
    datasets_["mnbs_b1"]     = DatasetInfo_t(4e6       ,     4e6, 0.,   1.,   rate_b1        , "mnbs0b1s5r0102", "dig.mu2e.NoPrimaryMix1BBTriggerable.MDC2025af_best_v1_1.art"              );
    datasets_["rmc_ext-flat"]= DatasetInfo_t(1.e9      ,  555992,70., 102., rate_rmc         , "fgam0b1s5r0102", "dig.mu2e.FlatGammaMix1BBTriggerable.MDC2025af_best_v1_1.art"              );
    datasets_["rmc_ext_c"]   = DatasetInfo_t(15.e6     , 5978198,70., 110., rate_rmc_c       , "fgamcb1s5r0101.s1v0", "dig.mu2e.FlatGammaCaloMix1BBTriggerable.MDC2025af_best_v1_1.art"     );
    datasets_["rpc_ext-flat"]= DatasetInfo_t(5e9       , 2982174, 0.,   1., rate_rpc         , "rpce0b0s5r0102", "dig.mu2e.RPCExternallOnSpill.MDC2025ap_best_v1_1.art"                     ); // N(gen) is true N(gen), not N(gen count) due to filter before it
    datasets_["rpc_ext-phys"]= DatasetInfo_t(5000000000,  344658, 0.,   1., rate_phys_rpc    , "rpce1b0s5r0102", "dig.mu2e.RPCExternalPhysicalOnSpill.MDC2025an_best_v1_1.art"              );
    datasets_["rpc_int-phys"]= DatasetInfo_t(125000000 , 1730738, 0.,   1., rate_phys_rpc_int, "rpci1b0s5r0102", "dig.mu2e.RPCInternalPhysicalOnSpill.MDC2025an_best_v1_1.art"              );

    // datasets_["mumem"]       = DatasetInfo_t(9995000  , 4885425, 0.,   1., rate_sig       , "cele0b1s5r0102", "dig.mu2e.CeEndpointMix1BBTriggerable.MDC2025af_best_v1_1.art"             );
    // datasets_["mumep"]       = DatasetInfo_t(9998000  , 3553349, 0.,   1., rate_sig       , "cpos0b1s5r0102", "dig.mu2e.CePlusEndpointMix1BBTriggerable.MDC2025af_best_v1_1.art"         );
    datasets_["mumem"]       = DatasetInfo_t(1e7      , 4756424, 0.,   1., rate_sig       , "cele1b1s5r0102", "dig.mu2e.CeMLEadingLogMix1BBTriggerable.MDC2025af_best_v1_1.art"             );
    // FIXME: Used LO for LL generation
    datasets_["mumep"]       = DatasetInfo_t(9960000  , 3546334, 0.,   1., rate_sig       , "cpos1b1s5r0102", "dig.mu2e.CePlusEndpointMix1BBTriggerable.MDC2025af_best_v1_1.art"         );

    // Assign input options for some datasets
    datasets_["dio"]     = (physical_dio_)     ? datasets_["dio-phys"]     : datasets_["dio-flat"];
    datasets_["rpc_ext"] = (physical_rpc_ext_) ? datasets_["rpc_ext-phys"] : datasets_["rpc_ext-flat"];
    datasets_["rpc_int"] = (physical_rpc_int_) ? datasets_["rpc_int-phys"] : datasets_["rpc_int-flat"];
    datasets_["rmc_ext"] = (physical_rmc_ext_) ? datasets_["rmc_ext-phys"] : datasets_["rmc_ext-flat"];
    datasets_["rmc_int"] = (physical_rmc_int_) ? datasets_["rmc_int-phys"] : datasets_["rmc_int-flat"];

    // MDC2020

    // Background
    // datasets_["dio"]     = DatasetInfo_t(1e7      ,  678584, 0., 105., rate_dio       , "fele0b1s5r0001", "dig.mu2e.FlateMinusMix1BBTriggerable.MDC2020aj_best_v1_3.art"             );
    // datasets_["dio"]     = DatasetInfo_t(751400   ,  241850, 0.,   1., rate_dio_95    , "dio00b1s5r0001", "dig.mu2e.DIOtail_95Mix1BBTriggerable.MDC2020ar_best_v1_3.art"             );
    // datasets_["cosmic"]  = DatasetInfo_t(1.1e7    , 4950180, 0.,   1.,       1.       , "cry4ab1s5r0001", "dig.mu2e.CosmicCRYSignalAllMix1BBTriggerable.MDC2020ap_best_v1_3.art"     );
    // datasets_["rpc_ext"] = DatasetInfo_t(394706.e4,  280484, 0.,   1., rate_rpc       , "rpce0b0s5r0000", "dig.mu2e.RPCExternalOnSpillTriggerable.MDC2020aq_best_v1_3.art"           );
    // datasets_["rpc_int"] = DatasetInfo_t(35196.e3 ,  255863, 0.,   1., rate_rpc_int   , "rpci0b0s5r0000", "dig.mu2e.RPCInternalOnSpillTriggerable.MDC2020aq_best_v1_3.art"           );
    // datasets_["rpc_ext"] = DatasetInfo_t(1998000000, 811956, 0.,   1., rate_rpc       , "rpce1b0s5r0001", "dig.mu2e.RPCExternalOnSpillTriggerable.MDC2020az_best_v1_3.art"           );
    // datasets_["rpc_int"] = DatasetInfo_t(100000000, 7366447, 0.,   1., rate_rpc_int   , "rpci1b0s5r0001", "dig.mu2e.RPCInternalOnSpillTriggerable.MDC2020az_best_v1_3.art"           );
    datasets_["pbar"]    = DatasetInfo_t(1e6      ,  172747, 0.,   1., rate_pbar      , "pbar1b2s5r0000", "dig.mu2e.PbarSTGunMix2BBTriggerable.MDC2020ar_best_v1_3.art"              );
    // datasets_["rmc_int"] = DatasetInfo_t(7416000  ,  502958, 0., 105., rate_rmc_int   , "fpos0b1s5r0000", "dig.mu2e.FlatePlusMix1BBTriggerable.MDC2020aj_best_v1_3.art"              );
    // datasets_["rmc_ext"] = DatasetInfo_t(99.6e6   , 5831595,80., 102.,rate_rmc_conv   , "rmcg0b0s6r0000", "dig.mu2e.RMCFlatGammaResamplingOnSpillTriggerable.MDC2020ar_best_v1_3.art");
    // datasets_["rmc_int"] = DatasetInfo_t(5e8      , 2359836, 0.,   1., rate_rmc_85_int, "rmci0b0s5r0000", "dig.mu2e.RMCInternalOnSpillTriggerable.MDC2020au_best_v1_3.art"           );
    // datasets_["rmc_ext"] = DatasetInfo_t(4e9      ,  403007, 0.,   1., rate_rmc_85    , "rmce0b0s5r0000", "dig.mu2e.RMCExternalCatOnSpillTriggerable.MDC2020au_best_v1_3.art"        );
    datasets_["ipa_dio"] = DatasetInfo_t(1e7      ,  254647, 0.,   1., rate_ipa_dio   , "ipad0b0s5r0001", "dig.mu2e.IPAMuminusMichelOnSpillTriggerable.MDC2020au_best_v1_3.art"      );

    // Signal
    // datasets_["mumem"]   = DatasetInfo_t(4000000  , 1495720, 0.,   1., rate_sig       , "cele1b1s5r0001", "dig.mu2e.CeMLeadingLogMix1BBTriggerable.MDC2020am_best_v1_3.art"          );
    // datasets_["mumep"]   = DatasetInfo_t(4000000  , 1117816, 0.,   1., rate_sig       , "cpos0b1s5r0000", "dig.mu2e.CePlusEndpointMix1BBTriggerable.MDC2020ae_best_v1_3.art"         );
    // Leading order signal
    // datasets_["mumem"]   = DatasetInfo_t(4000000  , 1541335, 0.,   1., rate_sig       , "cele0b1s5r0000", "dig.mu2e.CeEndpointMix1BBTriggerable.MDC2020ae_best_v1_3.art"             );


    // Data
    datasets_["data_mds1d"] = DatasetInfo_t(1.    , 1705987, 0.,   1., 1., "mds1db0s5r0000", "dig.mu2e.ensembleMDS1dOnSpillTriggerable.MDC2020aq_best_v1_3.art"    );
    datasets_["data_mds1f"] = DatasetInfo_t(1.    ,   95961, 0.,   1., 1., "mds1fb1s5r0000", "dig.mu2e.ensembleMDS1fMix1BBTriggerable.MDC2020ai_perfect_v1_3.art"  );
    datasets_["data_mds1g"] = DatasetInfo_t(1.    ,   69803, 0.,   1., 1., "mds1gb0s5r0000", "dig.mu2e.ensembleMDS1gOnSpillTriggerable.MDC2020aq_best_v1_3.art"    );
    datasets_["data_mds2a"] = DatasetInfo_t(1.    , 2617404, 0.,   1., 1., "mds2ab0s6r0000", "mcs.sophie.MDS2aTriggered.v0.art"                                    );
    datasets_["data_mds2b"] = DatasetInfo_t(1.    , 5687636, 0.,   1., 1., "mds2bb0s5r0001", "dig.mu2e.ensembleMDS2bOnSpillTriggerable.MDC2020az_best_v1_3.art"    );
    datasets_["data_mds2c"] = DatasetInfo_t(1.    , 5597655, 0.,   1., 1., "mds2cb1s5r0001", "dig.mu2e.ensembleMDS2cMix1BBTriggerable.MDC2020ba_best_v1_3.art"     );
    datasets_["data_mds3c"] = DatasetInfo_t(1.    , 4863968, 0.,   1., 1., "mds3cb1s5r0102", "dig.mu2e.ensembleMDS3cMix1BB.MDC2025am_best_v1_1.art"                );
  } else { // EventNtuple inputs
    // datasets_["mumem"]   = DatasetInfo_t(1e7       , 4009075, 0.,   1., rate_sig                     , "cele1b0s5r0100", "nts.mu2e.CeMLeadingLogOnSpill-reco-ntuple.MDC2025-002.root"  );
    // datasets_["cosmic"]  = DatasetInfo_t(t_cry4ab1 , 4120241, 0.,   1.,       1.                     , "cry4ab0s5r0100", "nts.mu2e.CosmicSignalOnSpill-reco-ntuple.MDC2025-002.root"   );
    // datasets_["dio"]     = DatasetInfo_t(25e6      , 8780533, 0.,   1., rate_dio_95                  , "dio00b0s5r0100", "nts.mu2e.DIOtail95OnSpill-reco-ntuple.MDC2025-002.root"      );
    datasets_["mumem"]      = DatasetInfo_t(1e7       , 4137592, 0.,   1., rate_sig                     , "cele1b1s5r0100", "nts.mu2e.CeMLeadingLogMix1BB.MDC2025ar_best_v1_1.root"  );
    datasets_["cosmic"]     = DatasetInfo_t(t_cry4ab1 , 4152355, 0.,   1.,       1.                     , "cry4ab1s5r0100", "nts.mu2e.CosmicSignalMix1BB.MDC2025ar_best_v1_1.root"   );
    datasets_["dio"]        = DatasetInfo_t(25e6      , 9348862, 0.,   1., rate_dio_95                  , "dio00b1s5r0100", "nts.mu2e.DIOtail95Mix1BB.MDC2025ar_best_v1_1.root"      );
    datasets_["rmc_ext_0n"] = DatasetInfo_t(7000000000, 4952890, 0.,   1., rate_rmc_ps_0n_80            , "rmce0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NExternalMix1BB.MDC2025ar_best_v1_1.root");
    datasets_["rmc_ext_1n"] = DatasetInfo_t(7000000000, 2965072, 0.,   1., rate_rmc_ps_1n_80            , "rmce1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NExternalMix1BB.MDC2025ar_best_v1_1.root");
    datasets_["rmc_int_0n"] = DatasetInfo_t(  50000000, 1211753, 0.,   1., rate_rmc_ps_0n_80*rmc_int_br_, "rmci0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NInternalMix1BB.MDC2025ar_best_v1_1.root");
    datasets_["rmc_int_1n"] = DatasetInfo_t(  50000000,  507641, 0.,   1., rate_rmc_ps_1n_80*rmc_int_br_, "rmci1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NInternalMix1BB.MDC2025ar_best_v1_1.root");
    datasets_["rpc_ext"]    = DatasetInfo_t(5000000000,  458818, 0.,   1., rate_phys_rpc                , "rpce1b1s5r0100", "nts.mu2e.RPCExternalPhysicalMix1BB.MDC2025au_best_v1_1.root");
    datasets_["rpc_int"]    = DatasetInfo_t( 125000000, 1899806, 0.,   1., rate_phys_rpc_int            , "rpci1b1s5r0100", "nts.mu2e.RPCExternalPhysicalMix1BB.MDC2025au_best_v1_1.root");

    datasets_["data_mds3c"] = DatasetInfo_t(1.    , 4921433, 0.,   1., 1., "mds3cb1s5r0100", "nts.mu2e.ensembleMDS3cMix1BB.MDC2025ar_best_v1_1.root"              );
  }
}

DatasetInfo_t get_dataset_info(TString name) {
  // Initialize the datasets if not already
  if(datasets_.size() == 0) init_dataset_info();
  if(datasets_.count(name) != 0) return datasets_[name];
  cout << __func__ << ": No dataset with name " << name.Data() << " found!\n";
  return DatasetInfo_t();
}

//---------------------------------------------------------------------------------------------------------------------------
void set_style(const TString name, TString& title, int& color) {
  title = name;
  color = kRed;
  if(name == "signal") {
    color = kBlue;
    title = "Signal";
  } else if(name == "mumem") {
    color = kBlue;
    title = "#mu^{-}#rightarrowe^{-}";
  } else if(name == "mumep") {
    color = kBlue;
    title = "#mu^{-}#rightarrowe^{+}";
  } else if(name == "rmc") { // signal version
    title = "RMC";
    color = kBlue;
  } else if(name == "dio") {
    title = "DIO";
    color = kRed-7;
  } else if(name == "ipa_dio") {
    title = "IPA DIO";
    color = kRed-3;
  } else if(name.BeginsWith("cosmic")) {
    title = "Cosmic ray";
    color = kOrange;
  } else if(name == "pbar") {
    title = "Antiproton";
    color = kGreen-6;
  } else if(name == "rpc_ext") {
    title = (combine_rpc_) ? "RPC" : "RPC (external)";
    color = kMagenta-10;
  } else if(name == "rpc_int") {
    title = (combine_rpc_) ? "RPC" : "RPC (internal)";
    color = kMagenta+1;
  } else if(name.BeginsWith("rmc_ext")) {
    title = (combine_rmc_) ? "RMC" : "RMC (external)";
    color = kAtlantic+2;
  } else if(name.BeginsWith("rmc_int")) {
    title = (combine_rmc_) ? "RMC" : "RMC (internal)";
    color = kAtlantic;
  } else if(name.BeginsWith("mnbs")) {
    title = "Pileup";
    color = kRed;
  }
}

#endif
