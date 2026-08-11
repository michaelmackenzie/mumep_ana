#ifndef __mumep_ana_ana_EventHist_t_hh
#define __mumep_ana_ana_EventHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace mumep_ana {

  struct EventHist_t {
    TH1F* fInstLumi[2];
    TH1F* fInstLumiApr;    // only filled if apr triggered event
    TH1F* fInstLumiCpr;    // only filled if cpr triggered event
    TH1F* fInstLumiAprCpr; // filled if either apr or cpr triggered event
    TH1F* fEventWeight[2];
    TH1F* fNAprHelices;
    TH1F* fNCprHelices;
    TH1F* fNHelices;
    TH1F* fNMatchedHelices;
    TH1F* fNAprTracks;
    TH1F* fNCprTracks;
    TH1F* fNTracks;
    TH1F* fNUeTracks;
    TH1F* fNDmuTracks;
    TH1F* fNUmuTracks;
    TH1F* fNCRVClusters;
    TH1F* fNGoodCRVClusters;
    TH1F* fNonCRVVetoID;
    TH1F* fTrackerHits;
    TH1F* fCaloHits;

    TH1F* fTrigBits[2];
    TH1F* fTrigPaths[2];
    TH2F* fTrigOverlap[2];
    TH1F* fNTriggerable; // N(triggerable sim particles)
    TH1D* fTriggered;

    // Primary process info
    TH1F* fPrimaryCode;
    TH1F* fPrimaryType;
    TH1F* fPrimaryGenE;

    // Process specific info
    TH1F* fRMCEnergy; // RMC photon energy

    // TriggerAna info
    TH2F* fAprPOTVsP;
    TH2F* fCprPOTVsP;
    TH2F* fAprHlxPOTVsP;
    TH2F* fCprHlxPOTVsP;
  };

} // namespace mumep_ana
#endif
