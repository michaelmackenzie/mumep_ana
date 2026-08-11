#ifndef __mumep_ana_ana_EventPar_t__
#define __mumep_ana_ana_EventPar_t__

namespace mumep_ana {

  struct EventPar_t {

    float fWeight = 1.f;

    float fInstLum = 1.f;
    int fNAprHelices = 0;
    int fNCprHelices = 0;
    int fNOfflineHelices = 0;
    int fNAprTracks = 0; // for trigger analysis
    int fNCprTracks = 0; // for trigger analysis
    int fNAprUeTracks = 0;
    int fNTracks = 0;
    int fNGoodTracks = 0;
    int fNUeTracks = 0;
    int fNUmuTracks = 0;
    int fNDmuTracks = 0;
    int fNCRVClusters = 0;
    int fNGoodCRVClusters = 0;
    int fTrackerHits = 0;
    int fCaloHits = 0;

    bool fPassedAprPath = false;
    bool fPassedCprPath = false;
    bool fTriggered = false;
    int fNTriggerable = 0;

    int fNonCRVVetoID = 0;

    // Process-specific information
    float fRMCEnergy = 0.f;
  };
} // namespace mumep_ana
#endif
