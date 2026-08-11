#ifndef __mumep_ana_ana_TrackHist_t_hh
#define __mumep_ana_ana_TrackHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace mumep_ana {

  struct TrackHist_t {
    TH1F* fP[2];
    TH1F* fObs; // observable being fit
    TH1F* fPt;
    TH1F* fPCorr; // attempting to correct the momentum
    TH1F* fPCenter[2];
    TH1F* fPExit;
    TH1F* fPST[2];
    TH1F* fPSTDiff;
    TH1F* fPSTApproxDiff;
    TH1F* fPExitDiff;
    TH1F* fPTrkFront; // tracker front surface explicitly
    TH1F* fT0;
    TH1F* fT0Err;
    TH1F* fD0;
    TH1F* fDP;
    TH2F* fDPvsP;
    TH2F* fDPvsNH;
    TH1F* fDPCorr;
    TH1F* fDPCenter;
    TH1F* fChi2NDof;
    TH1F* fFitCons[2];
    TH1F* fFitMomErr; // estimated uncertainty
    TH1F* fTanDip;
    TH1F* fCosTheta;
    TH1F* fRadius;
    TH1F* fRMax;
    TH1F* fNActive;
    TH1F* fNActiveFrac;
    TH1F* fTrkQual[2]; // 0: Defined in TStnTrack; 1: Local evaluation
    TH1F* fPID[2];     // 0: Normal PID; 1: Tracker-only PID
    TH1F* fCosmicID;
    TH1F* fClusterE;
    TH1D* fClusterDisk;
    TH1F* fDt;
    TH1F* fEp;
    TH1F* fTZSlope;
    TH1F* fTZSlopeSig;
    TH1F* fTZSlopeRatio;
    TH1F* fBestAlg;
    TH1F* fAlgMask;
    TH1F* fSTBoundary;
    TH1F* fSTInters;
    TH1F* fIPAInters;
    TH1F* fOPAInters;
    TH1F* fTrackID;
    TH1F* fExlTrackID; // flagging what is exclusively cut by a given ID bit

    // Matched CRV cluster info
    TH1F* fCRVDeltaT;        // time difference after corrections
    TH1F* fCRVDeltaTCRV;     // time difference from stub time
    TH1F* fCRVDeltaTST;      // assuming path from ST
    TH1F* fCRVDeltaTCalo[2]; // assuming path from Calo: 0: electron; 1: muon
    TH1F* fCRVDeltaTExtrap;  // assuming path from extrapolation
    TH1F* fCRVMinDeltaT;     // minimum delta T between ST and Calo paths
    TH1F* fCRVExtrapZ;
    TH2F* fCRVXZ;
    TH2F* fCRVYZ;
    TH2F* fCRVdTZ;       // using the corrected time
    TH2F* fCRVdTZCRV;    // using the stub time
    TH2F* fCRVdTZExtrap; // using the extrapolated time

    // Matched upstream track info
    TH1F* fUpstreamDt;
    TH1F* fUpstreamDp;

    // MC truth information
    TH1F* fMCPFront; // P(front of tracker)
    TH1F* fMCPStOut; // P(ST exit)
    TH1F* fMCPCenterDiff;
    TH1F* fMCPStDiff;
    TH1F* fMCPStDiffDiff;       // Reco - MC (P(ST exit) - P(Tracker Front))
    TH1F* fMCApproxPStDiffDiff; // Reco - MC (P(ST exit) - P(Tracker Front))
    TH1F* fMCGenE;              // E(sim particle) at generation
    TH1F* fMCPGenEDiff;         // P(front of tracker) - Gen(energy)
    TH1F* fMCPSig;              // error / uncertainty
    TH1F* fMCPdg[2];
    TH1F* fMCStrawHits;
    TH1F* fMCGoodHits;
    TH1F* fMCTrajectory;
    TH1F* fMCSimProc;
  };

} // namespace mumep_ana
#endif
