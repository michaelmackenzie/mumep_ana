#ifndef __mumep_ana_ana_CRVHist_t_hh
#define __mumep_ana_ana_CRVHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace mumep_ana {

  struct CRVClusterHist_t {
    TH1F* fSector;
    TH1F* fFirstBar; // bar # of the first pulse
    TH1F* fNPulses;
    TH1F* fNPe;   // N(PE) - apparently, the sum
    TH1F* fNPePP; // N(PE) per pulse
    TH1F* fStartTime;
    TH1F* fEndTime;
    TH1F* fWidth;
    TH2F* fXVsZ;
    TH2F* fYVsZ;
    TH1F* fCorrTime;
    TH1F* fCorrTimeProp;
    TH1F* fCorrTimeToF;
    TH1F* fApproxTimeST;
    TH1F* fApproxTimeCalo;
    TH1F* fApproxTimeExtrap;
    TH1F* fApproxTimeSTToFront;
    TH1F* fApproxTimeCaloToFront;
    TH1F* fApproxTimeExtrapToFront;
    TH1F* fBarsOneEnd;
    TH1F* fCrvPropdT;
    TH1F* fNSectors;
    TH1F* fBarsTwoEnd;
    TH1F* fNDiffLSectors;
    TH1F* fStubSlope;
    TH1F* fStubSlopeChi2;
    TH1F* fStubSlopeDelta;
    TH1F* fStubQN;
    TH1F* fStubSlopeMCProduct;
    TH1F* fMCTime;
    TH1F* fMCdT[2];
  };

  struct CRVHist_t {
    CRVClusterHist_t fClusterHist;
  };
} // namespace mumep_ana
#endif
