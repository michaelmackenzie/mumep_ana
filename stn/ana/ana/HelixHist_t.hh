#ifndef __mumep_ana_ana_HelixHist_t_hh
#define __mumep_ana_ana_HelixHist_t_hh

#include "TH1.h"

namespace mumep_ana {

  struct HelixHist_t {
    TH1F* fNHits;
    TH1F* fHelicity;
    TH1F* fClusterTime;
    TH1F* fClusterEnergy;
    TH1F* fRadius; // fabs(1/omega)
    TH1F* fRMax;
    TH1F* fMom;
    TH1F* fPt;
    TH1F* fGenMom;
    TH1F* fGenID;
    TH2F* fGenRZ;
    TH1F* fDp;
    TH1F* fDpT;
    TH1F* fLambda; // dz/dphi
    TH1F* fTanDip;
    TH1F* fT0;
    TH1F* fT0Err;
    TH1F* fD0;
    TH1F* fBestAlg;
    TH1F* fAlgMask;
    TH1F* fChi2XY;
    TH1F* fChi2ZPhi;
    TH1F* fTZSlope;
    TH1F* fTZSlopeErr;
    TH1F* fTZSlopeSig;
    TH1F* fTZSlopeSigCDF[2]; // two CDF directions
  };

} // namespace mumep_ana
#endif
