#ifndef __mumep_ana_ana_ClusterHist_t_hh
#define __mumep_ana_ana_ClusterHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace mumep_ana {

  struct ClusterHist_t {
    TH1D* fDiskID;
    TH1F* fEnergy;
    TH1F* fT0;
    TH1F* fRow;
    TH1F* fCol;
    TH1F* fX;
    TH1F* fY;
    TH1F* fZ;
    TH1F* fR;
    TH1F* fNCr0; // all clustered
    TH1F* fNCr1; // above 1MeV
    TH1F* fYMean;
    TH1F* fZMean;
    TH1F* fSigY;
    TH1F* fSigZ;
    TH1F* fSigR;
    TH1F* fFrE1;
    TH1F* fFrE2;
    TH1F* fSigE1;
    TH1F* fSigE2;
    TH1F* fTimeRMS;
    TH1F* fMaxR;
    TH1F* fE9OverE;
    TH1F* fE25OverE;
    TH1F* fRingEOverE;
    TH1F* fRingEOverE1;
    TH1F* fOutRingE;
    TH1F* fOutRingEOverE;
    TH1F* fNCoreCrystals;
    TH1F* fCoreEnergy;
    TH1F* fCoreEnergyFrac;

    TH1F* fMCSimEDep;
    TH1F* fMCSimMomIn;
    TH1F* fMCSimPdg;
    TH1F* fMCSimPdgName;
    TH1F* fMCSimEStart;
    TH1F* fMCEDep;
    TH1F* fMCTime;
    TH1F* fMC_dE;
    TH1F* fMC_dt;
    TH1F* fMC_dGenE;
    TH1F* fMCGenE;
    TH2F* fMCEvsE;
    TH2F* fMCGenEvsE;
  };

} // namespace mumep_ana
#endif
