#ifndef __mumep_ana_ana_CRVStubPar_t__
#define __mumep_ana_ana_CRVStubPar_t__

#include "Stntuple/ana/ParBase_t.hh"
#include "Stntuple/obj/TCrvCoincidenceCluster.hh"

namespace mumep_ana {
  class CRVStubPar_t : public stntuple::ParBase_t {
  public:
    TCrvCoincidenceCluster* fCluster = nullptr;
    float fTime;
    float fZ;

    float fCorrTime;
    float fCorrTimeTof;
    float fCorrTimeProp;
    float fExtrapZ; // extrapolating along the stub slope to the y-axis

    // approximate times just using the cluster location + time to standard points
    float fApproxTimeST;
    float fApproxTimeCalo;
    float fApproxTimeExtrap;
    float fApproxTimeSTToFront;     // time at tracker front assuming path from ST
    float fApproxTimeCaloToFront;   // time at tracker front assuming path from Calo
    float fApproxTimeExtrapToFront; // time at tracker front using stub slope + extrapolation

    int fSector;   // sector number (not type)
    int fFirstBar; // bar # of the first pulse
    int fTwoEndBars;
    int fTotalBars;
    int fNSectors;
    int fNDiffLSectors;
    int fStubQN;

    float fXCorrected;
    float fNPePP;

    float fTCorrAana; // time with all analysis-based corrections (evolving)

    float fStubDYDZ;
    float fStubSlopeChi2;
    float fStubDYDZMC;
    float fStubSlopeMCProduct;

    void reset() {
      fCluster = nullptr;
      fTime = 0;
      fZ = 0;
      fCorrTime = 0;
      fCorrTimeTof = 0;
      fCorrTimeProp = 0;
      fExtrapZ = 0;
      fApproxTimeST = 0;
      fApproxTimeCalo = 0;
      fApproxTimeExtrap = 0;
      fApproxTimeSTToFront = 0;
      fApproxTimeCaloToFront = 0;
      fApproxTimeExtrapToFront = 0;
      fSector = 0;
      fFirstBar = 0;
      fTwoEndBars = 0;
      fTotalBars = 0;
      fNSectors = 0;
      fNDiffLSectors = 0;
      fStubQN = 0;
      fXCorrected = 0;
      fNPePP = 0;
      fTCorrAana = 0;
      fStubDYDZ = 0;
      fStubSlopeChi2 = 0;
      fStubDYDZMC = 0;
      fStubSlopeMCProduct = 0;
    }
  };
} // namespace mumep_ana
#endif
