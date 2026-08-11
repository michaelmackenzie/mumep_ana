// Information to store in an output tree
#ifndef __mumep_ana_ana_Tree_t__
#define __mumep_ana_ana_Tree_t__

#include "TTree.h"

#include <vector>

namespace mumep_ana {

  struct Tree_t {

    // Event info
    int fRun;
    int fSubRun;
    int fEvent;
    float fWeight = 1.f;
    float fTrain;

    // Track information
    float fTrkP;
    float fTrkT0;
    float fTrkD0;
    float fTrkTanDip;
    float fTrkCosTheta;
    float fTrkCluster;
    float fTrkEP;
    float fTrkDt;
    float fTrkRMax;
    float fTrkUsDt;
    float fTrkActiveRatio;
    float fTrkNullRatio;
    float fTrkFitCon;
    float fTrkLogFitCon;
    float fTrkTZSlope;
    float fTrkTZSlopeSig;
    float fTrkTZSlopeRatio;
    float fTrkPExitDiff;
    float fTrkQual;
    float fTrkPID;
    float fTrkOnlyPID;
    float fTrkCosmicID;
    float fTrkCharge;
    float fTrkMCDp;
    float fTrkMCPDG;

    // CRV information
    float fCRVZ;
    float fCRVDeltaT;
    float fCRVNPulses;
    float fCRVNPe;
    float fCRVNPePP;

    //----------------------------------------------------------
    // Track quality info

    float fTrkQual_nactive;
    float fTrkQual_activehitsfraction;
    float fTrkQual_nullhitsfraction;
    float fTrkQual_activematsitesfraction;
    float fTrkQual_fitcons;
    float fTrkQual_momerr;
    float fTrkQual_t0err;

    //----------------------------------------------------------
    // Specific for TReflectionAna

    // Upstream track info
    float fRefl_p_us;
    float fRefl_pt_us;
    float fRefl_cos_us;
    float fRefl_d0_us;
    float fRefl_r_us;
    float fRefl_rmax_us;
    float fRefl_rmin_us;
    float fRefl_ep_us;
    int fRefl_nhits_us;
    int fRefl_q_us;
    int fRefl_stboundary_us;
    int fRefl_stinters_us;
    int fRefl_triggered_us;

    // Downstream track info
    float fRefl_p_ds;
    float fRefl_pt_ds;
    float fRefl_cos_ds;
    float fRefl_d0_ds;
    float fRefl_r_ds;
    float fRefl_rmax_ds;
    float fRefl_rmin_ds;
    float fRefl_ep_ds;
    int fRefl_nhits_ds;
    int fRefl_q_ds;
    int fRefl_stboundary_ds;
    int fRefl_stinters_ds;
    int fRefl_triggered_ds;

    // Comparison info
    float fRefl_dp;
    float fRefl_dt0;
    float fRefl_dtfront;

    // Photon CNN info
    std::vector<float> fCrystalE;
    std::vector<float> fCrystalT;

    // TTree
    TTree* fTree = nullptr;

    Tree_t() { reset(); }

    void reset() {

      fRun = 0;
      fSubRun = 0;
      fEvent = 0;
      fWeight = 1.f;
      fTrain = 0.f;

      fTrkP = 0.f;
      fTrkT0 = 0.f;
      fTrkD0 = 0.f;
      fTrkTanDip = 0.f;
      fTrkCluster = 0.f;
      fTrkEP = 0.f;
      fTrkDt = 0.f;
      fTrkRMax = 0.f;
      fTrkActiveRatio = 0.f;
      fTrkNullRatio = 0.f;
      fTrkFitCon = 0.f;
      fTrkLogFitCon = 0.f;
      fTrkTZSlope = 0.f;
      fTrkTZSlopeSig = 0.f;
      fTrkTZSlopeRatio = 0.f;
      fTrkMCDp = 0.f;
      fTrkMCPDG = 0.f;

      fCRVZ = 0.f;
      fCRVDeltaT = 0.f;
      fCRVNPulses = 0.f;
      fCRVNPe = 0.f;
      fCRVNPePP = 0.f;

      fTrkQual_nactive = 0.f;
      fTrkQual_activehitsfraction = 0.f;
      fTrkQual_nullhitsfraction = 0.f;
      fTrkQual_activematsitesfraction = 0.f;
      fTrkQual_fitcons = 0.f;
      fTrkQual_momerr = 0.f;
      fTrkQual_t0err = 0.f;

      fRefl_p_us = 0.f;
      fRefl_pt_us = 0.f;
      fRefl_cos_us = 0.f;
      fRefl_d0_us = 0.f;
      fRefl_r_us = 0.f;
      fRefl_rmax_us = 0.f;
      fRefl_rmin_us = 0.f;
      fRefl_ep_us = 0.f;
      fRefl_nhits_us = 0;
      fRefl_q_us = 0;
      fRefl_stboundary_us = 0;
      fRefl_stinters_us = 0;
      fRefl_triggered_us = 0;
      fRefl_p_ds = 0.f;
      fRefl_pt_ds = 0.f;
      fRefl_cos_ds = 0.f;
      fRefl_d0_ds = 0.f;
      fRefl_r_ds = 0.f;
      fRefl_rmax_ds = 0.f;
      fRefl_rmin_ds = 0.f;
      fRefl_ep_ds = 0.f;
      fRefl_nhits_ds = 0;
      fRefl_q_ds = 0;
      fRefl_stboundary_ds = 0;
      fRefl_stinters_ds = 0;
      fRefl_triggered_ds = 0;
      fRefl_dp = 0.f;
      fRefl_dt0 = 0.f;
      fRefl_dtfront = 0.f;

      fCrystalE.clear();
      fCrystalT.clear();
    }
  };
} // namespace mumep_ana
#endif
