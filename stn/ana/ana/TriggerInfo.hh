///////////////////////////////////////////////////////////////////////////////
// Trigger path name information from a trigger bit
///////////////////////////////////////////////////////////////////////////////
#ifndef __mumep_ana_ana_TriggerInfo_hh__
#define __mumep_ana_ana_TriggerInfo_hh__

// c++ includes
#include <iostream>
#include <map>

// ROOT includes
#include "TString.h"

namespace mumep_ana {
  class TriggerInfo {
  public:
    static TString BitToName(const int bit) {
      switch(bit) {
      case 100:
        return "tpr_TrkDe_80m70p_D0200";
      case 110:
        return "tpr_TrkDe_80m70p";
      case 125:
        return "tpr_TrkDe_50_D0200";
      case 130:
        return "tpr_HlxDe_70m50p";
      case 131:
        return "tpr_HlxUe_50m30p";
      case 140:
        return "tpr_HlxDe_30p_IPA";
      case 141:
        return "tpr_HlxDe_30p_IPAPhi";
      case 150:
        return "cpr_TrkDe_80m70p_D0200";
      case 151:
        return "cpr_TrkDe_75m70p_D0200";
      case 152:
        return "cpr_TrkDe_75_D0200";
      case 153:
        return "cpr_TrkDe_70_D0200";
      case 160:
        return "cpr_TrkDe_80m70p";
      case 161:
        return "cpr_TrkDe_75m70p";
      case 162:
        return "cpr_TrkDe_75";
      case 163:
        return "cpr_TrkDe_70";
      case 175:
        return "cpr_TrkDe_50_D0200";
      case 180:
        return "cpr_HlxDe_50";
      case 181:
        return "cpr_HlxUe_40";
      case 200:
        return "apr_TrkDe_80m70p_D0200";
      case 201:
        return "apr_TrkDe_75m70p_D0200";
      case 202:
        return "apr_TrkDe_75_D0200";
      case 203:
        return "apr_TrkDe_70_D0200";
      case 210:
        return "apr_TrkDe_80m70p";
      case 211:
        return "apr_TrkDe_75m70p";
      case 212:
        return "apr_TrkDe_75";
      case 213:
        return "apr_TrkDe_70";
      case 230:
        return "apr_TrkUe_80m70p";
      case 240:
        return "apr_TwoTrkDe_80m70p_D0200";
      case 250:
        return "apr_TwoTrkDe_50";
      case 255:
        return "apr_Hlx_50_Hlx_30";
      case 260:
        return "apr_TrkDe_50_D0200";
      case 265:
        return "apr_Hlx_50";
      case 266:
        return "apr_Hlx_30";
      case 275:
        return "mpr_TrkDe_80m70p_D0200";
      case 400:
        return "calo_photon";
      case 401:
        return "calo_NNCE";
      case 402:
        return "calo_cosmic";
      case 420:
        return "calo_RMC";
      case 425:
        return "calo_cluster_50";
      case 426:
        return "calo_cluster_60";
      case 427:
        return "calo_cluster_70";
      case 428:
        return "calo_cluster_75";
      case 429:
        return "calo_cluster_80";
      case 500:
        return "cst_TimeCluster";
      case 510:
        return "apr_TC";
      case 511:
        return "apr_TC_calo";
      case 520:
        return "cst_CosmicTrackSeed";
      case 600:
        return "minBias_SDCount";
      case 610:
        return "minBias_CDCount";
      case 700:
        return "calo_N0Source";
      case 800:
        return "lumiStream";

        // Old bits
        // case 210: return "mprDe_highP_stopTarg";
        // case 100: return "tprDe_highP_stopTarg";
        // case 101: return "tprDe_highP";
        // case 110: return "tprDe_lowP_stopTarg";
        // case 120: return "tprHelixDe_ipa";
        // case 121: return "tprHelixDe_ipa_phiScaled";
        // case 130: return "tprHelixDe";
        // case 131: return "tprHelixUe";
        // case 150: return "cprDe_highP_stopTarg";
        // case 151: return "cprDe_highP";
        // case 160: return "cprDe_lowP_stopTarg";
        // case 170: return "cprHelixDe";
        // case 171: return "cprHelixUe";
        // case 180: return "apr_highP_stopTarg";
        // case 181: return "apr_highP";
        // case 185: return "apr_ue_highP";
        // case 190: return "apr_lowP_stopTarg";
        // case 191: return "apr_lowP_stopTarg_multiTrk";
        // case 192: return "apr_lowP_multiHelix";
        // case 195: return "aprHelix";
        // case 200: return "caloFast_photon";
        // case 201: return "caloFast_NNCE";
        // case 202: return "caloFast_cosmic";
        // case 220: return "caloFast_RMC";
        // case 310: return "cstTimeCluster";
        // case 320: return "cstCosmicTrackSeed";
        // case 400: return "minBias_SDCount";
        // case 410: return "minBias_CDCount";
        // case 500: return "caloHitRec_N0Source";
      }
      return Form("Unknown-%i", bit);
    }

    static double BitToPrescale(const int bit) { // for testing prescale levels
      switch(bit) {
      case 100:
        return -1; // "tpr_TrkDe_80m70p_D0200"
      case 110:
        return -1; // "tpr_TrkDe_80m70p"
      case 125:
        return -1; // "tpr_TrkDe_50_D0200"
      case 130:
        return -1; // "tpr_HlxDe_70m50p"
      case 131:
        return -1; // "tpr_HlxUe_50m30p"
      case 140:
        return -1; // "tpr_HlxDe_30p_IPA"
      case 141:
        return -1; // "tpr_HlxDe_30p_IPAPhi"
      case 150:
        return 1; // "cpr_TrkDe_80m70p_D0200"
      case 151:
        return 1; // "cpr_TrkDe_75m70p_D0200"
      case 152:
        return 1; // "cpr_TrkDe_75_D0200"
      case 153:
        return 10; // "cpr_TrkDe_70_D0200"
      case 160:
        return 1; // "cpr_TrkDe_80m70p"
      case 161:
        return 1; // "cpr_TrkDe_75m70p"
      case 162:
        return 1; // "cpr_TrkDe_75"
      case 163:
        return 10; // "cpr_TrkDe_70"
      case 175:
        return 1000; // "cpr_TrkDe_50_D0200"
      case 180:
        return 5000; // "cpr_HlxDe_50"
      case 181:
        return 2000; // "cpr_HlxUe_40"
      case 200:
        return 1; // "apr_TrkDe_80m70p_D0200"
      case 201:
        return 1; // "apr_TrkDe_75m70p_D0200"
      case 202:
        return 1; // "apr_TrkDe_75_D0200"
      case 203:
        return 10; // "apr_TrkDe_70_D0200"
      case 210:
        return 1; // "apr_TrkDe_80m70p"
      case 211:
        return 1; // "apr_TrkDe_75m70p"
      case 212:
        return 1; // "apr_TrkDe_75"
      case 213:
        return 10; // "apr_TrkDe_70"
      case 230:
        return 1; // "apr_TrkUe_80m70p"
      case 240:
        return 1; // "apr_TwoTrkDe_80m70p_D0200"
      case 250:
        return 1; // "apr_TwoTrkDe_50"
      case 260:
        return 10000; // "apr_TrkDe_50_D0200"
      case 265:
        return 20000; // "apr_Hlx_50"
      case 266:
        return 30000; // "apr_Hlx_30"
      case 275:
        return -1; // "mpr_TrkDe_80m70p_D0200"
      case 400:
        return 20000; // "calo_photon"
      case 401:
        return 1; // "calo_NNCE"
      case 402:
        return 10; // "calo_cosmic"
      case 420:
        return 10; // "calo_RMC"
      case 425:
        return 100000; // "calo_cluster_50"
      case 426:
        return 50000; // "calo_cluster_60"
      case 427:
        return 5000; // "calo_cluster_70"
      case 428:
        return 500; // "calo_cluster_75"
      case 429:
        return 100; // "calo_cluster_80"
      case 500:
        return -1; // "cst_TimeCluster"
      case 520:
        return -1; // "cst_CosmicTrackSeed"
      case 600:
        return -1; // "minBias_SDCount"
      case 610:
        return -1; // "minBias_CDCount"
      case 700:
        return -1; // "calo_N0Source"
      case 800:
        return -1; // "lumiStream"
      }
      return -1;
    }
  };
} // namespace mumep_ana
#endif
