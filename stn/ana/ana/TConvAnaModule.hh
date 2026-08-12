///////////////////////////////////////////////////////////////////////////////
// Module to search for muon converion signals
///////////////////////////////////////////////////////////////////////////////
#ifndef mumep_ana_ana_TConvAnaModule_hh
#define mumep_ana_ana_TConvAnaModule_hh

// local includes
#include "mumep_ana/stn/ana/CRVHist_t.hh"
#include "mumep_ana/stn/ana/EventHist_t.hh"
#include "mumep_ana/stn/ana/EventPar_t.hh"
#include "mumep_ana/stn/ana/MVATools.hh"
#include "mumep_ana/stn/ana/SysHist_t.hh"
#include "mumep_ana/stn/ana/Systematics.hh"
#include "mumep_ana/stn/ana/TAnaModule.hh"
#include "mumep_ana/stn/ana/TrackHist_t.hh"
#include "mumep_ana/stn/ana/TrackPar_t.hh"
#include "mumep_ana/stn/ana/Tree_t.hh"

// Stntuple includes
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvCoincidence.hh"
#include "Stntuple/obj/TCrvCoincidenceCluster.hh"
#include "Stntuple/obj/TGenParticle.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimParticle.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStnCluster.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnHelix.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrack.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnTriggerBlock.hh"

// mu2e includes
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "ArtAnalysis/TrkDiag/inc/TrackPID_v1.hxx"

// ROOT includes
#include "TRandom3.h"
#include "TTree.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

namespace TMVA_SOFIE_TrackPID_v1 {
  class Session;
}

namespace mumep_ana {
  class TConvAnaModule : public TAnaModule {
  public:
    enum { kNHistSets = 10000, kMaxCRVStubs = 200, kMaxTracks = 50 };
    enum { kCRVVetoOffset = 1000, kTimeCutOffset = 2000 }; // control region offsets
    enum { kBackground, kSignal, kData };                  // data types
    enum { kDIO, kCosmic, kExternalRPC, kInternalRPC, kExternalRMC, kInternalRMC, kSimplePBar }; // background types

    struct Hist_t {
      mumep_ana::EventHist_t* fEvent[kNHistSets];
      mumep_ana::GenpHist_t* fGenp[kNHistSets];
      mumep_ana::SimpHist_t* fSimp[kNHistSets];
      mumep_ana::TrackHist_t* fTrack[kNHistSets];
      mumep_ana::HelixHist_t* fHelix[kNHistSets];
      mumep_ana::ClusterHist_t* fCluster[kNHistSets];
      mumep_ana::CRVHist_t* fCRV[kNHistSets];
      mumep_ana::SysHist_t* fSys[kNHistSets];
      TTree* fTree[kNHistSets];
      Hist_t() {
        for(int i = 0; i < kNHistSets; ++i) {
          fEvent[i] = nullptr;
          fGenp[i] = nullptr;
          fSimp[i] = nullptr;
          fTrack[i] = nullptr;
          fHelix[i] = nullptr;
          fCluster[i] = nullptr;
          fCRV[i] = nullptr;
          fSys[i] = nullptr;
          fTree[i] = nullptr;
        }
      }
    };

    // structure for mapping corresponding helices and tracks
    struct HelixPair_t {
      TStnHelix* fAprHelix = nullptr;
      TStnHelix* fCprHelix = nullptr;
      TStnHelix* fOfflineHelix = nullptr;
      TStnTrack* fAprTrack = nullptr;
      TStnTrack* fCprTrack = nullptr;
      TStnTrack* fOfflineTrack = nullptr;
    };

    //-----------------------------------------------------------------------------
    //  data members
    //-----------------------------------------------------------------------------
  public:
    TGenpBlock* fGenpBlock = nullptr;
    TString fGenpBlockName;
    TSimpBlock* fSimpBlock = nullptr;
    TString fSimpBlockName;
    TCrvClusterBlock* fCRVBlock = nullptr;
    TString fCRVBlockName;
    TStnClusterBlock* fClusterBlock = nullptr;
    TString fClusterBlockName;

    TStnHelixBlock* fAprHelixBlock = nullptr;
    TStnHelixBlock* fCprDeHelixBlock = nullptr;
    TStnHelixBlock* fOfflineDeHelixBlock = nullptr;
    TStnHelixBlock* fOfflineUeHelixBlock = nullptr;
    TStnTrackBlock* fAprTrackBlock = nullptr;
    TStnTrackBlock* fCprTrackBlock = nullptr;
    TStnTrackBlock* fOfflineDeTrackBlock = nullptr;
    TStnTrackBlock* fOfflineUeTrackBlock = nullptr;
    TStnTrackBlock* fOfflineDmuTrackBlock = nullptr;
    TStnTrackBlock* fOfflineUmuTrackBlock = nullptr;
    TString fAprHelixBlockName;
    TString fCprDeHelixBlockName;
    TString fOfflineDeHelixBlockName;
    TString fOfflineUeHelixBlockName;
    TString fAprTrackBlockName;
    TString fCprTrackBlockName;
    TString fOfflineDeTrackBlockName;
    TString fOfflineUeTrackBlockName;
    TString fOfflineDmuTrackBlockName;
    TString fOfflineUmuTrackBlockName;

    TStnTrack* fTrack = nullptr;
    TStnTrack* fOfflineTrack = nullptr;
    TStnHelix* fHelix = nullptr;
    TStnCluster* fCluster = nullptr;
    TGenParticle* fGen = nullptr;
    TSimParticle* fSim = nullptr;
    CosmicVetoData_t fCosmicVetoData;

    mumep_ana::EventPar_t fEvtPar;
    mumep_ana::TrackPar_t fTrkPar;
    mumep_ana::HelixPar_t fHlxPar;
    mumep_ana::ClusterPar_t fClusterPar;
    mumep_ana::CRVStubPar_t fCRVStubPar[kMaxCRVStubs];
    mumep_ana::CRVStubPar_t* fMCCRVStubPar = nullptr; // Stub corresponding to the actual cosmic
    Hist_t fHist;
    Tree_t fTreeData;
    Systematics fSystematics;
    TRandom3 fRand;
    int fBookTrees = 1;
    int fEvaluateMVAs = 0;
    TMVA::Reader* fTrkQual = nullptr;
    int fTrkQualVersion = 0;
    TMVA::Reader* fPID = nullptr;
    int fPIDVersion = 0;
    TMVA::Reader* fTrkPID = nullptr;
    int fTrkPIDVersion = 0;
    TMVA::Reader* fCosmicID = nullptr;
    int fCosmicIDVersion = 0;
    TMVA_SOFIE_TrackPID_v1::Session* fOfflinePID = nullptr;
    int fDoPIDTrees = 0;
    int fDoCosmicIDTrees = 0;

    int fIgnoreCRV = 0;      // ignore CRV info
    int fBatchModeSim = 1;   // intensity mode the dataset is generated with
    int fBatchModeUse = 1;   // intensity mode to model
    int fUseBeamWeights = 1; // apply beam-intensity reweighting
    int fUsePionWeights = 1; // apply pion lifetime weights to RPC events
    int fRMCVersion = kPlestid; // RMC spectrum to use
    int fIncludeSys = 1;     // include systematic histograms
    int fFillVerboseSys = 0; // include debug systematic histograms

    TString fDataset = ""; // dataset name
    int fDataType = kBackground;
    int fBackgroundType = -1; // specific background sample
    int fApplyID = 0;         // apply track ID when filling histograms
    float fPMin = -1.f;       // track selection window
    float fPMax = -1.f;

    //-----------------------------------------------------------------------------
    //  functions
    //-----------------------------------------------------------------------------
  public:
    TConvAnaModule(const char* name = "mumep_ana_mumep_ana", const char* title = "mumep_ana");
    ~TConvAnaModule();

    //-----------------------------------------------------------------------------
    // accessors
    //-----------------------------------------------------------------------------
    Hist_t* GetHist() { return &fHist; }

    //-----------------------------------------------------------------------------
    // overloaded methods of TStnModule
    //-----------------------------------------------------------------------------
    int BeginJob();
    int BeginRun();
    int Event(int ientry);
    int EndJob();
    //-----------------------------------------------------------------------------
    // other methods
    //-----------------------------------------------------------------------------
    void BookEventHistograms(mumep_ana::EventHist_t* Hist, const char* Folder);
    void BookHelixHistograms(mumep_ana::HelixHist_t* Hist, const char* Folder);
    void BookSystematicHistograms(mumep_ana::SysHist_t* Hist, const char* Folder);
    void BookTree(TTree* tree);
    void BookHistograms();
    void FillEventHistograms(mumep_ana::EventHist_t* Hist, mumep_ana::EventPar_t* EvtPar, float Weight = 1.f);
    void FillHelixHistograms(mumep_ana::HelixHist_t* Hist, mumep_ana::HelixPar_t* HlxPar, float Weight = 1.f);
    void FillSystematicHistograms(mumep_ana::SysHist_t* Hist, mumep_ana::TrackPar_t* TrkPar, float Weight = 1.f);
    void FillTree(TTree* tree);
    void FillAllHistograms(const int index);
    void FillHistograms();

    //-----------------------------------------------------------------------------
    // custom functions
    //-----------------------------------------------------------------------------

    float GetEventWeight();

    bool PassesCuts(mumep_ana::TrackPar_t* TrkPar) {
      if(!TrkPar)
        return false;
      bool passed = true;
      passed &= !fApplyID || TrackIDNoCRs(*TrkPar) == 0;
      passed &= fPMin < 0.f || (TrkPar->fTrack && TrkPar->fTrack->P() > fPMin);
      passed &= fPMax < 0.f || (TrkPar->fTrack && TrkPar->fTrack->P() < fPMax);
      return passed;
    }

    bool CompareHelices(TStnHelix* h1, TStnHelix* h2);
    TStnTrack* GetMatchingTrack(TStnHelix* h, int h_index, TStnTrackBlock* block);
    int HelixID(TStnHelix* h, HelixPar_t* hpar);
    int TrackID(const mumep_ana::TrackPar_t& tp);
    int TrackIDNoCRs(const mumep_ana::TrackPar_t& tp) {
      const int ID = TrackID(tp);
      const int id_no_crv = ID & (~(1 << kCRV)); // ID without the CRV coincidence cluster considered
      const int id_no_time = ID & (~(1 << kT0)); // ID with a looser timing cut
      const int id_no_crv_time = id_no_crv & id_no_time;
      return id_no_crv_time;
    }

    void InitTrackPar(TStnTrack* Trk, mumep_ana::TrackPar_t* TrkPar, TStnHelix* Hlx = nullptr) {
      TAnaModule::InitTrackPar(Trk, TrkPar, Hlx);
      if(!TrkPar)
        return;
      if(!Trk)
        return;

      // Use momentum as the observable
      TrkPar->fObs = Trk->fP; // + TrkPar->fApproxDpST;

      // Check for an associated calo cluster
      if(fClusterBlock && Trk->fClusterE > 0.f) {
        for(int icluster = 0; icluster < fClusterBlock->NClusters(); ++icluster) {
          auto cluster = fClusterBlock->Cluster(icluster);
          if(std::fabs(cluster->Energy() - Trk->fClusterE) < 1.f) { // FIXME: Check by cluster indices
            TrkPar->fCluster = cluster;
          }
        }
      }

      TrkPar->fCRVStubPar = MatchCRVToTrack(Trk, fCRVBlock, fCRVStubPar, kMaxCRVStubs);
      TStnTrack* ue_trk = MatchUpstreamTrack(Trk, fOfflineUeTrackBlock);
      TStnTrack* um_trk = MatchUpstreamTrack(Trk, fOfflineUmuTrackBlock);
      if(!ue_trk)
        TrkPar->fUpstreamTrack = um_trk;
      else if(!um_trk)
        TrkPar->fUpstreamTrack = ue_trk;
      else { // pick best track if both are defined
        // try to see if one passes the cuts
        const bool ue_pass = ue_trk->fFitCons > 1.e-5 && (ue_trk->TrkQual() < -1. || ue_trk->TrkQual() > 0.01);
        const bool um_pass = um_trk->fFitCons > 1.e-5 && (um_trk->TrkQual() < -1. || um_trk->TrkQual() > 0.01);
        if(ue_pass && !um_pass)
          TrkPar->fUpstreamTrack = ue_trk;
        else if(um_pass && !ue_pass)
          TrkPar->fUpstreamTrack = um_trk;
        else { // take closest in time
          const float ue_dt = TrkPar->fTrack->fT0 - ue_trk->fT0;
          const float um_dt = TrkPar->fTrack->fT0 - um_trk->fT0;
          if(ue_dt < um_dt)
            TrkPar->fUpstreamTrack = ue_trk;
          else
            TrkPar->fUpstreamTrack = um_trk;
        }
      }

      // Check for alternate fit hypotheses
      // if(fOfflineDeTrackBlock) {
      //   for(int itrk = 0; itrk < fOfflineDeTrackBlock->NTracks(); ++itrk) {
      //     TStnTrack* alt = fOfflineDeTrackBlock->Track(itrk);
      //     if(!alt || alt == Trk) continue;
      //     if(alt->fHelixIndex == Trk->fHelixIndex) {
      //       TrkPar->fAltHypotheses[TrkPar->fNAlt] = alt;
      //       ++TrkPar->fNAlt;
      //       break;
      //     }
      //   }
      // }
      if(fOfflineUeTrackBlock) {
        for(int itrk = 0; itrk < fOfflineUeTrackBlock->NTracks(); ++itrk) {
          if(TrkPar->fNAlt >= kMaxTrackFits)
            break;
          TStnTrack* alt = fOfflineUeTrackBlock->Track(itrk);
          if(!alt || alt == Trk)
            continue;
          if(!alt || alt->fHelixIndex == Trk->fHelixIndex) {
            TrkPar->fAltHypotheses[TrkPar->fNAlt] = alt;
            ++TrkPar->fNAlt;
            break;
          }
        }
      }
      // if(fOfflineDmuTrackBlock) {
      //   for(int itrk = 0; itrk < fOfflineDmuTrackBlock->NTracks(); ++itrk) {
      //     TStnTrack* alt = fOfflineDmuTrackBlock->Track(itrk);
      //     if(!alt || alt == Trk) continue;
      //     if(alt->fHelixIndex == Trk->fHelixIndex) {
      //       TrkPar->fAltHypotheses[TrkPar->fNAlt] = alt;
      //       ++TrkPar->fNAlt;
      //       break;
      //     }
      //   }
      // }
      // if(fOfflineUmuTrackBlock) {
      //   for(int itrk = 0; itrk < fOfflineUmuTrackBlock->NTracks(); ++itrk) {
      //     TStnTrack* alt = fOfflineUmuTrackBlock->Track(itrk);
      //     if(alt == Trk) continue;
      //     if(alt->fHelixIndex == Trk->fHelixIndex) {
      //       TrkPar->fAltHypotheses[TrkPar->fNAlt] = alt;
      //       ++TrkPar->fNAlt;
      //       break;
      //     }
      //   }
      // }

      // Evaluate MVAs
      fTreeData.fTrkQual_nactive = Trk->NActive();
      fTreeData.fTrkQual_activehitsfraction = Trk->NActive() * 1. / Trk->NHits();
      fTreeData.fTrkQual_nullhitsfraction = (Trk->NActive() > 0) ? Trk->NHitsAmbZero() * 1. / Trk->NActive() : -1.f;
      fTreeData.fTrkQual_activematsitesfraction = (Trk->NActive() > 0) ? Trk->NMatActive() * 1. / Trk->NActive() : -1.f;
      fTreeData.fTrkQual_fitcons = Trk->fFitCons;
      fTreeData.fTrkQual_momerr = Trk->fFitMomErr;
      fTreeData.fTrkQual_t0err = Trk->fT0Err;

      fTreeData.fTrkEP = (std::isfinite(Trk->fEp)) ? Trk->fEp : -1.f;
      fTreeData.fTrkDt = Trk->fDt;

      fTreeData.fTrkFitCon = (Trk->fFitCons > 0.) ? Trk->fFitCons : -100.f;
      fTreeData.fTrkLogFitCon = (Trk->fFitCons > 0.) ? log10(Trk->fFitCons) : -100.f;
      fTreeData.fTrkActiveRatio = Trk->NActive() * 1.f / Trk->NHits();
      fTreeData.fTrkNullRatio = Trk->NHitsAmbZero() * 1.f / Trk->NHits();
      fTreeData.fTrkTZSlope = TrkPar->fTZSlope;
      fTreeData.fTrkTZSlopeSig = TrkPar->TZSlopeSig();
      fTreeData.fTrkTZSlopeRatio = TrkPar->TZSlopeRatio();

      fTreeData.fTrkD0 = Trk->fD0;
      fTreeData.fTrkTanDip = Trk->fTanDip;
      fTreeData.fTrkCosTheta = TrkPar->CosTheta();
      fTreeData.fTrkRMax = TrkPar->fRMax;

      ValidateVariable(fTreeData.fTrkQual_nactive, "TrkQual_nactive");
      ValidateVariable(fTreeData.fTrkQual_activehitsfraction, "TrkQual_activehitsfraction");
      ValidateVariable(fTreeData.fTrkQual_nullhitsfraction, "TrkQual_nullhitsfraction");
      ValidateVariable(fTreeData.fTrkQual_activematsitesfraction, "TrkQual_activematsitesfraction");
      ValidateVariable(fTreeData.fTrkQual_fitcons, "TrkQual_fitcons");
      ValidateVariable(fTreeData.fTrkQual_momerr, "TrkQual_momerr");
      ValidateVariable(fTreeData.fTrkQual_t0err, "TrkQual_t0err");
      ValidateVariable(fTreeData.fTrkEP, "TrkEP");
      ValidateVariable(fTreeData.fTrkDt, "TrkDt");
      ValidateVariable(fTreeData.fTrkFitCon, "TrkFitCon");
      ValidateVariable(fTreeData.fTrkLogFitCon, "TrkLogFitCon");
      ValidateVariable(fTreeData.fTrkActiveRatio, "TrkActiveRatio");
      ValidateVariable(fTreeData.fTrkNullRatio, "TrkNullRatio");
      ValidateVariable(fTreeData.fTrkTZSlope, "TrkTZSlope");
      ValidateVariable(fTreeData.fTrkTZSlopeSig, "TrkTZSlopeSig");
      ValidateVariable(fTreeData.fTrkTZSlopeRatio, "TrkTZSlopeRatio");
      ValidateVariable(fTreeData.fTrkD0, "TrkD0");
      ValidateVariable(fTreeData.fTrkTanDip, "TrkTanDip");
      ValidateVariable(fTreeData.fTrkCosTheta, "TrkCosTheta");
      ValidateVariable(fTreeData.fTrkRMax, "TrkRMax");

      if(fEvaluateMVAs) {
        fWatch->SetTime("MVAs");
        TrkPar->fTrkQual = (fTrkQual) ? fTrkQual->EvaluateMVA("TrkQual") : -999.f;
        TrkPar->fPID = (fPID) ? fPID->EvaluateMVA("PID") : -999.f;
        TrkPar->fTrkPID = (fTrkPID) ? fTrkPID->EvaluateMVA("TrkPID") : -999.f;
        TrkPar->fCosmicID = (fCosmicID) ? fCosmicID->EvaluateMVA("CosmicID") : -999.f;
        TrkPar->fOfflinePID = -999.f;
        if(fOfflinePID) {
          // Only evaluate if the features are valid
          if(Trk->fClusterE > 0.) {
            float p_closest = Trk->fPTrackerExit;
            if(p_closest <= 0.) p_closest = Trk->fPTrackerMiddle;
            if(p_closest <= 0.) p_closest = Trk->fP;
            const float de = Trk->fClusterE - p_closest;
            const float dt = fTreeData.fTrkDt + 1.15; // add offset
            const float r  = 0.; // Radial position of the cluster
            // const float dx = Trk->fDx;
            // const float dy = Trk->fDy;
            // const float dr = std::sqrt(dx*dx + dy*dy);
            const float cz = 1.; // FIXME: Can't get cos(theta) without cluster position and 3D momentum vector
            if(de < 5.f) { // E > P + 5 = cosmic
              std::array<float,4> features{de, r, cz, dt};
              const auto mvaout = fOfflinePID->infer(features.data());
              TrkPar->fOfflinePID = mvaout[0];
            }
          }
        }
        fWatch->StopTime("MVAs");
      }

      // Apply the IDs
      TrkPar->fIDWord[0] = TrackID(*TrkPar);

      // MC CRV info
      TrkPar->fMCCRVStubPar = fMCCRVStubPar;

      // Try to find the associated sim particle
      if(fSimpBlock && Trk) {
        for(int isimp = 0; isimp < fSimpBlock->NParticles(); ++isimp) {
          const auto simp = fSimpBlock->Particle(isimp);
          if(Trk->SimID() == int(simp->GetUniqueID())) {
            TrkPar->fGenE = simp->fStartMom.E();
            TrkPar->fSimp = simp;
            break;
          }
        }
      }
    }

    void ValidateVariable(float var, const char* name) {
      if(!std::isfinite(var)) {
        auto event = GetEvent();
        printf(">>> Event %5i/%5i/%6i: Variable %s is non-finite = %f\n", event->fRunNumber, event->fSectionNumber, event->fEventNumber, name, var);
      }
    }

    void Debug();
  };
} // namespace mumep_ana
#endif
