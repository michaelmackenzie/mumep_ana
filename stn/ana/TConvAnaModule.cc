//------------------------------------------------------------------------------
// analysis module used to study mu- --> e+ analysis
//-----------------------------------------------------------------------------
#include "mumep_ana/stn/ana/TConvAnaModule.hh"

namespace mumep_ana {

  //-----------------------------------------------------------------------------
  TConvAnaModule::TConvAnaModule(const char* name, const char* title) : TAnaModule(name, title), fRand(90) {
    fAprHelixBlockName = "HelixBlockAprHighP";
    fCprDeHelixBlockName = "HelixBlockCprDeHighP";
    fOfflineDeHelixBlockName = "HelixBlockDe";
    fOfflineUeHelixBlockName = "HelixBlockUe";
    fAprTrackBlockName = "TrackBlockAprHighP";
    fCprTrackBlockName = "TrackBlockCprDeHighP";
    fOfflineDeTrackBlockName = "TrackBlockDe";
    fOfflineUeTrackBlockName = "TrackBlockUe";
    fOfflineDmuTrackBlockName = "TrackBlockDmu";
    fOfflineUmuTrackBlockName = "TrackBlockUmu";
    fCRVBlockName = "CrvClusterBlock";
    fTriggerBlockName = "TriggerBlock";
    fGenpBlockName = "GenpBlock";
    fSimpBlockName = "SimpBlock";
    fClusterBlockName = "ClusterBlock";

    // RMC model info
    fKMax = 90.1f;
    fKinematicLimit = 101.866f;
    fRMCSpectra = nullptr;
    fSpectrum = kClosure;
    fIntSpectrum = 0;

    // Additional fields
    fDebugLevel = 0;
  }

  //-----------------------------------------------------------------------------
  TConvAnaModule::~TConvAnaModule() {}

  //-----------------------------------------------------------------------------
  int TConvAnaModule::BeginJob() {

    //-----------------------------------------------------------------------------
    // register data blocks
    //-----------------------------------------------------------------------------
    RegisterDataBlock(fAprHelixBlockName, "TStnHelixBlock", &fAprHelixBlock);
    RegisterDataBlock(fCprDeHelixBlockName, "TStnHelixBlock", &fCprDeHelixBlock);
    RegisterDataBlock(fOfflineDeHelixBlockName, "TStnHelixBlock", &fOfflineDeHelixBlock);
    RegisterDataBlock(fOfflineUeHelixBlockName, "TStnHelixBlock", &fOfflineUeHelixBlock);
    RegisterDataBlock(fAprTrackBlockName, "TStnTrackBlock", &fAprTrackBlock);
    RegisterDataBlock(fCprTrackBlockName, "TStnTrackBlock", &fCprTrackBlock);
    RegisterDataBlock(fOfflineDeTrackBlockName, "TStnTrackBlock", &fOfflineDeTrackBlock);
    RegisterDataBlock(fOfflineUeTrackBlockName, "TStnTrackBlock", &fOfflineUeTrackBlock);
    RegisterDataBlock(fOfflineDmuTrackBlockName, "TStnTrackBlock", &fOfflineDmuTrackBlock);
    RegisterDataBlock(fOfflineUmuTrackBlockName, "TStnTrackBlock", &fOfflineUmuTrackBlock);
    RegisterDataBlock(fCRVBlockName, "TCrvClusterBlock", &fCRVBlock);
    RegisterDataBlock(fClusterBlockName, "TStnClusterBlock", &fClusterBlock);
    RegisterDataBlock(fTriggerBlockName, "TStnTriggerBlock", &fTriggerBlock);
    RegisterDataBlock(fGenpBlockName, "TGenpBlock", &fGenpBlock);
    RegisterDataBlock(fSimpBlockName, "TSimpBlock", &fSimpBlock);

    //-----------------------------------------------------------------------------
    // book histograms/normalization tree
    //-----------------------------------------------------------------------------
    BookHistograms();
    BookNormTree();

    // Add the CRV data pointers
    fCosmicVetoData.fHelixBlockDe = fOfflineDeHelixBlock;
    fCosmicVetoData.fHelixBlockUe = fOfflineUeHelixBlock;
    fCosmicVetoData.fTrackBlockDe = fOfflineDeTrackBlock;
    fCosmicVetoData.fTrackBlockUe = fOfflineUeTrackBlock;

    // Initialize the MVA models

    fTrkQual = new TMVA::Reader("!Color:!Silent");
    MVATools::InitializeVariables(*fTrkQual, "TrkQual", fTreeData, fTrkQualVersion);
    fTrkQual->BookMVA("TrkQual", "mumep_ana/data/trkqual_MLP.weights.xml");

    fPID = new TMVA::Reader("!Color:!Silent");
    MVATools::InitializeVariables(*fPID, "PID", fTreeData, fPIDVersion);
    fPID->BookMVA("PID", "mumep_ana/data/pid_MLP.weights.xml");

    fTrkPID = new TMVA::Reader("!Color:!Silent");
    MVATools::InitializeVariables(*fTrkPID, "TrkPID", fTreeData, fTrkPIDVersion);
    fTrkPID->BookMVA("TrkPID", "mumep_ana/data/trkpid_MLP.weights.xml");

    fCosmicID = new TMVA::Reader("!Color:!Silent");
    MVATools::InitializeVariables(*fCosmicID, "CosmicID", fTreeData, fCosmicIDVersion);
    fCosmicID->BookMVA("CosmicID", "mumep_ana/data/cosmicid_MLP.weights.xml");

    // fOfflinePID = new TMVA_SOFIE_TrackPID_v1::Session("ArtAnalysis/TrkDiag/data/TrackPID_v1.dat");

    return 0;
  }

  //-----------------------------------------------------------------------------
  void TConvAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) { TAnaModule::BookEventHistograms(Hist, Folder); }

  //-----------------------------------------------------------------------------
  void TConvAnaModule::BookHelixHistograms(HelixHist_t* Hist, const char* Folder) { TAnaModule::BookHelixHistograms(Hist, Folder); }

  //-----------------------------------------------------------------------------
  void TConvAnaModule::BookSystematicHistograms(SysHist_t* Hist, const char* Folder) {
    if(!fIncludeSys)
      return;
    if(!Hist) {
      std::cout << "TConvAnaModule::" << __func__ << ": Undefined histogram object in Folder " << Folder << std::endl;
      return;
    }
    for(int isys = 0; isys < kMaxSystematics; ++isys) {
      // check if the systematic is defined
      TString name = fSystematics.GetName(isys);
      if(name == "")
        continue;
      HBook1F(Hist->fObs[isys], Form("obs_%i", isys), Form("%s: Systematic %s", Folder, name.Data()), 300, 80., 110., Folder); // FIXME: This should inherit from the nominal observable binning
      // For debug investigations
      if(fFillVerboseSys) {
        HBook1F(Hist->fDeltaObs[isys], Form("delta_obs_%i", isys), Form("%s: Systematic %s: #DeltaObs", Folder, name.Data()), 100, -2., 2., Folder);
        HBook1F(Hist->fWeight[isys], Form("weight_%i", isys), Form("%s: Systematic %s: Weight", Folder, name.Data()), 100, 0., 2., Folder);
        HBook1F(Hist->fDeltaWeight[isys], Form("delta_weight_%i", isys), Form("%s: Systematic %s: #DeltaWt/Wt", Folder, name.Data()), 100, -1., 1., Folder);
      }
    }
  }

  //-----------------------------------------------------------------------------
  void TConvAnaModule::BookTree(TTree* tree) {
    if(!tree) {
      return;
    }
    tree->Branch("trkp", &fTreeData.fTrkP);
    tree->Branch("trkt0", &fTreeData.fTrkT0);
    tree->Branch("trkd0", &fTreeData.fTrkD0);
    tree->Branch("trktandip", &fTreeData.fTrkTanDip);
    tree->Branch("trkcostheta", &fTreeData.fTrkCosTheta);
    tree->Branch("trkfitcon", &fTreeData.fTrkFitCon);
    tree->Branch("trklogfitcon", &fTreeData.fTrkLogFitCon);
    tree->Branch("trkrmax", &fTreeData.fTrkRMax);
    tree->Branch("trkcluster", &fTreeData.fTrkCluster);
    tree->Branch("trkep", &fTreeData.fTrkEP);
    tree->Branch("trkdt", &fTreeData.fTrkDt);
    tree->Branch("trkusdt", &fTreeData.fTrkUsDt);
    tree->Branch("trkactiveratio", &fTreeData.fTrkActiveRatio);
    tree->Branch("trknullratio", &fTreeData.fTrkNullRatio);
    tree->Branch("trktzslope", &fTreeData.fTrkTZSlope);
    tree->Branch("trktzslopesig", &fTreeData.fTrkTZSlopeSig);
    tree->Branch("trktzsloperatio", &fTreeData.fTrkTZSlopeRatio);
    tree->Branch("trkpexitdiff", &fTreeData.fTrkPExitDiff);
    tree->Branch("trkqual", &fTreeData.fTrkQual);
    tree->Branch("trkpid", &fTreeData.fTrkPID);
    tree->Branch("trkonlypid", &fTreeData.fTrkOnlyPID);
    tree->Branch("trkcsmid", &fTreeData.fTrkCosmicID);
    tree->Branch("trkq", &fTreeData.fTrkCharge);
    tree->Branch("trkmcdp", &fTreeData.fTrkMCDp);

    tree->Branch("crvz", &fTreeData.fCRVZ);
    tree->Branch("crvdeltat", &fTreeData.fCRVDeltaT);
    tree->Branch("crvnpulses", &fTreeData.fCRVNPulses);
    tree->Branch("crvnpe", &fTreeData.fCRVNPe);
    tree->Branch("crvnpepp", &fTreeData.fCRVNPePP);

    tree->Branch("weight", &fEvtPar.fWeight);

    tree->Branch("train", &fTreeData.fTrain);
  }

  //-----------------------------------------------------------------------------
  void TConvAnaModule::FillEventHistograms(EventHist_t* Hist, EventPar_t* EvtPar, float Weight) { TAnaModule::FillEventHistograms(Hist, EvtPar, Weight); }

  //-----------------------------------------------------------------------------
  void TConvAnaModule::FillHelixHistograms(HelixHist_t* Hist, HelixPar_t* HlxPar, float Weight) { TAnaModule::FillHelixHistograms(Hist, HlxPar, Weight); }

  //-----------------------------------------------------------------------------
  void TConvAnaModule::FillSystematicHistograms(SysHist_t* Hist, TrackPar_t* TrkPar, float Weight) {
    if(!fIncludeSys)
      return;
    if(!Hist) {
      throw std::runtime_error(Form("TConvAnaModule::%s: Uninitialized Histogram set\n", __func__));
    }
    if(!TrkPar) {
      throw std::runtime_error(Form("TConvAnaModule::%s: Uninitialized TrkPar\n", __func__));
    }
    if(!TrkPar->fTrack) {
      throw std::runtime_error(Form("TConvAnaModule::%s: Uninitialized Track\n", __func__));
    }

    // Store parameters that may be changed
    const float o_trk_p(TrkPar->fTrack->P());
    const bool o_passed(PassesCuts(TrkPar));

    for(int isys = 0; isys < kMaxSystematics; ++isys) {
      // check if the systematic is defined
      TString name = fSystematics.GetName(isys);
      if(name == "")
        continue;
      const bool is_up = fSystematics.IsUp(isys);
      float weight(Weight); // start with the nominal weight
      bool changed = false;
      if(name == "Scale") {
        TrkPar->fTrack->fP = (is_up) ? TrkPar->fPUp : TrkPar->fPDown;
        changed = true;
      }

      // Re-evaluate variables given shifted kinematics
      if(changed) {
        InitTrackPar(TrkPar->fTrack, TrkPar, TrkPar->fHelix);
        fCluster = fTrkPar.fCluster;
        InitClusterPar(fTrkPar.fCluster, &fClusterPar);
      }

      // Re-apply the ID after systematic variation
      const bool passed = (changed) ? PassesCuts(TrkPar) : o_passed;
      // Fill the systematically varied histograms
      if(passed) {
        Hist->fObs[isys]->Fill(TrkPar->fObs, weight);
        if(fFillVerboseSys) {
          const float delta_obs = TrkPar->fTrack->fP - o_trk_p;
          const float delta_weight = (Weight != 0.) ? (weight - Weight) / Weight : 0.f;
          Hist->fDeltaObs[isys]->Fill(delta_obs, weight);
          Hist->fWeight[isys]->Fill(weight, weight); // FIXME: Should this be weighted or unweighted?
          Hist->fDeltaWeight[isys]->Fill(delta_weight, weight);
        }
      }

      // Restore variables after the shifts
      if(changed) {
        TrkPar->fTrack->fP = o_trk_p;
        InitTrackPar(TrkPar->fTrack, TrkPar, TrkPar->fHelix);
        fCluster = fTrkPar.fCluster;
        InitClusterPar(fTrkPar.fCluster, &fClusterPar);
      }
    }
  }

  //-----------------------------------------------------------------------------
  void TConvAnaModule::FillTree(TTree* tree) {
    if(!tree)
      return;
    if(!fTrack)
      return;

    // Initialize tree variables
    fTreeData.fTrkP = fTrack->fP;
    fTreeData.fTrkT0 = fTrack->fT0;
    fTreeData.fTrkD0 = fTrack->fD0;
    fTreeData.fTrkTanDip = fTrack->fTanDip;
    fTreeData.fTrkCosTheta = fTrkPar.CosTheta();
    fTreeData.fTrkFitCon = fTrack->fFitCons;
    fTreeData.fTrkLogFitCon = (fTrack->fFitCons > 0.) ? log10(fTrack->fFitCons) : -100.f;
    fTreeData.fTrkRMax = fTrkPar.fRMax;
    fTreeData.fTrkCluster = fTrack->fClusterE;
    fTreeData.fTrkEP = fTrack->fEp;
    fTreeData.fTrkDt = fTrack->fDt;
    fTreeData.fTrkUsDt = (fTrkPar.fUpstreamTrack) ? fTrack->fT0 - fTrkPar.fUpstreamTrack->fT0 : 0.;
    fTreeData.fTrkActiveRatio = fTrack->NActive() * 1.f / fTrack->NHits();
    fTreeData.fTrkNullRatio = fTrack->NHitsAmbZero() * 1.f / fTrack->NHits();
    fTreeData.fTrkTZSlope = fTrkPar.fTZSlope;
    fTreeData.fTrkTZSlopeSig = fTrkPar.TZSlopeSig();
    fTreeData.fTrkTZSlopeRatio = fTrkPar.TZSlopeRatio();
    fTreeData.fTrkPExitDiff = fTrack->fP - fTrack->fPTrackerExit;
    fTreeData.fTrkQual = fTrkPar.fTrkQual;
    fTreeData.fTrkPID = fTrkPar.fPID;
    fTreeData.fTrkOnlyPID = fTrkPar.fTrkPID;
    fTreeData.fTrkCosmicID = fTrkPar.fCosmicID;
    fTreeData.fTrkCharge = fTrack->Charge();
    fTreeData.fTrkMCDp = fTrack->fP - fTrack->fPFront;
    fTreeData.fTrkMCPDG = fTrack->PDGCode();

    if(fTrkPar.fCRVStubPar) {
      fTreeData.fCRVZ = fTrkPar.fCRVStubPar->fZ;
      fTreeData.fCRVDeltaT = fTrkPar.fCRVStubPar->fCorrTime - fTrack->fT0;
      fTreeData.fCRVNPulses = fTrkPar.fCRVStubPar->fCluster->NPulses();
      fTreeData.fCRVNPe = fTrkPar.fCRVStubPar->fCluster->NPe();
      fTreeData.fCRVNPePP = fTrkPar.fCRVStubPar->fNPePP;
    } else {
      fTreeData.fCRVZ = 0.f;
      fTreeData.fCRVDeltaT = 0.f;
      fTreeData.fCRVNPulses = 0.f;
      fTreeData.fCRVNPe = 0.f;
      fTreeData.fCRVNPePP = 0.f;
    }

    fTreeData.fWeight = fEvtPar.fWeight;

    // FIXME: For now just doing 50% splitting by ID
    const int event = GetEvent()->fEventNumber;
    fTreeData.fTrain = (event % 2 == 0) ? 1.f : -1.f;

    if(fDebugLevel > 1 || GetDebugBit(8)) {
      printf("TConvAnaModule::%s: Tree info:\n", __func__);
      const int nwords = int(((float*)&fTreeData.fTrain) - &fTreeData.fTrkP) + 1;
      float* words = &fTreeData.fTrkP;
      for(int iword = 0; iword < nwords; ++iword)
        printf(" %f\n", words[iword]);
    }
    // Fill the tree
    tree->Fill();
  }

  //_____________________________________________________________________________
  void TConvAnaModule::BookHistograms() {

    TFolder* fol;
    TFolder* hist_folder;
    char folder_name[200];
    const char* folder_title;

    DeleteHistograms();
    hist_folder = (TFolder*)GetFolder()->FindObject("Hist");

    //-----------------------------------------------------------------------------
    // book histogram selections
    //-----------------------------------------------------------------------------
    struct hist_info_t {
      TString _dsc;
      bool _trk;
      bool _hlx;
      bool _smp;
      bool _gnp;
      bool _crv;
      bool _cls;
      bool _sys;
      bool _crs; // control regions included
      bool _trs; // trees
      hist_info_t(TString dsc = "", bool trk = false, bool hlx = false, bool smp = false, bool gnp = false, bool crv = false, bool cls = false, bool sys = false, bool crs = false, bool trs = false)
          : _dsc(dsc), _trk(trk), _hlx(hlx), _smp(smp), _gnp(gnp), _crv(crv), _cls(cls), _sys(sys), _crs(crs), _trs(trs) {}
    };

    hist_info_t* hist_sets[kNHistSets];
    for(int i = 0; i < kNHistSets; i++) {
      hist_sets[i] = nullptr;
    }

    // clang-format off
    //                                 description                         trk    hlx    simp   genp   crv    cls    sys    crs    trs
    hist_sets[  0] = new hist_info_t("All events, Offline track"        ,  true,  true,  true,  true,  true,  true, false, false, false);
    hist_sets[  1] = new hist_info_t("All events, APR track"            ,  true,  true,  true,  true, false, false, false, false, false);
    hist_sets[  2] = new hist_info_t("All events, CPR track"            ,  true,  true,  true,  true, false, false, false, false, false);
    hist_sets[  3] = new hist_info_t("No weights, Offline track"        ,  true,  true,  true,  true, false, false, false, false, false);
    hist_sets[  4] = new hist_info_t("e- and p > 95"                    ,  true, false,  true, false, false, false, false, false, false);
    hist_sets[  5] = new hist_info_t("e+ and p > 75"                    ,  true, false,  true, false, false, false, false, false, false);
    hist_sets[  6] = new hist_info_t("Gen(E) > 95"                      ,  true, false,  true, false, false, false, false, false, false);
    hist_sets[  7] = new hist_info_t("e+-, TrkID, no MC cut"            ,  true, false,  true, false, false, false, false, false, false);
    hist_sets[  8] = new hist_info_t("Events with a track"              ,  true, false,  true, false, false, false, false, false, false);
    hist_sets[  9] = new hist_info_t("Track + MC cuts"                  ,  true, false,  true, false, false, false, false, false, false);
    hist_sets[ 10] = new hist_info_t("Offline, event selection"         ,  true, false,  true, false, false,  true, false,  true, false);
    hist_sets[ 11] = new hist_info_t("e+/-: trigger"                    ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 12] = new hist_info_t("e+/-: failed trigger"             ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 13] = new hist_info_t("e+/-: negative"                   ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 14] = new hist_info_t("e+/-: positive"                   ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 15] = new hist_info_t("e+/-: no weights"                 ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 16] = new hist_info_t("e-: p > 95"                       ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 17] = new hist_info_t("e-: p > 95, no ECL"               ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 18] = new hist_info_t("e+/-: dp > 1"                     ,  true, false, false, false, false, false, false, false, false);
    hist_sets[ 20] = new hist_info_t("e-: full window"                  ,  true,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 21] = new hist_info_t("e-: narrow window"                ,  true,  true,  true,  true, false,  true, false,  true, false);
    hist_sets[ 22] = new hist_info_t("e-: full window, no weights"      ,  true,  true,  true,  true, false, false, false,  true, false);
    hist_sets[ 23] = new hist_info_t("e-: high error"                   ,  true,  true,  true,  true, false, false, false,  true, false);
    hist_sets[ 24] = new hist_info_t("e-:  alt ID set"                  ,  true,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 25] = new hist_info_t("e-: !alt ID set"                  ,  true,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 27] = new hist_info_t("e-: full window, no cluster"      ,  true,  true,  true,  true,  true, false,  true,  true,  true);
    hist_sets[ 28] = new hist_info_t("e-: narrow window, no cluster"    ,  true,  true,  true,  true, false, false, false,  true, false);
    hist_sets[ 30] = new hist_info_t("e-: no CRV veto"                  ,  true, false, false, false,  true, false, false,  true, false);
    hist_sets[ 31] = new hist_info_t("e+-: no CRV veto or p cut"        ,  true, false, false, false,  true, false, false,  true, false);
    hist_sets[ 34] = new hist_info_t("e-: low dP(ST)"                   ,  true, false, false, false, false, false,  true,  true, false);
    hist_sets[ 35] = new hist_info_t("e-: high dP(ST)"                  ,  true, false, false, false, false, false,  true,  true, false);
    hist_sets[ 36] = new hist_info_t("e-: < 2 IPA int"                  ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 37] = new hist_info_t("e-: >= 2 IPA int"                 ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 38] = new hist_info_t("e-: tdip < 0.9"                   ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 39] = new hist_info_t("e-: tdip > 0.9"                   ,  true, false, false, false, false, false, false,  true, false);
    hist_sets[ 40] = new hist_info_t("e+: full window"                  ,  true,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 41] = new hist_info_t("e+: narrow window"                ,  true,  true,  true,  true, false,  true, false,  true, false);
    hist_sets[ 42] = new hist_info_t("e+: full window, no weights"      ,  true,  true,  true,  true, false, false, false,  true, false);
    hist_sets[ 50] = new hist_info_t("e+: no CRV veto"                  ,  true, false, false, false,  true, false, false,  true, false);
    hist_sets[ 60] = new hist_info_t("CRV veto + MC cut"                ,  true, false, false, false, false, false, false, false, false);

    hist_sets[ 70] = new hist_info_t("Cut-set training"                 ,  true, false, false, false, false, false, false, false,  true);
    hist_sets[ 71] = new hist_info_t("Cut-set result"                   ,  true,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 72] = new hist_info_t("Cut-set result"                   ,  true,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 73] = new hist_info_t("Cut-set result"                   ,  true,  true,  true,  true,  true,  true,  true,  true,  true);

    // CRV studies histograms
    hist_sets[ 80] = new hist_info_t("CRV: 1"                           ,  true, false, false, false,  true, false, false, false, false);
    hist_sets[ 81] = new hist_info_t("CRV: Ue/Umu tracks"               ,  true, false, false, false,  true, false, false, false, false);
    hist_sets[ 82] = new hist_info_t("CRV: calo cluster with track"     ,  true, false, false, false,  true, false, false, false, false);
    hist_sets[ 83] = new hist_info_t("CRV: dt < -60 ns"                 ,  true, false, false, false,  true, false, false, false, false);
    hist_sets[ 84] = new hist_info_t("CRV: MC electrons"                ,  true, false, false, false,  true, false, false, false,  true);
    hist_sets[ 85] = new hist_info_t("CRV: MC muons"                    ,  true, false, false, false,  true, false, false, false,  true);
    hist_sets[ 86] = new hist_info_t("CRV: Correct cluster"             ,  true, false, false, false,  true, false, false, false, false);
    hist_sets[ 87] = new hist_info_t("CRV: Upstream veto"               ,  true, false, false, false,  true, false, false, false, false);
    hist_sets[ 88] = new hist_info_t("CRV: MC downstream"               ,  true, false, false, false,  true, false, false, false, false);
    hist_sets[ 89] = new hist_info_t("CRV: MC upstream"                 ,  true, false, false, false,  true, false, false, false, false);

    // PID histograms
    hist_sets[150] = new hist_info_t("PID: ds e- + cluster"             ,  true, false, false, false, false, false, false, false, false);
    hist_sets[151] = new hist_info_t("PID: ds mu- + cluster"            ,  true, false, false, false, false, false, false, false, false);
    hist_sets[160] = new hist_info_t("PID: ds e- no cluster"            ,  true, false, false, false, false, false, false, false, fDoPIDTrees);
    hist_sets[161] = new hist_info_t("PID: ds mu- no cluster"           ,  true, false, false, false, false, false, false, false, fDoPIDTrees);
    hist_sets[170] = new hist_info_t("PID: all ds e-"                   ,  true, false, false, false, false, false, false, false, fDoPIDTrees);
    hist_sets[171] = new hist_info_t("PID: all ds mu-"                  ,  true, false, false, false, false, false, false, false, fDoPIDTrees);

    // Cosmic vs. ST-origin histograms
    hist_sets[190] = new hist_info_t("CsmID: e+-"                       ,  true, false, false, false, false, false, false, false, fDoCosmicIDTrees);

    // By process histograms
    hist_sets[200] = new hist_info_t("All DIO"                          ,  true, false, false, false, false, false, false, false, false);
    hist_sets[201] = new hist_info_t("All Cosmic"                       ,  true, false, false, false, false, false, false, false, false);
    hist_sets[202] = new hist_info_t("All RPC-ext"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[203] = new hist_info_t("All RPC-int"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[204] = new hist_info_t("All RMC-ext"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[205] = new hist_info_t("All RMC-int"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[206] = new hist_info_t("All CE"                           ,  true, false, false, false, false, false, false, false, false);
    hist_sets[207] = new hist_info_t("All IPA DIO"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[210] = new hist_info_t("ID: DIO"                          ,  true, false, false, false, false, false, false, false, false);
    hist_sets[211] = new hist_info_t("ID: Cosmic"                       ,  true, false, false, false, false, false, false, false, false);
    hist_sets[212] = new hist_info_t("ID: RPC-ext"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[213] = new hist_info_t("ID: RPC-int"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[214] = new hist_info_t("ID: RMC-ext"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[215] = new hist_info_t("ID: RMC-int"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[216] = new hist_info_t("ID: CE"                           ,  true, false, false, false, false, false, false, false, false);
    hist_sets[217] = new hist_info_t("ID: IPA DIO"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[220] = new hist_info_t("MC: DIO"                          ,  true, false, false, false, false, false, false, false, false);
    hist_sets[221] = new hist_info_t("MC: Cosmic"                       ,  true, false, false, false, false, false, false, false, false);
    hist_sets[222] = new hist_info_t("MC: RPC-ext"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[223] = new hist_info_t("MC: RPC-int"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[224] = new hist_info_t("MC: RMC-ext"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[225] = new hist_info_t("MC: RMC-int"                      ,  true, false, false, false, false, false, false, false, false);
    hist_sets[226] = new hist_info_t("MC: CE"                           ,  true, false, false, false, false, false, false, false, false);
    hist_sets[227] = new hist_info_t("MC: IPA DIO"                      ,  true, false, false, false, false, false, false, false, false);
    // clang-format on

    for(int i = 0; i < kNHistSets; i++) {
      if(hist_sets[i % 1000]) {
        if(i >= 1000 && !hist_sets[i % 1000]->_crs)
          continue;
        if(i >= 4000)
          break; // Control regions above 4000 not yet implemented
        TString desc = hist_sets[i % 1000]->_dsc;
        //-----------------------------------------------------------------------------
        // book event histograms
        //-----------------------------------------------------------------------------
        sprintf(folder_name, "evt_%i", i);
        fol = (TFolder*)hist_folder->FindObject(folder_name);
        folder_title = desc.Data();
        if(!fol)
          fol = hist_folder->AddFolder(folder_name, folder_title);
        fHist.fEvent[i] = new EventHist_t;
        BookEventHistograms(fHist.fEvent[i], Form("Hist/%s", folder_name));
        //-----------------------------------------------------------------------------
        // book genp histograms
        //-----------------------------------------------------------------------------
        if(hist_sets[i % 1000]->_gnp) {
          sprintf(folder_name, "gen_%i", i);
          fol = (TFolder*)hist_folder->FindObject(folder_name);
          folder_title = desc.Data();
          if(!fol)
            fol = hist_folder->AddFolder(folder_name, folder_title);
          fHist.fGenp[i] = new GenpHist_t;
          BookGenpHistograms(fHist.fGenp[i], Form("Hist/%s", folder_name));
        }
        //-----------------------------------------------------------------------------
        // book simp histograms
        //-----------------------------------------------------------------------------
        if(hist_sets[i % 1000]->_smp) {
          sprintf(folder_name, "sim_%i", i);
          fol = (TFolder*)hist_folder->FindObject(folder_name);
          folder_title = desc.Data();
          if(!fol)
            fol = hist_folder->AddFolder(folder_name, folder_title);
          fHist.fSimp[i] = new SimpHist_t;
          BookSimpHistograms(fHist.fSimp[i], Form("Hist/%s", folder_name));
        }
        //-----------------------------------------------------------------------------
        // book track histograms
        //-----------------------------------------------------------------------------
        if(hist_sets[i % 1000]->_trk) {
          sprintf(folder_name, "trk_%i", i);
          fol = (TFolder*)hist_folder->FindObject(folder_name);
          folder_title = desc.Data();
          if(!fol)
            fol = hist_folder->AddFolder(folder_name, folder_title);
          fHist.fTrack[i] = new TrackHist_t;
          BookTrackHistograms(fHist.fTrack[i], Form("Hist/%s", folder_name));
        }
        //-----------------------------------------------------------------------------
        // book helix histograms
        //-----------------------------------------------------------------------------
        if(hist_sets[i % 1000]->_hlx) {
          sprintf(folder_name, "hlx_%i", i);
          fol = (TFolder*)hist_folder->FindObject(folder_name);
          folder_title = desc.Data();
          if(!fol)
            fol = hist_folder->AddFolder(folder_name, folder_title);
          fHist.fHelix[i] = new HelixHist_t;
          BookHelixHistograms(fHist.fHelix[i], Form("Hist/%s", folder_name));
        }
        //-----------------------------------------------------------------------------
        // book CRV histograms
        //-----------------------------------------------------------------------------
        if(hist_sets[i % 1000]->_crv) {
          sprintf(folder_name, "crv_%i", i);
          fol = (TFolder*)hist_folder->FindObject(folder_name);
          folder_title = desc.Data();
          if(!fol)
            fol = hist_folder->AddFolder(folder_name, folder_title);
          fHist.fCRV[i] = new CRVHist_t;
          BookCRVClusterHistograms(&(fHist.fCRV[i]->fClusterHist), Form("Hist/%s", folder_name));
        }
        //-----------------------------------------------------------------------------
        // book cluster histograms
        //-----------------------------------------------------------------------------
        if(hist_sets[i % 1000]->_cls) {
          sprintf(folder_name, "cls_%i", i);
          fol = (TFolder*)hist_folder->FindObject(folder_name);
          folder_title = desc.Data();
          if(!fol)
            fol = hist_folder->AddFolder(folder_name, folder_title);
          fHist.fCluster[i] = new ClusterHist_t;
          BookClusterHistograms(fHist.fCluster[i], Form("Hist/%s", folder_name));
        }
        //-----------------------------------------------------------------------------
        // book systematic histograms
        //-----------------------------------------------------------------------------
        if(fIncludeSys && hist_sets[i % 1000]->_sys) {
          sprintf(folder_name, "sys_%i", i);
          fol = (TFolder*)hist_folder->FindObject(folder_name);
          folder_title = desc.Data();
          if(!fol)
            fol = hist_folder->AddFolder(folder_name, folder_title);
          fHist.fSys[i] = new SysHist_t;
          BookSystematicHistograms(fHist.fSys[i], Form("Hist/%s", folder_name));
        }
        //-----------------------------------------------------------------------------
        // book trees
        //-----------------------------------------------------------------------------
        if(fBookTrees && hist_sets[i % 1000]->_trs) {
          sprintf(folder_name, "trs_%i", i);
          fol = (TFolder*)hist_folder->FindObject(folder_name);
          folder_title = desc.Data();
          if(!fol)
            fol = hist_folder->AddFolder(folder_name, folder_title);
          fHist.fTree[i] = new TTree("Tree", "Tree");
          fol->Add(fHist.fTree[i]);
          BookTree(fHist.fTree[i]);
        }
      }
    }
  }

  //_____________________________________________________________________________
  int TConvAnaModule::BeginRun() {
    TAnaModule::BeginRun();
    return 0;
  }

  //_____________________________________________________________________________
  void TConvAnaModule::FillAllHistograms(const int index) {
    fWatch->SetTime("FillHistogram");
    // Apply IDs to account for systematic buffers on fill calls
    const bool passed = PassesCuts(&fTrkPar);

    if(passed) {
      if(fHist.fEvent[index])
        FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
      else
        throw std::runtime_error(Form("Histogram set %i is not initialized!", index));

      // histogram all gen/sim particles
      if(fHist.fGenp[index] && fGenpBlock) {
        for(int igen = 0; igen < fGenpBlock->NParticles(); ++igen) {
          FillGenpHistograms(fHist.fGenp[index], fGenpBlock->Particle(igen), fEvtPar.fWeight);
        }
      }
      if(fHist.fSimp[index] && fSimpBlock) {
        for(int isim = 0; isim < fSimpBlock->NParticles(); ++isim) {
          FillSimpHistograms(fHist.fSimp[index], fSimpBlock->Particle(isim), fEvtPar.fWeight);
        }
      }
      // histogram the selected track/helix
      if(fTrack && fHist.fTrack[index])
        FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
      if(fHelix && fHist.fHelix[index])
        FillHelixHistograms(fHist.fHelix[index], &fHlxPar, fEvtPar.fWeight);
      if(fCluster && fHist.fCluster[index])
        FillClusterHistograms(fHist.fCluster[index], &fClusterPar, fEvtPar.fWeight);

      // histogram CRV clusters
      if(fHist.fCRV[index] && fCRVBlock && fIgnoreCRV < 1) {
        const int nstubs = std::min(fCRVBlock->NClusters(), (int)kMaxCRVStubs);

        for(int is = 0; is < nstubs; is++) {
          if(fDebugLevel > 1)
            printf("TAnaModule::%s: Index %i: Filling CRV cluster histogram for cluster %i\n", __func__, index, is);
          CRVStubPar_t* sp = &fCRVStubPar[is];
          FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), sp, fEvtPar.fWeight);
        }
      }
    }

    // systematic histograms, allow events failing the track ID
    if(fIncludeSys && fHist.fSys[index])
      FillSystematicHistograms(fHist.fSys[index], &fTrkPar, fEvtPar.fWeight);

    // Output trees
    if(passed && fHist.fTree[index])
      FillTree(fHist.fTree[index]);

    fWatch->StopTime("FillHistogram");
  }

  //_____________________________________________________________________________
  void TConvAnaModule::FillHistograms() {

    // Reset selection
    fPMin = -1.f;
    fPMax = -1.f;

    // Event selection
    const int event_id =
        (0 + ((fEvtPar.fNTracks == 0) ? 1 << 0 : 0)  // must have a track
         + ((fEvtPar.fNGoodTracks > 1) ? 1 << 1 : 0) // no more than two good tracks
        );

    // FIXME: Add in helices
    fHelix = nullptr;
    InitHelixPar(fHelix, &fHlxPar);

    // No selection, Offline track
    fTrack = fOfflineTrack;
    InitTrackPar(fTrack, &fTrkPar, nullptr);
    fCluster = fTrkPar.fCluster;
    InitClusterPar(fTrkPar.fCluster, &fClusterPar);

    // All events
    FillAllHistograms(0);
    if(fPrimary && fPrimary->fStartMom.E() > 95.f)
      FillAllHistograms(6); // relevant region for MDS studies

    // No selection, first APR track
    fTrack = (fEvtPar.fNAprTracks > 0) ? fAprTrackBlock->Track(0) : nullptr;
    InitTrackPar(fTrack, &fTrkPar, nullptr);
    fCluster = fTrkPar.fCluster;
    InitClusterPar(fTrkPar.fCluster, &fClusterPar);
    FillAllHistograms(1);

    // No selection, first CPR track
    fTrack = (fEvtPar.fNCprTracks > 0) ? fCprTrackBlock->Track(0) : nullptr;
    InitTrackPar(fTrack, &fTrkPar, nullptr);
    fCluster = fTrkPar.fCluster;
    InitClusterPar(fTrkPar.fCluster, &fClusterPar);
    FillAllHistograms(2);

    // No selection, no weights, Offline track
    const float nominal_weight(fEvtPar.fWeight);
    fEvtPar.fWeight = 1.f;
    fTrack = fOfflineTrack;
    if(fTrack && fTrack->fHelixIndex >= 0 && fTrack->fHelixIndex < fOfflineDeHelixBlock->NHelices()) {
      fHelix = fOfflineDeHelixBlock->Helix(fTrack->fHelixIndex);
      InitHelixPar(fHelix, &fHlxPar);
    } else {
      fHelix = nullptr;
      InitHelixPar(fHelix, &fHlxPar);
    }
    InitTrackPar(fTrack, &fTrkPar, fHelix);
    fCluster = fTrkPar.fCluster;
    InitClusterPar(fTrkPar.fCluster, &fClusterPar);
    FillAllHistograms(3);
    fEvtPar.fWeight = nominal_weight;

    // Fill basic histograms per Offline De track
    for(int itrk = 0; itrk < fOfflineDeTrackBlock->NTracks(); ++itrk) {
      fTrack = fOfflineDeTrackBlock->Track(itrk);
      InitTrackPar(fTrack, &fTrkPar, nullptr);

      // Removing low energy events that are only partially modeled in the MDS
      if(fTrack->Charge() < 0 && fTrack->P() > 95.f)
        FillAllHistograms(4);
      if(fTrack->Charge() > 0 && fTrack->P() > 75.f)
        FillAllHistograms(5);
    }

    // Initialize the track parameters
    fTrack = fOfflineTrack;
    if(fTrack && fTrack->fHelixIndex >= 0 && fTrack->fHelixIndex < fOfflineDeHelixBlock->NHelices()) {
      fHelix = fOfflineDeHelixBlock->Helix(fTrack->fHelixIndex);
      InitHelixPar(fHelix, &fHlxPar);
    } else {
      fHelix = nullptr;
      InitHelixPar(fHelix, &fHlxPar);
    }
    InitTrackPar(fTrack, &fTrkPar, fHelix);
    fCluster = fTrkPar.fCluster;
    InitClusterPar(fTrkPar.fCluster, &fClusterPar);

    // Track IDs
    const int ID = TrackID(fTrkPar);
    const int id_no_crv = ID & (~(1 << kCRV));         // ID without the CRV coincidence cluster considered
    const int id_no_time = ID & (~(1 << kT0));         // ID with a looser timing cut
    const int id_no_csm_id = ID & (~(1 << kCosmicID)); // ID without the cosmic ID cut
    const int id_no_crv_time = id_no_crv & id_no_time; // no CRV or time cuts
    const int id_no_mc_cut = ID & (~(1 << kMC));       // ID without the MC cuts
    const bool crv_veto = (ID & (1 << kCRV)) == 0;     // only check the CRV bit
    const bool pid_cut = (ID & (1 << kPID)) == 0;      // only check the PID bit
    const bool us_cut  = (ID & (1 << kUpstream)) == 0; // only check the Upstream veto bit
    const bool mc_cut = (ID & (1 << kMC)) == 0;        // only check the MC bit
    if(fDebugLevel > 1)
      std::cout << "[TConvAnaModule::" << __func__ << "] Event IDs:\n"
                << " ID             = " << std::hex << ID << std::dec << std::endl
                << " id_no_crv      = " << std::hex << id_no_crv << std::dec << std::endl
                << " id_no_time     = " << std::hex << id_no_time << std::dec << std::endl
                << " id_no_csm_id   = " << std::hex << id_no_csm_id << std::dec << std::endl
                << " id_no_crv_time = " << std::hex << id_no_crv_time << std::dec << std::endl
                << " id_no_mc_cut   = " << std::hex << id_no_mc_cut << std::dec << std::endl
                << " crv_veto       = " << crv_veto << std::endl
                << " pid_cut        = " << pid_cut << std::endl
                << " us_cut         = " << us_cut << std::endl
                << " mc_cut         = " << mc_cut << std::endl;

    // Offset to control regions
    int set_offset(0);
    if(ID & (1 << kCRV))
      set_offset += kCRVVetoOffset;
    if(ID & (1 << kT0))
      set_offset += kTimeCutOffset;

    // Selection parameters
    const float eminus_pmin(100.f), eminus_pmax(110.f);
    const float eplus_pmin(84.5f), eplus_pmax(97.f);
    const float sys_buffer((fIncludeSys) ? 0.5f : 0.f);


    // Fill basic histograms per Offline track
    if(fTrack) {
      // Looking at a looser momentum selection without the MC filtering
      if(event_id == 0 && id_no_mc_cut == 0 && fTrack->P() > 75.f)
        FillAllHistograms(7);

      // Events with an Offline track
      FillAllHistograms(8);

      // Only the MC cut
      if(mc_cut)
        FillAllHistograms(9);
    }

    if(event_id == 0 && id_no_crv_time == 0) { // Offline event selection
      if(GetDebugBit(12) && fTrkPar.fSimp && fTrkPar.fSimp->CreationCode() == mu2e::ProcessCode::mu2eMuonDecayAtRest) {
        auto event = GetEvent();
        printf(">>> DIO track: Event %5i/%5i/%6i: ID = %x\n", event->fRunNumber, event->fSectionNumber, event->fEventNumber, ID);
      }
      FillAllHistograms(10 + set_offset);
      if(fEvtPar.fTriggered) {              // triggered events
        FillAllHistograms(11 + set_offset); // all trigger, accepted events
        if(fTrack->Charge() < 0)
          FillAllHistograms(13 + set_offset);
        if(fTrack->Charge() > 0)
          FillAllHistograms(14 + set_offset);
        fEvtPar.fWeight = 1.f;
        FillAllHistograms(15 + set_offset); // no weights
        fEvtPar.fWeight = nominal_weight;
        if(fTrack->Charge() < 0 && fTrack->P() > 95.f) {
          if(fTrack->fClusterE > 0.)
            FillAllHistograms(16 + set_offset);
          else
            FillAllHistograms(17 + set_offset);
        }

        // High track momentum error
        const float dp = fTrack->fP - fTrack->fPFront;
        if(dp > 1.f)
          FillAllHistograms(18);

        //-------------------------------------
        // mu- --> e- selections

        if(fTrack->Charge() < 0) {
          // full momentum window
          if(fTrack->P() > eminus_pmin - sys_buffer && fTrack->P() < eminus_pmax + sys_buffer) {
            fPMin = eminus_pmin;
            fPMax = eminus_pmax;
            if(GetDebugBit(10)) {
              auto event = GetEvent();
              printf(
                  ">>> Event %5i/%5i/%6i: Accepted in e- event window, p = %.1f, ecl = %.1f, "
                  "weight "
                  "= %.3e\n",
                  event->fRunNumber, event->fSectionNumber, event->fEventNumber, fTrack->P(), fTrack->fClusterE, fEvtPar.fWeight
              );
            }

            if(fTrack->fClusterE > 0.) { // cluster associated
              FillAllHistograms(20 + set_offset);
              // narrow momentum window
              if(fTrack->P() > 103.6 && fTrack->P() < 104.9)
                FillAllHistograms(21 + set_offset);
              fEvtPar.fWeight = 1.f;
              FillAllHistograms(22 + set_offset);
              fEvtPar.fWeight = nominal_weight;
              if(fTrack->P() - fTrack->fPFront > 1.f)
                FillAllHistograms(23 + set_offset);

              // alternate ID set on top of standard ID (~75% signal eff)
              const bool alt_id = (fTrkPar.CosTheta() > 0.5 && fTrkPar.CosTheta() < 0.646805 &&
                                   fTrack->fT0 > 594.049 &&
                                   fTrkPar.fRMax > 482.031 &&
                                   fTrack->fEp < 0.960571 &&
                                   fTrkPar.fTrkPID > 0.140625);
              if(alt_id) FillAllHistograms(24 + set_offset);
              else       FillAllHistograms(25 + set_offset);


              // Categories based on the amount of material traversed between the target and the
              // tracker
              const float pst_diff = fTrack->fPSTBack - fTrack->fP;
              const bool low_material = fTrack->NIPAIntersections() < 2;
              if(pst_diff < 0.4f) {
                FillAllHistograms(34 + set_offset);
              } else {
                FillAllHistograms(35 + set_offset);
              }
              if(low_material) {
                FillAllHistograms(36 + set_offset);
              } else {
                FillAllHistograms(37 + set_offset);
              }
              if(std::fabs(fTrack->TanDip()) < 0.9) { // cosmics typically have high tan dip
                FillAllHistograms(38 + set_offset);
              } else {
                FillAllHistograms(39 + set_offset);
              }
            } else { // no calorimeter cluster
              if(fTrkPar.TZSlopeSig() < 0.) {
                auto event = GetEvent();
                printf(
                    "TConvAnaModule::%s: Event %5i:%5i:%5i: TZ slope < 0 = %.3g, sig = %.2g, ratio "
                    "= "
                    "%.2f, ID = %x or %x\n",
                    __func__, event->fRunNumber, event->fSectionNumber, event->fEventNumber, fTrkPar.fTZSlope, fTrkPar.TZSlopeSig(), fTrkPar.TZSlopeRatio(), TrackID(fTrkPar), ID
                );
              }
              FillAllHistograms(27 + set_offset);
              // narrow momentum window
              if(fTrack->P() > 103.6 && fTrack->P() < 104.9)
                FillAllHistograms(28 + set_offset);
            }
          }
          fPMin = -1.f;
          fPMax = -1.f;
        }

        //-------------------------------------
        // mu- --> e+ selections

        if(fTrack->Charge() > 0) {
          // full momentum window
          if(fTrack->P() > eplus_pmin && fTrack->P() < eplus_pmax) {
            if(GetDebugBit(11)) {
              auto event = GetEvent();
              printf(">>> Event %5i/%5i/%6i: Accepted in e+ event window, p = %.1f, weight = %.3e\n", event->fRunNumber, event->fSectionNumber, event->fEventNumber, fTrack->P(), fEvtPar.fWeight);
            }
            FillAllHistograms(40 + set_offset);
            fEvtPar.fWeight = 1.f;
            FillAllHistograms(42 + set_offset);
            fEvtPar.fWeight = nominal_weight;
          }
          // narrow momentum window
          if(fTrack->P() > 91.0 && fTrack->P() < 92.3)
            FillAllHistograms(41 + set_offset);
        }
      } else { // failed trigger
        FillAllHistograms(12 + set_offset);
      }
    }

    if(event_id == 0 && id_no_crv == 0 && fEvtPar.fTriggered) {
      if(fTrack->Charge() < 0 && fTrack->P() > 100. && fTrack->P() < 110. && fTrack->fClusterE > 0.) {
        FillAllHistograms(30);
      }
      if(fTrack->fClusterE > 0.)
        FillAllHistograms(31); // loose momentum window and no charge selection
      if(fTrack->Charge() > 0 && fTrack->P() > eplus_pmin && fTrack->P() < eplus_pmax && fTrack->fClusterE > 0.) {
        FillAllHistograms(50);
      }
    }

    if(event_id == 0 && crv_veto && mc_cut && pid_cut) { // only apply CRV veto, MC, and PID IDs
      FillAllHistograms(60);
    }

    // loose selection to train cut-sets
    // select on momentum and loose veto on cosmics/RPC
    fTrack = fOfflineTrack;
    if(fTrack) {
      InitTrackPar(fTrack, &fTrkPar, fHelix);
      if(event_id == 0 && fEvtPar.fTriggered &&
         ((fTrack->Charge() < 0 && fTrack->P() > 100. && fTrack->P() < 110.) ||
          (fTrack->Charge() > 0 && fTrack->P() >  88. && fTrack->P() <  98.))
         && fTrkPar.STBoundary() && fTrack->T0() > 500.
         && crv_veto && pid_cut && fTrack->fClusterE > 0.) {
        FillAllHistograms(70);
      }
      // Test selection using an optimized selection here
      const bool cut_opt_id = (event_id == 0 && fEvtPar.fTriggered && fTrack->Charge() < 0
                               && fTrack->P() > 100. - sys_buffer && fTrack->P() < 110. + sys_buffer
                               && fTrkPar.STBoundary() != 0
                               && pid_cut && fTrack->fClusterE > 0.
                               && us_cut
                               && fTrkPar.CosTheta() > 0.4 && fTrkPar.CosTheta() < 0.649 // start optimized cuts
                               && fTrack->T0() > 600. && fTrack->T0() < 1650.
                               && fTrack->RMax() > 400. && fTrack->RMax() < 657.8
                               && fTrack->fDt > -4. && fTrack->fDt < 2.625
                               && fTrkPar.fTrkQual > 0.148438
                               && fTrkPar.fTrkPID > 0.078125
                               && fTrkPar.fCosmicID > 0.887779);
      if(cut_opt_id) {
        fPMin = 100.;
        fPMax = 110.;
        int cut_opt_offset = 0;
        if(!crv_veto) cut_opt_offset += kCRVVetoOffset;
        FillAllHistograms(71 + cut_opt_offset);
        fPMin = -1.;
        fPMax = -1.f;
      }
      // Test selection with the Asimov significance as the metric
      const bool asm_opt_id = (event_id == 0 && fEvtPar.fTriggered && fTrack->Charge() < 0
                               && fTrack->P() > 100. - sys_buffer && fTrack->P() < 110. + sys_buffer
                               && fTrkPar.STBoundary() != 0
                               && pid_cut && fTrack->fClusterE > 0.
                               && us_cut
                               && fTrkPar.CosTheta() > 0.517 && fTrkPar.CosTheta() < 0.649 // start optimized cuts
                               && fTrack->RMax() > 482. && fTrack->RMax() < 658.
                               && fTrack->T0() > 540. && fTrack->T0() < 1650.
                               && fTrkPar.fTrkQual > 0.248047
                               && fTrkPar.fTrkPID > 0.078125
                               && fTrkPar.fCosmicID > 0.8819);
      if(asm_opt_id) {
        fPMin = 100.;
        fPMax = 110.;
        int cut_opt_offset = 0;
        if(!crv_veto) cut_opt_offset += kCRVVetoOffset;
        FillAllHistograms(72 + cut_opt_offset);
        fPMin = -1.;
        fPMax = -1.f;
      }
      // Test selection with no CosmicID
      const bool no_csm_opt_id = (event_id == 0 && fEvtPar.fTriggered && fTrack->Charge() < 0
                               && fTrack->P() > 100. - sys_buffer && fTrack->P() < 110. + sys_buffer
                               && fTrkPar.STBoundary() != 0
                               && pid_cut && fTrack->fClusterE > 0.
                               && us_cut
                               && fTrkPar.CosTheta() > 0.525 && fTrkPar.CosTheta() < 0.649 // start optimized cuts
                               && fTrack->RMax() > 482. && fTrack->RMax() < 642.6
                               && fTrack->T0() > 540. && fTrack->T0() < 1650.
                               && fTrkPar.fTrkQual > 0.235352
                               && fTrkPar.fTrkPID > 0.078125);
      if(no_csm_opt_id) {
        fPMin = 100.;
        fPMax = 110.;
        int cut_opt_offset = 0;
        if(!crv_veto) cut_opt_offset += kCRVVetoOffset;
        FillAllHistograms(73 + cut_opt_offset);
        fPMin = -1.;
        fPMax = -1.f;
      }
    }

    //-------------------------------------
    // CRV study selections

    auto crv_stub_par = fTrkPar.fCRVStubPar;
    const float trkqual = (fTrack) ? fTrack->TrkQual() : -999.f;
    const float pid = fTrkPar.fTrkPID;
    // const float trkqual = fTrkPar.fTrkQual;
    if(fTrack && crv_stub_par && fTrack->fFitCons > 1.e-5 && trkqual > 0.2 && (pid < -10.f || pid > 0.5f) && fTrack->P() > 80.0 && fTrack->P() < 130.0 && fTrack->fT0 > 600. && fTrack->fT0 < 1650.) {
      const int trk_id = TrackID(fTrkPar);
      // const int id_no_crv = trk_id & (~(1 << kCRV)); //ID without the CRV coincidence cluster
      // considered const int id_no_crv_no_us = id_no_crv & (~(1 << kUpstream)); //ID without the
      // CRV veto or Upstream veto

      int index = 80;
      FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
      FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
      FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), crv_stub_par, fEvtPar.fWeight);

      bool upstream_candidate(false);
      for(int itrk = 0; itrk < fEvtPar.fNUeTracks; ++itrk) {
        if(upstream_candidate)
          break;
        auto trk = fOfflineUeTrackBlock->Track(itrk);
        const double deltat = fTrack->fT0 - trk->fT0;
        upstream_candidate |= deltat > 50. && deltat < 300.;
      }
      for(int itrk = 0; itrk < fEvtPar.fNUmuTracks; ++itrk) {
        if(upstream_candidate)
          break;
        auto trk = fOfflineUmuTrackBlock->Track(itrk);
        const double deltat = fTrack->fT0 - trk->fT0;
        upstream_candidate |= deltat > 50. && deltat < 300.;
      }
      if(upstream_candidate) {
        index = 81;
        FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
        FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
        FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), crv_stub_par, fEvtPar.fWeight);
      }

      // track passes upstream veto
      if(!upstream_candidate) {
        if(fTrkPar.fUpstreamTrack && GetDebugBit(11)) {
          auto event = GetEvent();
          printf(
              ">>> Event %5i/%5i/%6i: Accepted CRV upstream veto event with upstream track: qual = "
              "%.4f, fitcon = %.5f, ID = %x\n",
              event->fRunNumber, event->fSectionNumber, event->fEventNumber, fTrkPar.fUpstreamTrack->TrkQual(), fTrkPar.fUpstreamTrack->fFitCons, trk_id
          );
        }
        index = 87;
        FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
        FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
        FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), crv_stub_par, fEvtPar.fWeight);
      }
      if(fTrack->fClusterE > 0.) { // cluster associated with the track
        index = 82;
        FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
        FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
        FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), crv_stub_par, fEvtPar.fWeight);
      }
      const double deltat = fTrack->fT0 - crv_stub_par->fCorrTime;
      if(deltat < -60.) {
        index = 83;
        FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
        FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
        FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), crv_stub_par, fEvtPar.fWeight);
      }
      const int mc_pdg = std::abs(fTrack->PDGCode());
      if(mc_pdg == 11) {
        index = 84;
        FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
        FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
        FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), crv_stub_par, fEvtPar.fWeight);
      }
      if(mc_pdg == 13) {
        index = 85;
        FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
        FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
        FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), crv_stub_par, fEvtPar.fWeight);
      }
      if(fMCCRVStubPar) {
        index = 86;
        fTrkPar.fCRVStubPar = fMCCRVStubPar;
        FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
        FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
        FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), fMCCRVStubPar, fEvtPar.fWeight);
        fTrkPar.fCRVStubPar = crv_stub_par;
      }
      if(fTrack->fMcDirection > 0) {
        index = 88;
        FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
        FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
        FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), crv_stub_par, fEvtPar.fWeight);
      }
      if(fTrack->fMcDirection < 0) {
        index = 89;
        FillEventHistograms(fHist.fEvent[index], &fEvtPar, fEvtPar.fWeight);
        FillTrackHistograms(fHist.fTrack[index], &fTrkPar, fEvtPar.fWeight);
        FillCRVClusterHistograms(&(fHist.fCRV[index]->fClusterHist), crv_stub_par, fEvtPar.fWeight);
      }
    }

    //-------------------------------------
    // PID selections

    if(fTrack && fTrack->P() > 75.f && fTrack->P() < 150.f) {
      InitTrackPar(fTrack, &fTrkPar, fHelix);
      fCluster = fTrkPar.fCluster;
      InitClusterPar(fTrkPar.fCluster, &fClusterPar);
      fEvtPar.fWeight = 1.f; // Ignore weights
      const int mc_pdg = std::abs(fTrack->PDGCode());
      const bool has_cluster = fTrack->fClusterE > 0.f;
      const bool mc_downstream = fTrack->fMcDirection > 0;
      if(mc_downstream) {
        const float tzslope = fTrkPar.fTZSlope;
        if(tzslope == 0.f) {
          auto event = GetEvent();
          printf(
              ">>> Missing TZSlope: Event %5i/%5i/%6i: Helix = %o, Helix slope = %.3g\n", event->fRunNumber, event->fSectionNumber, event->fEventNumber, fHelix != nullptr,
              (fHelix) ? fHelix->TZSlope() : 0.f
          );
        }
        if(has_cluster) { // PID can include calo
          if(mc_pdg == 11)
            FillAllHistograms(150);
          if(mc_pdg == 13)
            FillAllHistograms(151);
        } else {
          if(mc_pdg == 11)
            FillAllHistograms(160);
          if(mc_pdg == 13)
            FillAllHistograms(161);
        }
        if(mc_pdg == 11)
          FillAllHistograms(170);
        if(mc_pdg == 13)
          FillAllHistograms(171);
        if(GetDebugBit(11) && mc_pdg == 13 && fTrkPar.fCRVStubPar && fTrkPar.CRVSTDeltaT() < -50.) {
          auto event = GetEvent();
          printf(">>> Event %5i/%5i/%6i: Muon with trk delta t (ST) < -50 ns = %.1f\n", event->fRunNumber, event->fSectionNumber, event->fEventNumber, fTrkPar.CRVSTDeltaT());
        }
      }
      fEvtPar.fWeight = nominal_weight;
    }

    //-------------------------------------
    // Cosmic ID selections

    if(fTrack && fTrack->P() > 80.f && fTrack->P() < 110.f && fTrack->fFitCons >= 1.e-5) {
      InitTrackPar(fTrack, &fTrkPar, fHelix);
      fCluster = fTrkPar.fCluster;
      InitClusterPar(fTrkPar.fCluster, &fClusterPar);
      fEvtPar.fWeight = 1.f; // Ignore weights
      const int mc_pdg = std::abs(fTrack->PDGCode());
      // const bool has_cluster = fTrack->fClusterE > 0.f;
      // const bool mc_downstream = fTrack->fMcDirection > 0;
      // const float tzslope = fTrkPar.fTZSlope;
      if(mc_pdg == 11 && fTrack->Charge() < 0)
        FillAllHistograms(190);
      fEvtPar.fWeight = nominal_weight;
    }

    //-------------------------------------
    // By process selections
    if(fTrack && fPrimary) {
      TString type = (IsIPADIO(fPrimary)) ? "IPA DIO" : ProcessGroup(fPrimary->CreationCode());
      // No selection
      if(type == "DIO")
        FillAllHistograms(200);
      else if(type == "Primary")
        FillAllHistograms(201);
      else if(type == "Ext-RPC")
        FillAllHistograms(202);
      else if(type == "Int-RPC")
        FillAllHistograms(203);
      else if(type == "Ext-RMC")
        FillAllHistograms(204);
      else if(type == "Int-RMC")
        FillAllHistograms(205);
      else if(type.BeginsWith("CE"))
        FillAllHistograms(206);
      else if(type == "IPA DIO")
        FillAllHistograms(207);
      // Wide momentum window selection
      if(ID == 0 && ((fTrack->Charge() < 0 && fTrack->fP > eminus_pmin && fTrack->fP < eminus_pmax) || (fTrack->Charge() > 0 && fTrack->fP > eplus_pmin && fTrack->fP < eplus_pmax))) {
        if(type == "DIO")
          FillAllHistograms(210);
        else if(type == "Primary")
          FillAllHistograms(211);
        else if(type == "Ext-RPC")
          FillAllHistograms(212);
        else if(type == "Int-RPC")
          FillAllHistograms(213);
        else if(type == "Ext-RMC")
          FillAllHistograms(214);
        else if(type == "Int-RMC")
          FillAllHistograms(215);
        else if(type.BeginsWith("CE"))
          FillAllHistograms(216);
        else if(type == "IPA DIO")
          FillAllHistograms(217);
      }
      // Just the MC region selection
      if(mc_cut && ((fTrack->Charge() < 0 && fTrack->fP > 95. && fTrack->fP < 150.) || (fTrack->Charge() > 0 && fTrack->fP > 80. && fTrack->fP < 150.))) {
        if(type == "DIO")
          FillAllHistograms(220);
        else if(type == "Primary")
          FillAllHistograms(221);
        else if(type == "Ext-RPC")
          FillAllHistograms(222);
        else if(type == "Int-RPC")
          FillAllHistograms(223);
        else if(type == "Ext-RMC")
          FillAllHistograms(224);
        else if(type == "Int-RMC")
          FillAllHistograms(225);
        else if(type.BeginsWith("CE"))
          FillAllHistograms(226);
        else if(type == "IPA DIO")
          FillAllHistograms(227);
      }
    }
  }

  //-----------------------------------------------------------------------------
  int TConvAnaModule::Event(int ientry) {
    fWatch->Increment("Total");
    fWatch->SetTime("Event");
    ++fNormInfo.nseen_;

    // get entry for the ith event for each data block
    fWatch->SetTime("GetBlocks");
    fAprHelixBlock->GetEntry(ientry);
    fCprDeHelixBlock->GetEntry(ientry);
    fOfflineDeHelixBlock->GetEntry(ientry);
    fOfflineUeHelixBlock->GetEntry(ientry);
    fAprTrackBlock->GetEntry(ientry);
    fCprTrackBlock->GetEntry(ientry);
    fOfflineDeTrackBlock->GetEntry(ientry);
    fOfflineUeTrackBlock->GetEntry(ientry);
    fOfflineDmuTrackBlock->GetEntry(ientry);
    fOfflineUmuTrackBlock->GetEntry(ientry);
    fTriggerBlock->GetEntry(ientry);
    if(fClusterBlock)
      fClusterBlock->GetEntry(ientry);
    if(fGenpBlock)
      fGenpBlock->GetEntry(ientry);
    if(fSimpBlock)
      fSimpBlock->GetEntry(ientry);
    if(fIgnoreCRV < 3 && fCRVBlock)
      fCRVBlock->GetEntry(ientry);
    fWatch->StopTime("GetBlocks");

    InitEventInfo();

    // get/set event parameters
    fGen = (fGenpBlock && fGenpBlock->NParticles() > 0) ? fGenpBlock->Particle(0) : nullptr;
    fSim = (fSimpBlock && fSimpBlock->NParticles() > 0) ? fSimpBlock->Particle(0) : nullptr;

    // Find the "primary" in the event and initialize other event info
    fPrimary = nullptr;
    if(fSimpBlock) {
      // Can use IsPrimary() check
      const bool check_primary = (fSimpBlock->NParticles() > 0) ? fSimpBlock->Particle(0)->IsA()->GetClassVersion() > 3 : false;
      for(int isim = 0; isim < fSimpBlock->NParticles(); ++isim) {
        auto sim = fSimpBlock->Particle(isim);
        if(check_primary) {
          if(!sim->IsPrimary()) continue; // keep looking for the primary
          // Keep this, unless there are two primaries and the other is higher energy
          if(fPrimary == nullptr || fPrimary->fStartMom.E() < sim->fStartMom.E()) {
            fPrimary = sim;
            break;
          } else continue;
        }
        if(IsSignal(sim->CreationCode())) { // prioritize signal label over background label
          fPrimary = sim;
          break;
        }
        // Check for a cosmic
        if(sim->CreationCode() == mu2e::ProcessCode::mu2ePrimary && abs(sim->PDGCode()) == 13) {
          fPrimary = sim;
          break;
        }
        // Other backgrounds
        if(IsBackground(sim->CreationCode())) {
          fPrimary = sim;
          if(sim->CreationCode() == mu2e::ProcessCode::mu2eFlatPhoton) { // ensure we don't use the e-/e+ daughters as the
                                                                         // "primary"
            break;
          }
        }
      }
      // If failed to find the primary, try again looking for a cosmic primary (proton, neutron,
      // etc.)
      if(!fPrimary) {
        for(int isim = 0; isim < fSimpBlock->NParticles(); ++isim) {
          auto sim = fSimpBlock->Particle(isim);
          // Take the first "mu2ePrimary" seen
          if(sim->CreationCode() == mu2e::ProcessCode::mu2ePrimary) {
            fPrimary = sim;
            break;
          }
        }
      }

      // Initialize event-level MC info
      fEvtPar.fRMCEnergy = 0.;
      for(int isim = 0; isim < fSimpBlock->NParticles(); ++isim) {
        auto sim = fSimpBlock->Particle(isim);
        const int sim_code = sim->CreationCode();
        // Check for RMC photons
        if(sim_code == mu2e::ProcessCode::mu2eExternalRMC || sim_code == mu2e::ProcessCode::mu2eFlatPhoton) {
          const float gen_energy = sim->fStartMom.E();
          fEvtPar.fRMCEnergy = gen_energy;
        }
        if(sim_code == mu2e::ProcessCode::mu2eInternalRMC) {
          const float gen_energy = sim->fStartMom.E();
          fEvtPar.fRMCEnergy += gen_energy; // sum of the two photons
        }
      }
    }

    fEvtPar.fInstLum = GetHeaderBlock()->fInstLum;
    fEvtPar.fNAprHelices = fAprHelixBlock->NHelices();
    fEvtPar.fNCprHelices = fCprDeHelixBlock->NHelices();
    fEvtPar.fNOfflineHelices = fOfflineDeHelixBlock->NHelices();
    fEvtPar.fNAprTracks = fAprTrackBlock->NTracks();
    fEvtPar.fNCprTracks = fCprTrackBlock->NTracks();
    fEvtPar.fNTracks = fOfflineDeTrackBlock->NTracks();
    fEvtPar.fNUeTracks = fOfflineUeTrackBlock->NTracks();
    fEvtPar.fNDmuTracks = fOfflineDmuTrackBlock->NTracks();
    fEvtPar.fNUmuTracks = fOfflineUmuTrackBlock->NTracks();
    fEvtPar.fPassedCprPath = fTriggerBlock->PathPassed(160); // path 150 = cprDe_TrkDe_80m70p_D0200, 160 = cprDe_TrkDe_80m70p
    fEvtPar.fPassedAprPath = fTriggerBlock->PathPassed(210); // path 200 = apr_TrkDe_80m70p_D0200, 210 = apr_TrkDe_80m70p
    fEvtPar.fTriggered = fEvtPar.fPassedCprPath || fEvtPar.fPassedAprPath;

    // Initialize CRV info
    fMCCRVStubPar = nullptr;
    if(fCRVBlock && fCRVBlock->NClusters() > kMaxCRVStubs)
      throw std::runtime_error(Form("TConvAnaModule::%s: More CRV clusters than allows (%i seen)\n", __func__, fCRVBlock->NClusters()));

    fEvtPar.fNCRVClusters = fCRVBlock ? fCRVBlock->NClusters() : 0;
    if(fCRVBlock && fIgnoreCRV < 2) {
      for(int i = 0; i < kMaxCRVStubs; ++i)
        fCRVStubPar[i].reset();
      InitCRVStubPar(fCRVBlock, fCRVStubPar, kMaxCRVStubs);
    }
    // count good CRV clusters
    fEvtPar.fNGoodCRVClusters = 0;
    if(fCRVBlock && fIgnoreCRV < 1) {
      for(int icl = 0; icl < fEvtPar.fNCRVClusters; ++icl) {
        auto cl = fCRVBlock->Cluster(icl);
        if(!cl)
          continue;
        if(cl->EndTime() < 1700. && cl->StartTime() > 450.) {
          ++fEvtPar.fNGoodCRVClusters;
        }
      }
    }
    // look for the MC stub
    for(int icrv = 0; icrv < fEvtPar.fNCRVClusters; ++icrv) {
      auto stub = &fCRVStubPar[icrv];
      if(!stub->fCluster)
        continue;
      const int sim_id = stub->fCluster->fSimID;
      if(sim_id < 0)
        continue;
      if(fSimpBlock) {
        for(int isim = 0; isim < fSimpBlock->NParticles(); ++isim) {
          auto sim = fSimpBlock->Particle(isim);
          if(int(sim->GetUniqueID()) == sim_id) {
            if(fDebugLevel > 0)
              printf("TConvAnaModule::%s: Found CRV cluster for SIM ID %i\n", __func__, sim_id);
            fMCCRVStubPar = stub;
            break;
          }
        }
      }
      if(fMCCRVStubPar)
        break;
    }
    if(!fMCCRVStubPar && fDebugLevel > 0)
      printf("TConvAnaModule::%s: Did not find CRV cluster!\n", __func__);

    fEvtPar.fNGoodTracks = 0;
    fOfflineTrack = nullptr;
    for(int itrk = 0; itrk < fEvtPar.fNTracks; ++itrk) {
      fTrack = fOfflineDeTrackBlock->Track(itrk);
      if(!fOfflineTrack)
        fOfflineTrack = fTrack; // if no good track, just take a track
      if(fTrack->fHelixIndex >= 0 && fTrack->fHelixIndex < fOfflineDeHelixBlock->NHelices()) {
        fHelix = fOfflineDeHelixBlock->Helix(fTrack->fHelixIndex);
        InitHelixPar(fHelix, &fHlxPar);
      } else {
        fHelix = nullptr;
        InitHelixPar(fHelix, &fHlxPar);
      }
      InitTrackPar(fTrack, &fTrkPar, fHelix);
      fCluster = fTrkPar.fCluster;
      InitClusterPar(fTrkPar.fCluster, &fClusterPar);
      const int ID = TrackID(fTrkPar);
      if((ID & (~(1 << kMC))) == 0) { // don't include MC selections in this cut
        ++fEvtPar.fNGoodTracks;
        if(ID == 0 || !fOfflineTrack)
          fOfflineTrack = fTrack; // store a good selected track preferentially
      }
    }

    fEvtPar.fNonCRVVetoID = NonCRVCosmicVeto(&fCosmicVetoData);
    fEvtPar.fRMCEnergy = 0.f;

    fEvtPar.fWeight = GetEventWeight();
    if(GetDebugBit(9)) {
      GetHeaderBlock()->Print(Form(" Weight = %8.2g, RMC energy = %6.2f, Background type = %2i",
                                   fEvtPar.fWeight, fEvtPar.fRMCEnergy, fBackgroundType));
    }
    if(!std::isfinite(fEvtPar.fWeight) || fEvtPar.fWeight < 0.) {
      auto event = GetEvent();
      printf(">>> Undefined or negative event weight: Event %5i/%5i/%6i: weight = %f\n", event->fRunNumber, event->fSectionNumber, event->fEventNumber, fEvtPar.fWeight);
      fEvtPar.fWeight = 0.f;
    }

    Debug();

    fWatch->SetTime("Histograms");
    FillHistograms();
    fWatch->StopTime("Histograms");

    fWatch->StopTime("Event");

    return 0;
  }

  //_____________________________________________________________________________
  float TConvAnaModule::GetEventWeight() {
    float weight(1.f); // default to no weight

    // Don't apply weights to data (for now)
    if(fDataType == kData)
      return weight;

    // Identify the background type
    fBackgroundType = -1;


    //-------------------------------------------------
    // Spectrum weights

    if(fSimpBlock) {
      // Can use IsPrimary() check
      const bool check_primary = (fSimpBlock->NParticles() > 0) ? fSimpBlock->Particle(0)->IsA()->GetClassVersion() > 3 : false;

      // Determine the sample from the SIM input
      bool beam_process(false), is_rpc(false);
      const TSimParticle* pion_sim(nullptr);
      for(int isim = 0; isim < fSimpBlock->NParticles(); ++isim) {
        const auto sim = fSimpBlock->Particle(isim);
        if(std::abs(sim->PDGCode()) == 211) pion_sim = sim; // save for potential RPC weights
        const int sim_code = sim->CreationCode();
        bool inspect = !check_primary || sim->IsPrimary();
        inspect |= sim_code == mu2e::ProcessCode::mu2eAntiproton; // antiproton may not be the primary
        if(!inspect) continue;
        beam_process |= IsBeamProcess(sim_code);
        if(sim_code == mu2e::ProcessCode::mu2eFlateMinus) { // assume re-weighting to DIO
          fBackgroundType = kDIO;
          const float gen_energy = sim->fStartMom.E();
          // FIXME: Decide what order weight to use
          const float pdf = fStntuple->DioWeightAlFull(gen_energy); // TStntuple::DioWeightAl_LL(gen_energy);
          if(fDebugLevel > 0)
            printf("TConvAnaModule::%s: DIO: E = %.3f --> weight = %.3g\n", __func__, gen_energy, pdf);
          weight *= pdf;
          break;
        } else if(sim_code == mu2e::ProcessCode::mu2eFlatePlus) { // assume re-weighting to RMC
                                                                  // internal positron spectrum
          fBackgroundType = kInternalRMC;
          const float gen_energy = sim->fStartMom.E();
          // check if the spectrum weights have been initialized
          if(!fRMCSpectra) {
            std::cout << "TConvAnaModule::" << __func__ << ": Initializing RMC Spectrum weight information\n";
            if(fRMCVersion == kPlestid) fSpectrum = RMCSpectra::kPlestid;
            fRMCSpectra = new RMCSpectra(fKMax, fKinematicLimit, fSpectrum, 1, fIntSpectrum);
            fRMCSpectra->InitializeSpectrum();
            // Scale the input spectrum to be normalized above 57 MeV instead of over the entire photon spectrum
            const double scale = fRMCSpectra->fSpectrum_->Integral(0., 102.) / fRMCSpectra->fSpectrum_->Integral(57., 102.);
            fRMCSpectra->hSpectrum_->Scale(scale);
          }
          const float pdf = fRMCSpectra->Weight(gen_energy);
          if(fDebugLevel > 0)
            printf("TConvAnaModule::%s: RMC: E = %.3f --> weight = %.3g\n", __func__, gen_energy, pdf);
          weight *= pdf;
          break;
        } else if(sim_code == mu2e::ProcessCode::mu2eFlatPhoton) { // assume re-weighting to RMC
                                                                   // photon spectrum
          fBackgroundType = kExternalRMC;
          const float gen_energy = sim->fStartMom.E();
          fEvtPar.fRMCEnergy = gen_energy;
          weight *= RMCWeight(gen_energy, fRMCVersion);
          // weight *= (gen_energy > fKMax) ? 0. : TStntuple::RMC_ClosureAppxWeight(gen_energy, fKMax) / fKMax;
          break;
        } else if(sim_code == mu2e::ProcessCode::mu2eExternalRMC) { // assume physical spectrum
          fBackgroundType = kInternalRMC;
        } else if(sim_code == mu2e::ProcessCode::mu2eInternalRMC) { // assume physical spectrum
          fBackgroundType = kInternalRMC;
        } else if(sim_code == mu2e::ProcessCode::mu2eAntiproton) { // assume re-weighting the PbarSTGun
          fBackgroundType = kSimplePBar;
          const float z = sim->fStartPos.Z();
          const float time = sim->fStartPos.T();
          const float r = std::sqrt(std::pow(sim->fStartPos.X() + 3904., 2) + std::pow(sim->fStartPos.Y(), 2));
          const float pdf = PBarWeight(z, time, r);
          if(fDebugLevel > 0)
            printf("TConvAnaModule::%s: Pbar: z = %.1f, t = %.1f, r = %.1f --> weight = %.3g\n", __func__, z, time, r, pdf);
          weight *= pdf;
          break;
        } else if(sim_code == mu2e::ProcessCode::mu2eExternalRPC || sim_code == mu2e::ProcessCode::mu2eInternalRPC) {
          fBackgroundType = (sim_code == mu2e::ProcessCode::mu2eExternalRPC) ? kExternalRPC : kInternalRPC;
          is_rpc = true;
        }
      }

      // Apply RPC re-weighting if a pion is found and an RPC process was created
      if(is_rpc && pion_sim) {
        if(!fDataset.BeginsWith("rpce1") && !fDataset.BeginsWith("rpci1")) { // not physical RPC samples
          const float proper_time = pion_sim->EndProperTime(); // particle lifetime is included in the time already (time/tau)
          const float pion_weight = std::exp(-proper_time);
          if(fUsePionWeights > 0)
            weight *= pion_weight;
          else if(fUsePionWeights == -1)
            weight *= 0.179 / 83146.; // apply the average weight (to model MDS2)
          if(fDebugLevel > 0)
            printf(
                   "TConvAnaModule::%s: RPC: t_start = %.1f, t_end = %.1f, proper_time = %.1f --> "
                   "weight = %.3g\n",
                   __func__, pion_sim->StartPos()->T(), pion_sim->EndPos()->T(), proper_time, pion_weight
                   );
        }
      }

      //-------------------------------------------------
      // Beam process weights

      // Correct the intensity profile in beam-related simulations
      float beam_weight = 1.f;
      if(beam_process && fEvtPar.fInstLum > 1. && fUseBeamWeights)
        beam_weight = BeamProcessWeight(fEvtPar.fInstLum, fBatchModeSim);
      if(fDebugLevel > 2)
        printf(
            "TConvAnaModule::%s: UseBeamWeights = %i, beam_process = %o, N(POT) = %.2g, beam "
            "weight "
            "= %.3f\n",
            __func__, fUseBeamWeights, beam_process, fEvtPar.fInstLum, beam_weight
        );
      weight *= beam_weight;

      // If using the wrong profile simulation, correct this
      if(fBatchModeSim > 0 && fBatchModeUse > 0 && fBatchModeUse != fBatchModeSim) {
        if(fBatchModeSim > 2)
          throw std::runtime_error(Form("Unknown beam simulation mode %i!", fBatchModeSim));
        const int num = (fBatchModeSim == 1) ? 2 : 1;
        const int den = (fBatchModeSim == 1) ? 1 : 2;
        const float batch_weight = BatchModeWeight(fEvtPar.fInstLum, num) / BatchModeWeight(fEvtPar.fInstLum, den);
        if(!std::isfinite(batch_weight)) {
          auto event = GetEvent();
          printf(">>> Undefined batch weight weight: Event %5i/%5i/%6i: lum = %f\n", event->fRunNumber, event->fSectionNumber, event->fEventNumber, fEvtPar.fInstLum);
        } else
          weight *= batch_weight;
      }
    }

    return weight;
  }

  //_____________________________________________________________________________
  TStnTrack* TConvAnaModule::GetMatchingTrack(TStnHelix* h, int, TStnTrackBlock* block) {
    for(int t_index = 0; t_index < block->NTracks(); ++t_index) {
      auto track = block->Track(t_index);
      if(!track)
        continue;
      bool matched(true);
      matched &= std::fabs(h->D0() - track->D0()) < 20.;
      matched &= std::fabs(h->Pt() - track->Pt()) / (h->Pt() + track->Pt()) / 2. < 0.02;
      if(matched)
        return track;
    }
    return nullptr;
  }

  //_____________________________________________________________________________
  // Apply selection cuts to the helix
  int TConvAnaModule::HelixID(TStnHelix*, HelixPar_t*) {
    int ID(0);
    return ID;
  }

  //_____________________________________________________________________________
  // Apply selection cuts to a track
  int TConvAnaModule::TrackID(const mumep_ana::TrackPar_t& trkpar) {
    auto track = trkpar.fTrack;
    if(!track)
      return 0;
    const float trkqual = trkpar.fTrkQual;
    int ID(0);
    if(track->Charge() < 0) {
      if(track->P() < 85. || track->P() > 130.)
        ID += 1 << kP;
    } else {
      if(track->P() < 85. || track->P() > 130.)
        ID += 1 << kP;
    }
    if(trkpar.fRMax < 430. || trkpar.fRMax > 650.)
      ID += 1 << kRMax;
    else if(track->NOPAIntersections())
      ID += 1 << kRMax;
    if(trkqual > -10. && trkqual < 0.2)
      ID += 1 << kTrkQual;
    if(track->T0() < 600. || track->T0() > 1650.)
      ID += 1 << kT0;
    if(track->fFitCons < 1.e-5)
      ID += 1 << kFitCon;
    if(trkpar.STBoundary() == 0)
      ID += 1 << kD0; // consistent with ST
    // float d0_sign = track->fD0 * track->Charge();
    // if(d0_sign < -100. || d0_sign > 60.)                       ID += 1 << kD0; //consistent with
    // ST
    if(track->fTanDip < 0.5 || track->fTanDip > 2.0)
      ID += 1 << kTDip;
    // if(track->fTanDip < 0.5 || track->fTanDip > 1.5)           ID += 1 << kTDip;
    if(track->T0() < 500. || track->T0() > 1650.)
      ID += 1 << kT0Loose; // for control regions

    // upstream track rejection
    auto us_trk = trkpar.fUpstreamTrack;
    if(us_trk) {
      // minimal selection on the upstream track quality
      // if(us_trk->fFitCons > 1.e-5 && (us_trk->TrkQual() < -10. || us_trk->TrkQual() > 0.01))
      const double us_dt = track->T0() - us_trk->T0();
      if(us_dt > 60. && us_dt < 110.) // window of time on the upstream leg
        ID += 1 << kUpstream;
    }

    // PID requirements
    if(track->fClusterE <= 0.) {
      if(trkpar.fTrkPID < -100.f)
        ID += 1 << kClusterE; // require a cluster if no tracker-only PID
      else if(trkpar.fTrkPID < 0.5f)
        ID += 1 << kPID;
    } else {
      if(trkpar.fPID < -100.f && track->fEp < 0.65)
        ID += 1 << kPID; // strict cluster requirement if no PID
      else if(trkpar.fPID < 0.5f)
        ID += 1 << kPID;
    }

    // Kinematic cosmic ID
    if(trkpar.fCosmicID > -100.f && trkpar.fCosmicID < 0.85f)
      ID += 1 << kCosmicID;

    // alternate hypothesis (trajectory, PID) veto
    bool fail_fit_hyp = false;
    for(int ialt = 0; ialt < trkpar.fNAlt; ++ialt) {
      fail_fit_hyp = !ResolveAmbiguity(track, trkpar.fAltHypotheses[ialt]);
      if(fail_fit_hyp)
        break;
    }

    // without a cluster, add a stricter trajectory selection
    if(!fail_fit_hyp && track->fClusterE <= 0.) {
      fail_fit_hyp = trkpar.TZSlopeSig() < 0.;
    }
    if(fail_fit_hyp)
      ID += 1 << kFitHyp;

    // CRV rejection
    bool tagged(false);
    auto crv_stub_par = trkpar.fCRVStubPar;
    if(crv_stub_par) {
      // const float deltat_calo = fTrack->fT0 - crv_stub_par->fApproxTimeCaloToFront; //consider
      // both trajectory hypotheses
      const float deltat_st = trkpar.CRVSTDeltaT();                 // ST origin
      const float deltat_calo = trkpar.CRVCaloFrontDeltaT(false);   // for electron tracks from calo
      const float deltat_calo_mu = trkpar.CRVCaloFrontDeltaT(true); // for muon tracks from calo
      const float deltat_crv = fTrack->fT0 - crv_stub_par->fTime;   // effective for upstream misreconstruction
      tagged |= (deltat_st > -40.f && deltat_st < 50.f);            // ST hypothesis
      tagged |= (deltat_calo > -30.f && deltat_calo < 40.f);        // Calo hypothesis
      // tagged |= (deltat_crv > -10.f && deltat_crv < 10.f); // Upstream hypothesis
      if(track->fClusterE <= 0.) {
        tagged |= deltat_calo_mu > -40.f && deltat_calo_mu < 50.f; // Muon hypothesis if no cluster for PID
        tagged |= (deltat_crv > -10.f && deltat_crv < 10.f);       // Upstream hypothesis
      }
      if(tagged)
        ID += 1 << kCRV;

      // const float min_extrap_dt(-50.f), max_extrap_dt(60.f);
      // if((deltat_st   > min_extrap_dt && deltat_st   < max_extrap_dt) ||
      //    (deltat_calo > min_extrap_dt && deltat_calo < max_extrap_dt) ||
      //    (deltat_crv > -25.f && deltat_crv < 0.f))               ID += 1 << kCRV;
    }

    // MC selection to match MC cuts applied to the Mock Data Samples
    bool mc_selection = true;

    if(fDataType == kBackground) {
      if(fBackgroundType == kDIO && fPrimary && !fDataset.BeginsWith("ipad")) {
        mc_selection &= trkpar.fGenE <= 0. || trkpar.fGenE > 95.f; // only included DIO with E > 95 MeV
      }
    }

    // if(fPrimary
    //    && fPrimary->CreationCode() != mu2e::ProcessCode::mu2eFlatPhoton  // FIXME: don't cut on
    //    RMC photons for now
    //    && trkpar.fTrack->Charge() < 0) {  // Only apply this to negative tracks
    //   mc_selection &= trkpar.fGenE <= 0. || trkpar.fGenE > 95.f; // only included DIO with E > 95
    //   MeV
    //   // mc_selection &= trkpar.fGenE <= 0. || trkpar.fGenE > 80.f;
    // } else if(fPrimary) {
    //   mc_selection &= trkpar.fGenE <= 0. || trkpar.fGenE > 85.f; // only include positrons > 85
    //   MeV
    // }

    // MC-truth cut, don't include pileup in background-specific histogramming
    if(fDataType == kBackground) {
      bool find_match(false); // FIXME: Add MC relation field to TSimParticle to check track particle relation
      switch(fBackgroundType) {
      case kInternalRPC:
      case kExternalRPC:
      case kInternalRMC:
      case kExternalRMC:
      case kDIO:
        find_match = true;
        break;
      default:
        break;
      }
      if(find_match) {
        auto sim = trkpar.fSimp;
        const bool found = sim != nullptr;
        mc_selection &= found;
        if(found) {
          // check the found sim is consistent with the background type
          if(fBackgroundType != kDIO &&
             (sim->CreationCode() == mu2e::ProcessCode::mu2eMuonDecayAtRest ||
              sim->CreationCode() == mu2e::ProcessCode::DIO ||
              sim->CreationCode() == mu2e::ProcessCode::mu2eDIOLeadingLog)) {
            mc_selection = false;
          }
          // only consider IPA DIO tracks above 70 MeV/c if they're from the IPA DIO dataset
          const bool is_ipa_dio = IsIPADIO(sim);
          if(is_ipa_dio && sim->fStartMom.P() > 70. && !fDataset.BeginsWith("ipad"))
            mc_selection = false;
        }
      }
    }
    if(!mc_selection)
      ID += 1 << kMC;

    return ID;
  }

  //_____________________________________________________________________________
  void TConvAnaModule::Debug() {
    int n_high_p_apr_hel(0);
    for(int i = 0; i < fEvtPar.fNAprHelices; i++) {
      const auto helix = fAprHelixBlock->Helix(i);
      InitHelixPar(helix, &fHlxPar);
      if(helix->P() > 70.f)
        ++n_high_p_apr_hel;
    }

    bool print_event = GetDebugBit(0);                     // all events
    print_event |= GetDebugBit(1) && n_high_p_apr_hel > 1; // events with multiple high momentum helices
    if(GetDebugBit(2)) {                                   // high time events
      if(fSimpBlock) {
        for(int isim = 0; isim < fSimpBlock->NParticles(); ++isim) {
          const auto sim = fSimpBlock->Particle(isim);
          print_event |= sim->fStartPos.T() > 1000.;
        }
      }
    }
    if(print_event) { // print event info
      auto event = GetEvent();
      printf(">>> Event %5i/%5i/%6i: weight = %.3e\n", event->fRunNumber, event->fSectionNumber, event->fEventNumber, fEvtPar.fWeight);
      printf(" Helices: APR = %2i CPR = %2i Offline = %2i\n", fEvtPar.fNAprHelices, fEvtPar.fNCprHelices, fEvtPar.fNOfflineHelices);
      printf(" Tracks : APR = %2i CPR = %2i Offline = %2i\n", fEvtPar.fNAprTracks, fEvtPar.fNCprTracks, fEvtPar.fNTracks);
      if(fSimpBlock) {
        printf(" SIM particles:\n");
        fSimpBlock->Print();
      }
      printf(" APR helices:\n");
      fAprHelixBlock->Print();
      printf(" CPR helices:\n");
      fCprDeHelixBlock->Print();
      printf(" Offline De helices:\n");
      fOfflineDeHelixBlock->Print();
      printf(" Offline Ue helices:\n");
      fOfflineUeHelixBlock->Print();
      printf(" APR tracks:\n");
      fAprTrackBlock->Print();
      printf(" CPR tracks:\n");
      fCprTrackBlock->Print();
      printf(" Offline De tracks:\n");
      fOfflineDeTrackBlock->Print();
      printf(" Offline Ue tracks:\n");
      fOfflineUeTrackBlock->Print();
      if(fCRVBlock) {
        printf(" CRV info:\n");
        fCRVBlock->Print();
      }
    } else if(GetDebugBit(5) || (GetDebugBit(6) && !fPrimary)) { // print all sim blocks
      auto event = GetEvent();
      printf(">>> Event %5i/%5i/%6i\n", event->fRunNumber, event->fSectionNumber, event->fEventNumber);
      if(fSimpBlock)
        fSimpBlock->Print();
    }

    if(GetDebugBit(7)) { // bad APR events
      // Report events that triggered but are poorly reconstructed
      for(int itrack = 0; itrack < fAprTrackBlock->NTracks(); ++itrack) {
        const auto track = fAprTrackBlock->Track(itrack);
        const float dp = track->fP - track->fPFront;
        if(dp > 20.f) {
          GetHeaderBlock()->Print();
          track->Print();
        }
      }
    }
  }

  //_____________________________________________________________________________
  int TConvAnaModule::EndJob() {
    TAnaModule::EndJob();
    fWatch->StopTime("Total");
    fWatch->Print(std::cout);

    // Normalize the trigger overlap histograms to N(x && y) / N(x)
    for(int i = 0; i < kNHistSets; i++)
      NormalizeTriggerOverlap(fHist.fEvent[i]);

    // Print interesting cut-flows
    PrintCutFlow(fHist.fTrack[4] , "mumem");

    // Add spectrum and normalization output
    auto fol = (TFolder*)fFolder->FindObject("data");
    if(!fol)
      fol = fFolder->AddFolder("data", "Data");
    if(fNormTree) {
      fNormTree->Fill();
      fol->Add(fNormTree);
    }
    if(fRMCSpectra && fRMCSpectra->hSpectrum_) { // Add the RMC spectrum information to the output
      fol->Add(fRMCSpectra->hSpectrum_);
    }
    if(fStntuple && fStntuple->DioSpectrumHist()) { // Add the DIO spectrum information to the output
      fol->Add(fStntuple->DioSpectrumHist());
    }
    return 0;
  }
} // namespace mumep_ana
