#include "TCanvas.h"
#include "TEnv.h"
#include "TF1.h"
#include "TPad.h"
#include "TSystem.h"

#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/base/TStnDataset.hh"
#include "Stntuple/geom/TCrvNumerology.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

//------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
#include "ConvAna/ana/TAnaModule.hh"
ClassImp(ConvAna::TAnaModule)

    using std::vector;

namespace ConvAna {

  //-----------------------------------------------------------------------------
  TAnaModule::TAnaModule(const char* name, const char* title) : TStnModule(name, title) {

    //-----------------------------------------------------------------------------
    // TStntuple singleton
    //-----------------------------------------------------------------------------
    fStntuple = TStntuple::Instance();

    fRMCSpectra = nullptr;

    fWatch = new mu2e::StopWatch();
    fWatch->Calibrate();
  }

  //-----------------------------------------------------------------------------
  TAnaModule::~TAnaModule() {}

  //-----------------------------------------------------------------------------
  // common initializations
  //-----------------------------------------------------------------------------
  int TAnaModule::BeginJob() { return 0; }

  //-----------------------------------------------------------------------------
  // common initializations
  // assume that the input dataset is initialized.
  //-----------------------------------------------------------------------------
  int TAnaModule::BeginRun() {
    int rn = GetHeaderBlock()->RunNumber();
    TStntuple::Init(rn);
    return 0;
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::BookNormTree() {
    fNormTree = new TTree("Norm", "Normalization info");
    fNormTree->Branch("nseen", &fNormInfo.nseen_);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {

    HBook1F(Hist->fInstLumi[0], "inst_lumi", Form("%s: POT", Folder), 300, 0.0, 1.5e8, Folder);
    HBook1F(Hist->fInstLumi[1], "inst_lumi_1", Form("%s: POT", Folder), 300, 0.0, 1.5e8, Folder);
    HBook1F(Hist->fInstLumiApr, "inst_lumi_apr", Form("%s: POT", Folder), 300, 0.0, 1.5e8, Folder);
    HBook1F(Hist->fInstLumiCpr, "inst_lumi_cpr", Form("%s: POT", Folder), 300, 0.0, 1.5e8, Folder);
    HBook1F(Hist->fInstLumiAprCpr, "inst_lumi_apr_cpr", Form("%s: POT", Folder), 300, 0.0, 1.5e8, Folder);
    HBook1F(Hist->fEventWeight[0], "event_weight", Form("%s: Event weight", Folder), 100, 0., 2., Folder);
    HBook1F(Hist->fEventWeight[1], "event_weight_log", Form("%s: log10(Event weight)", Folder), 100, -20., 1., Folder);
    HBook1F(Hist->fNAprTracks, "nTracksApr", Form("%s: nTracksApr", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNCprTracks, "nTracksCpr", Form("%s: nTracksCpr", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNTracks, "nTracks", Form("%s: nTracks", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNUeTracks, "nUeTracks", Form("%s: nUeTracks", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNDmuTracks, "nDmuTracks", Form("%s: nUmuTracks", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNUmuTracks, "nUmuTracks", Form("%s: nUmuTracks", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNAprHelices, "nHelicesApr", Form("%s: nHelicesApr", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNCprHelices, "nHelicesCpr", Form("%s: nHelicesCpr", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNHelices, "nHelices", Form("%s: nHelices", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNCRVClusters, "nCRVClusters", Form("%s: N(CRV clusters)", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNGoodCRVClusters, "nGoodCRVClusters", Form("%s: N(Good CRV clusters)", Folder), 50, 0.0, 50.0, Folder);
    HBook1F(Hist->fNonCRVVetoID, "nonCRVVetoID", Form("%s: Non-CRV Veto ID", Folder), 30, 0., 30., Folder);
    HBook1F(Hist->fTrackerHits, "tracker_hits", Form("%s: N(tracker) hits", Folder), 100, 0., 12000., Folder);
    HBook1F(Hist->fCaloHits, "calo_hits", Form("%s: N(calo) hits", Folder), 100, 0., 4000., Folder);
    HBook1F(Hist->fTrigBits[0], "trig_bits", Form("%s: Trigger bits", Folder), 300, 100, 400, Folder);
    HBook1F(Hist->fTrigPaths[0], "trig_paths", Form("%s: Trigger paths", Folder), 150, 0, 150, Folder);
    HBook2F(Hist->fTrigOverlap[0], "trig_overlap", Form("%s: Trigger overlap: N(X and Y) / N(X)", Folder), 70, 0, 70, 70, 0, 70, Folder);
    HBook1F(Hist->fTrigBits[1], "trig_bits_ps", Form("%s: Trigger bits", Folder), 300, 100, 400, Folder);
    HBook1F(Hist->fTrigPaths[1], "trig_paths_ps", Form("%s: Trigger paths", Folder), 150, 0, 150, Folder);
    HBook1D(Hist->fTriggered, "triggered", Form("%s: Triggered", Folder), 2, -0.5, 1.5, Folder);
    HBook1F(Hist->fNTriggerable, "ntriggerable", Form("%s: N(triggerable sims)", Folder), 5, -0.5, 4.5, Folder);
    HBook1F(Hist->fPrimaryCode, "primary_code", Form("%s: Primary process code", Folder), 200, -0.5, 199.5, Folder);
    HBook1F(Hist->fPrimaryType, "primary_type", Form("%s: Primary process type", Folder), 50, 0., 50., Folder);
    HBook1F(Hist->fPrimaryGenE, "primary_gene", Form("%s: Primary gen energy", Folder), 500, 50., 150., Folder);
    HBook1F(Hist->fRMCEnergy, "rmc_energy", Form("%s: RMC photon energy", Folder), 440, 0., 110., Folder);

    // Initialize the trigger labels so they're stable
    fTrigBitToBin.clear();
    int trig_bin = 1; // bit -> histogram bin mapping
    for(int bit = 100; bit < 1000; ++bit) {
      TString path = TriggerInfo::BitToName(bit);
      if(path.BeginsWith("Unknown"))
        continue;
      fTrigBitToBin[bit] = trig_bin;
      // Hist->fTrigPaths[0]->SetBinLabel(bin, path.Data());
      // Hist->fTrigPaths[1]->SetBinLabel(bin, path.Data());
      Hist->fTrigPaths[0]->Fill(path.Data(), 0.);
      Hist->fTrigPaths[1]->Fill(path.Data(), 0.);
      Hist->fTrigOverlap[0]->Fill(path.Data(), path.Data(), 0.);
      ++trig_bin;
    }

    // Initialize the primary type axis labels so they're stable
    std::map<TString, bool> types;
    int bin = 1;
    for(int icode = 50; icode < 400; ++icode) {
      TString type = ProcessGroup(icode);
      if(type.BeginsWith("Unknown"))
        continue;
      if(types[type])
        continue;
      types[type] = true;
      // Hist->fPrimaryType->GetXaxis()->SetBinLabel(bin, type.Data());
      Hist->fPrimaryType->Fill(type.Data(), 0.);
      ++bin;
    }
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
    //   char name [200];
    //   char title[200];
    HBook1F(Hist->fP, "p", Form("%s: Momentum", Folder), 1000, 0, 200, Folder);
    HBook1F(Hist->fPdgCode[0], "pdg_code_0", Form("%s: PDG Code[0]", Folder), 200, -100, 100, Folder);
    HBook1F(Hist->fPdgCode[1], "pdg_code_1", Form("%s: PDG Code[1]", Folder), 500, -2500, 2500, Folder);
    HBook1F(Hist->fGenID, "gen_id", Form("%s: Generator ID", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fZ0, "z0", Form("%s: Z0", Folder), 500, 5400, 6400, Folder);
    HBook1F(Hist->fT0, "t0", Form("%s: T0", Folder), 200, 0, 2000, Folder);
    HBook1F(Hist->fR0, "r", Form("%s: R0", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fCosTh, "cos_th", Form("%s: Cos(Theta)", Folder), 200, -1., 1., Folder);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
    //  char name [200];
    //  char title[200];

    HBook1F(Hist->fPdgCode, "pdg", Form("%s: PDG code", Folder), 600, -300, 300, Folder);
    HBook1F(Hist->fCreationCode, "code", Form("%s: Creation code", Folder), 200, 0, 200, Folder);
    HBook1F(Hist->fNStrawHits, "nsth", Form("%s: n straw hits", Folder), 200, 0, 200, Folder);
    HBook1F(Hist->fMomTargetEnd, "ptarg", Form("%s: mom after Stopping Target", Folder), 400, 70, 110, Folder);
    HBook1F(Hist->fMomTrackerFront, "pfront", Form("%s: mom at the Tracker Front", Folder), 400, 70, 110, Folder);
    HBook1F(Hist->fMomStart, "pstart", Form("%s: mom at creation", Folder), 400, 0, 150, Folder);
    HBook1F(Hist->fTimeStart, "tstart", Form("%s: time at creation", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fTimeEnd, "tend", Form("%s: time at end", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fStartZ, "startz", Form("%s: z at creation", Folder), 500, 0, 20000, Folder);
    HBook1F(Hist->fStartR, "startr", Form("%s: r at creation", Folder), 200, 0, 800, Folder);
    HBook2F(Hist->fStartXY, "startxy", Form("%s: (x,y) at creation;x;y", Folder), 200, -800, 800, 200, -800, 800, Folder);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {

    HBook1F(Hist->fP[0], "p", Form("%s: track momentum", Folder), 600, -300.0, 300.0, Folder);
    HBook1F(Hist->fP[1], "p_2", Form("%s: track momentum", Folder), 600, 80., 110., Folder);
    HBook1F(Hist->fObs, "obs", Form("%s: Track momentum", Folder), 300, 80., 110., Folder); // fit histogram
    HBook1F(Hist->fPt, "pt", Form("%s: track transverse momentum", Folder), 300, 0., 300.0, Folder);
    HBook1F(Hist->fPCorr, "p_corr", Form("%s: corrected track momentum", Folder), 600, 80., 110., Folder);
    HBook1F(Hist->fPTrkFront, "p_trk_front", Form("%s: track momentum at TT_Front", Folder), 300, 0., 150., Folder);
    HBook1F(Hist->fPCenter[0], "pCenter", Form("%s: track momentum at tracker center", Folder), 600, -300.0, 300.0, Folder);
    HBook1F(Hist->fPCenter[1], "pCenter_2", Form("%s: track momentum at tracker center", Folder), 600, 80., 110., Folder);
    HBook1F(Hist->fPExit, "pExit", Form("%s: track momentum at tracker exit", Folder), 600, -300.0, 300.0, Folder);
    HBook1F(Hist->fPST[0], "pST", Form("%s: track momentum at ST exit", Folder), 600, -300.0, 300.0, Folder);
    HBook1F(Hist->fPST[1], "pST_2", Form("%s: track momentum at ST exit", Folder), 600, 80.0, 110.0, Folder);
    HBook1F(Hist->fPSTDiff, "pST_diff", Form("%s: track p(ST) - p(Front)", Folder), 400, -1.0, 9.0, Folder);
    HBook1F(Hist->fPSTApproxDiff, "p_approx_ST_diff", Form("%s: track p(ST) - p(Front)", Folder), 400, -1.0, 9.0, Folder);
    HBook1F(Hist->fPExitDiff, "pExit_diff", Form("%s: track p(Front) - p(Exit)", Folder), 400, -1.0, 4.0, Folder);
    HBook1F(Hist->fT0, "t0", Form("%s: track t_{0}", Folder), 400, 0., 2000., Folder);
    HBook1F(Hist->fT0Err, "t0err", Form("%s: track t_{0} uncertainty", Folder), 200, 0., 4., Folder);
    HBook1F(Hist->fD0, "d0", Form("%s: track d0", Folder), 400, -400.0, 400.0, Folder);
    HBook1F(Hist->fDP, "dP", Form("%s: track p_reco - p_mc", Folder), 400, -20.0, 20.0, Folder);
    HBook2F(Hist->fDPvsP, "dPvsP", Form("%s: track p_reco - p_mc", Folder), 20, 80, 120, 100, -5.0, 5.0, Folder);
    HBook2F(Hist->fDPvsNH, "dPvsNH", Form("%s: track p_reco - p_mc", Folder), 20, 10, 70, 100, -5.0, 5.0, Folder);
    HBook1F(Hist->fDPCorr, "dPCorr", Form("%s: track p_reco - p_mc", Folder), 400, -20.0, 20.0, Folder);
    HBook1F(Hist->fChi2NDof, "chi2NDof", Form("%s: track chi2/ndof", Folder), 200, 0.0, 10.0, Folder);
    HBook1F(Hist->fFitCons[0], "fitCons", Form("%s: track p(chi2,ndof)", Folder), 200, 0.0, 1.0, Folder);
    HBook1F(Hist->fFitCons[1], "fitCons_log", Form("%s: track log10(p(chi2,ndof))", Folder), 200, -6., 0.0, Folder);
    HBook1F(Hist->fFitMomErr, "fitMomErr", Form("%s: track momentum uncertainty", Folder), 100, 0.0, 1.0, Folder);
    HBook1F(Hist->fTanDip, "tanDip", Form("%s: track tanDip", Folder), 200, 0.0, 2.0, Folder);
    HBook1F(Hist->fCosTheta, "cosTheta", Form("%s: track cos(#theta)", Folder), 200, -1.0, 1.0, Folder);
    HBook1F(Hist->fRadius, "radius", Form("%s: track radius", Folder), 1000, 0.0, 1000, Folder);
    HBook1F(Hist->fRMax, "rMax", Form("%s: track rMax", Folder), 2000, 0.0, 2000, Folder);
    HBook1F(Hist->fNActive, "nActive", Form("%s: nHits used in fit", Folder), 150, 0.0, 150.0, Folder);
    HBook1F(Hist->fNActiveFrac, "nActiveFrac", Form("%s: nHits used in fit/nHits", Folder), 100, 0.0, 1.0, Folder);
    HBook1F(Hist->fTrkQual[0], "trkQual", Form("%s: track MVA score", Folder), 200, -1.0, 1.0, Folder);
    HBook1F(Hist->fTrkQual[1], "trkQual_1", Form("%s: track MVA score", Folder), 200, -1.0, 1.0, Folder);
    HBook1F(Hist->fPID[0], "pid", Form("%s: PID MVA score", Folder), 200, -1.0, 1.0, Folder);
    HBook1F(Hist->fPID[1], "trkpid", Form("%s: TrkPID MVA score", Folder), 200, -1.0, 1.0, Folder);
    HBook1F(Hist->fCosmicID, "cosmic_id", Form("%s: Cosmic MVA score", Folder), 200, -1.0, 1.0, Folder);
    HBook1F(Hist->fClusterE, "clusterE", Form("%s: track's cluster energy", Folder), 600, 0., 300., Folder);
    HBook1D(Hist->fClusterDisk, "clusterDisk", Form("%s: track's cluster energy", Folder), 3, -1., 2., Folder);
    HBook1F(Hist->fDt, "dt", Form("%s: track - cluster time", Folder), 200, -10., 10., Folder);
    HBook1F(Hist->fEp, "ep", Form("%s: cluster E / track P", Folder), 200, 0., 2., Folder);
    HBook1F(Hist->fTZSlope, "tzslope", Form("%s: TZ slope", Folder), 200, -0.1, 0.1, Folder);
    HBook1F(Hist->fTZSlopeSig, "tzslopesig", Form("%s: TZ slope significance", Folder), 200, -10., 10., Folder);
    HBook1F(Hist->fTZSlopeRatio, "tzsloperatio", Form("%s: TZ slope / expected slope", Folder), 200, -10., 10., Folder);
    HBook1F(Hist->fBestAlg, "bestAlg", Form("%s: Best fit algorithm", Folder), 10, 0., 10., Folder);
    HBook1F(Hist->fAlgMask, "algMask", Form("%s: Algorithm mask", Folder), 10, 0., 10., Folder);
    HBook1F(Hist->fSTBoundary, "st_boundary", Form("%s: Stopping target boundary", Folder), 2, 0., 2., Folder);
    HBook1F(Hist->fSTInters, "st_inters", Form("%s: Stopping target intersections", Folder), 10, 0., 10., Folder);
    HBook1F(Hist->fIPAInters, "ipa_inters", Form("%s: IPA intersections", Folder), 10, 0., 10., Folder);
    HBook1F(Hist->fOPAInters, "opa_inters", Form("%s: OPA intersections", Folder), 10, 0., 10., Folder);
    HBook1F(Hist->fTrackID, "track_id", Form("%s: Track ID bits", Folder), 33, 0., 33., Folder);
    HBook1F(Hist->fExlTrackID, "track_exl_id", Form("%s: Track ID bits for exclusive rejection", Folder), 32, 0., 32., Folder);

    // Initialize bin labels for the Track ID histograms
    Hist->fTrackID->GetXaxis()->SetBinLabel(1, "Passed");
    Hist->fExlTrackID->GetXaxis()->SetBinLabel(1, "Passed");
    int current_bin = 2;
    for(int bit = 0; bit <= Hist->fTrackID->GetNbinsX(); ++bit) {
      if(!kTrackIDNames.count(bit))
        continue;
      TString name(kTrackIDNames[bit]);
      Hist->fTrackID->GetXaxis()->SetBinLabel(current_bin, name.Data());
      Hist->fExlTrackID->GetXaxis()->SetBinLabel(current_bin, name.Data());
      ++current_bin;
    }

    HBook1F(Hist->fCRVDeltaT, "crv_deltat", Form("%s: Track t_{0} - CRV t_{0}", Folder), 200, -300., 300., Folder);
    HBook1F(Hist->fCRVDeltaTCRV, "crv_deltat_crv", Form("%s: Track t_{0} - CRV t_{0}", Folder), 200, -300., 300., Folder);
    HBook1F(Hist->fCRVDeltaTST, "crv_deltat_st", Form("%s: Track t_{0} - CRV through ST t_{0}", Folder), 200, -300., 300., Folder);
    HBook1F(Hist->fCRVDeltaTCalo[0], "crv_deltat_calo", Form("%s: Track t_{0} - CRV through Calo t_{0}", Folder), 200, -300., 300., Folder);
    HBook1F(Hist->fCRVDeltaTCalo[1], "crv_deltat_calo_mu", Form("%s: Track t_{0} - CRV through Calo t_{0}", Folder), 200, -300., 300., Folder);
    HBook1F(Hist->fCRVDeltaTExtrap, "crv_deltat_extrap", Form("%s: Track t_{0} - CRV through Extrapolated t_{0}", Folder), 200, -300., 300., Folder);
    HBook1F(Hist->fCRVMinDeltaT, "crv_min_deltat", Form("%s: Min #Deltat_{0} for ST and Calo paths", Folder), 200, -300., 300., Folder);
    HBook1F(Hist->fCRVExtrapZ, "crv_extrap_z", Form("%s: Extrapolation Z", Folder), 200, 0., 15000., Folder);
    HBook2F(Hist->fCRVXZ, "crv_x_vs_z", Form("%s: CRV X vs Z", Folder), 250, -5000, 25000, 200, -10000, 10000, Folder);
    HBook2F(Hist->fCRVYZ, "crv_y_vs_z", Form("%s: CRV Y vs Z", Folder), 250, -5000, 25000, 200, 0, 8000, Folder);
    HBook2F(Hist->fCRVdTZ, "crv_dt_vs_z", Form("%s: #Deltat vs CRV Z", Folder), 250, -5000, 25000, 200, -200, 200, Folder);
    HBook2F(Hist->fCRVdTZCRV, "crv_dtcrv_vs_z", Form("%s: #Deltat(CRV) vs CRV Z", Folder), 250, -5000, 25000, 200, -200, 200, Folder);
    HBook2F(Hist->fCRVdTZExtrap, "crv_dtextrap_vs_z", Form("%s: #Deltat(Extrapolation) vs CRV Z", Folder), 250, -5000, 25000, 200, -200, 200, Folder);

    HBook1F(Hist->fUpstreamDt, "us_dt", Form("%s: Upstream track #Deltat_{0}", Folder), 150, -50., 250., Folder);
    HBook1F(Hist->fUpstreamDp, "us_dp", Form("%s: Upstream track #Deltap", Folder), 100, -10., 5., Folder);

    HBook1F(Hist->fMCPFront, "MC_PFront", Form("%s: MC track P(tracker front)", Folder), 600, 0., 300., Folder);
    HBook1F(Hist->fMCPStOut, "MC_PSTOut", Form("%s: MC track P(ST exit)", Folder), 600, 0., 300., Folder);
    HBook1F(Hist->fMCPStDiff, "MC_PSTDiff", Form("%s: MC track P(ST exit) - P(tracker front)", Folder), 400, -20.0, 20.0, Folder);
    HBook1F(Hist->fMCPStDiffDiff, "MC_PSTDiffDiff", Form("%s: Reco - MC track (P(ST exit) - P(Tracker front))", Folder), 100, -5., 5., Folder);
    HBook1F(Hist->fMCApproxPStDiffDiff, "MC_ApproxPSTDiffDiff", Form("%s: Reco - MC track (P(ST exit) - P(Tracker front))", Folder), 100, -5., 5., Folder);
    HBook1F(Hist->fMCGenE, "MC_GenE", Form("%s: MC generated energy", Folder), 400, 0., 200., Folder);
    HBook1F(Hist->fMCPGenEDiff, "MC_PGenEDiff", Form("%s: P(reco front) - MC E(gen)", Folder), 500, -10., 5., Folder);
    HBook1F(Hist->fMCPSig, "MC_PSig", Form("%s: MC P(front) error / uncertainty", Folder), 400, -20., 20., Folder);
    HBook1F(Hist->fMCPdg[0], "MC_PDG_0", Form("%s: MC Particle PDG code", Folder), 40, -20., 20., Folder);
    HBook1F(Hist->fMCPdg[1], "MC_PDG_1", Form("%s: MC Particle |PDG code|", Folder), 220, 0., 2200., Folder);
    HBook1F(Hist->fMCStrawHits, "MC_strawhits", Form("%s: MC Particle N(straw hits)", Folder), 200, 0., 200., Folder);
    HBook1F(Hist->fMCGoodHits, "MC_goodhits", Form("%s: MC Particle N(good hits)", Folder), 100, 0., 100., Folder);
    HBook1F(Hist->fMCTrajectory, "MC_trajectory", Form("%s: MC track p_{z} trajectory", Folder), 3, -1.5, 1.5, Folder);
    HBook1F(Hist->fMCSimProc, "MC_simProc", Form("%s: MC Sim process code", Folder), 200, -0.5, 199.5, Folder);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::BookHelixHistograms(HelixHist_t* Hist, const char* Folder) {

    HBook1F(Hist->fNHits, "nhits", Form("%s: # of straw hits", Folder), 150, 0, 150, Folder);
    HBook1F(Hist->fHelicity, "hel", Form("%s: Helicity", Folder), 10, -5, 5, Folder);
    HBook1F(Hist->fClusterTime, "clusterTime", Form("%s: cluster time; t_{cluster}[ns]", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fClusterEnergy, "clusterE", Form("%s: cluster energy; E [MeV]      ", Folder), 400, 0, 200, Folder);
    HBook1F(Hist->fRadius, "radius", Form("%s: helix radius; r [mm]", Folder), 500, 0, 500, Folder);
    HBook1F(Hist->fRMax, "rmax", Form("%s: helix R(max); R(max) [mm]", Folder), 200, 0, 800, Folder);
    HBook1F(Hist->fMom, "p", Form("%s: momentum; p [MeV/c]", Folder), 400, 0, 200, Folder);
    HBook1F(Hist->fPt, "pT", Form("%s: pT; pT [MeV/c]", Folder), 600, 0, 150, Folder);
    HBook1F(Hist->fGenMom, "simMom", Form("%s: Sim particle P", Folder), 400, 0, 200, Folder);
    HBook1F(Hist->fGenID, "simPDG", Form("%s: Sim particle PDG type", Folder), 60, -30, 30, Folder);
    HBook2F(Hist->fGenRZ, "simRZ", Form("%s: Sim particle origin; Z (mm); R (mm)", Folder), 200, 5000., 15000., 200, 0, 1000., Folder);
    HBook1F(Hist->fDp, "dp", Form("%s: Reco p - MC p", Folder), 200, -10, 10, Folder);
    HBook1F(Hist->fDpT, "dpT", Form("%s: Reco - Sim particle P_{T}", Folder), 400, -50, 50, Folder);
    HBook1F(Hist->fLambda, "lambda", Form("%s: lambda; #lambda", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fTanDip, "tanDip", Form("%s: tanDip", Folder), 200, -2.0, 2.0, Folder);
    HBook1F(Hist->fT0, "t0", Form("%s: t0; t0[ns]", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fT0Err, "t0err", Form("%s: t0err; t0err [ns]", Folder), 100, 0, 10, Folder);
    HBook1F(Hist->fD0, "d0", Form("%s: D0; d0 [mm]", Folder), 1600, -400, 400, Folder);
    HBook1F(Hist->fAlgMask, "algmask", Form("%s: Algorithm Mask", Folder), 10, 0, 10, Folder);
    HBook1F(Hist->fBestAlg, "bestalg", Form("%s: Best Algorithm", Folder), 10, 0, 10, Folder);
    HBook1F(Hist->fChi2XY, "chi2xy", Form("%s: Chi2(XY)/DOF", Folder), 100, 0, 10, Folder);
    HBook1F(Hist->fChi2ZPhi, "chi2zphi", Form("%s: Chi2(ZPhi)/DOF", Folder), 100, 0, 10, Folder);
    HBook1F(Hist->fTZSlope, "tzslope", Form("%s: dz/dt", Folder), 100, -0.1, 0.1, Folder);
    HBook1F(Hist->fTZSlopeErr, "tzslopeerr", Form("%s: dz/dt error", Folder), 100, 0., 0.1, Folder);
    HBook1F(Hist->fTZSlopeSig, "tzslopesig", Form("%s: dz/dt significance", Folder), 200, -10, 10, Folder);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::BookClusterHistograms(ClusterHist_t* Hist, const char* Folder) {
    HBook1D(Hist->fDiskID, "disk_id", Form("%s: Disk ID", Folder), 2, 0, 2, Folder);
    HBook1F(Hist->fEnergy, "energy", Form("%s: Cluster Energy", Folder), 500, 0, 250, Folder);
    HBook1F(Hist->fT0, "t0", Form("%s: cluster T0", Folder), 200, 0, 2000, Folder);
    HBook1F(Hist->fRow, "row", Form("%s: cluster Row", Folder), 200, 0, 200, Folder);
    HBook1F(Hist->fCol, "col", Form("%s: cluster column", Folder), 200, 0, 200, Folder);
    HBook1F(Hist->fX, "x", Form("%s: cluster X", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fY, "y", Form("%s: cluster Y", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fZ, "z", Form("%s: cluster Z", Folder), 200, -10, 10, Folder);
    HBook1F(Hist->fR, "r", Form("%s: cluster Radius", Folder), 100, 300, 800, Folder);
    HBook1F(Hist->fYMean, "ymean", Form("%s: cluster YMean", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fZMean, "zmean", Form("%s: cluster ZMean", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fSigY, "sigy", Form("%s: cluster SigY", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fSigZ, "sigz", Form("%s: cluster SigZ", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fSigR, "sigr", Form("%s: cluster SigR", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fNCr0, "ncr0", Form("%s: cluster NCR[0]", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fNCr1, "ncr1", Form("%s: cluster NCR[1]", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fFrE1, "fre1", Form("%s: E1/Etot", Folder), 220, 0, 1.1, Folder);
    HBook1F(Hist->fFrE2, "fre2", Form("%s: (E1+E2)/Etot", Folder), 220, 0, 1.1, Folder);
    HBook1F(Hist->fSigE1, "sige1", Form("%s: SigmaE/Etot", Folder), 200, 0, 10, Folder);
    HBook1F(Hist->fSigE2, "sige2", Form("%s: SigmaE/Emean", Folder), 200, 0, 10, Folder);
    HBook1F(Hist->fTimeRMS, "time_rms", Form("%s: T(RMS)", Folder), 200, 0, 10, Folder);
    HBook1F(Hist->fMaxR, "maxr", Form("%s: Max R from main", Folder), 200, 0, 400, Folder);
    HBook1F(Hist->fE9OverE, "e9_over_e", Form("%s: E(3x3)/E", Folder), 220, 0, 1.1, Folder);
    HBook1F(Hist->fE25OverE, "e25_over_e", Form("%s: E(5x5)/E", Folder), 220, 0, 1.1, Folder);
    HBook1F(Hist->fRingEOverE, "ring_e_over_e", Form("%s: E(ring)/E", Folder), 220, 0, 1.1, Folder);
    HBook1F(Hist->fRingEOverE1, "ring_e_over_e1", Form("%s: E(ring)/E1", Folder), 200, 0, 3.0, Folder);
    HBook1F(Hist->fOutRingE, "out_ring_e", Form("%s: E(out ring)", Folder), 200, 0, 200, Folder);
    HBook1F(Hist->fOutRingEOverE, "out_ring_e_over_e", Form("%s: E(out ring)/E", Folder), 200, 0, 1.1, Folder);
    HBook1F(Hist->fNCoreCrystals, "n_core_crystals", Form("%s: N(core crystals)", Folder), 20, 0., 20., Folder);
    HBook1F(Hist->fCoreEnergy, "core_energy", Form("%s: Core Energy", Folder), 200, 0., 150., Folder);
    HBook1F(Hist->fCoreEnergyFrac, "core_energy_frac", Form("%s: Core Energy/Energy", Folder), 200, 0., 1.1, Folder);
    HBook1F(Hist->fMCSimEDep, "MC_sim_edep", Form("%s: MC Sim E(dep)", Folder), 250, 0, 250, Folder);
    HBook1F(Hist->fMCSimMomIn, "MC_sim_mom_in", Form("%s: MC Sim mom(in)", Folder), 250, 0, 250, Folder);
    HBook1F(Hist->fMCSimEStart, "MC_sim_estart", Form("%s: MC Sim E(start)", Folder), 250, 0, 250, Folder);
    HBook1F(Hist->fMCSimPdg, "MC_sim_pdg", Form("%s: MC Sim PDG ID", Folder), 200, -100, 100, Folder);
    HBook1F(Hist->fMCEDep, "MC_edep", Form("%s: MC E(dep)", Folder), 250, 0, 250, Folder);
    HBook1F(Hist->fMCTime, "MC_time", Form("%s: MC time", Folder), 200, 0, 2000, Folder);
    HBook1F(Hist->fMC_dE, "MC_dE", Form("%s: MC #DeltaE", Folder), 200, -20, 20, Folder);
    HBook1F(Hist->fMC_dt, "MC_dt", Form("%s: MC #Deltat", Folder), 200, -5, 5, Folder);
    HBook1F(Hist->fMC_dGenE, "MC_d_gen_e", Form("%s: Reco E - Gen E", Folder), 200, -40, 20, Folder);
    HBook1F(Hist->fMCGenE, "MC_gen_e", Form("%s: Gen E", Folder), 200, 0, 200, Folder);
    HBook2F(Hist->fMCEvsE, "MC_e_vs_e", Form("%s: MC E vs. E;MC E;Reco E", Folder), 150, 0, 150, 150, 0, 150, Folder);
    HBook2F(Hist->fMCGenEvsE, "MC_gen_e_vs_e", Form("%s: Gen E vs. E;Gen E;Reco E", Folder), 150, 0, 150, 150, 0, 150, Folder);
    HBook1F(Hist->fMCSimPdgName, "MC_sim_pdg_name", Form("%s: MC Sim Name", Folder), 20, 0, 20, Folder);

    // Initialize the bin labels
    vector<int> pdgs = {11, 13, 22, 2212, 2112, 211};
    for(auto pdg : pdgs)
      Hist->fMCSimPdgName->Fill(NameFromPDG(pdg).Data(), 0.);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::BookCRVClusterHistograms(CRVClusterHist_t* Hist, const char* Folder) {
    HBook1F(Hist->fSector, "sector", Form("%s: CRV sector", Folder), 30, 0, 30, Folder);
    HBook1F(Hist->fFirstBar, "fbar", Form("%s: first pulse bar#", Folder), 600, 0, 6000, Folder);
    HBook1F(Hist->fNPulses, "npulses", Form("%s: N(pulses)", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fNPe, "npe", Form("%s: N(PE)", Folder), 500, 0, 5000, Folder);
    HBook1F(Hist->fNPePP, "npepp", Form("%s: N(PE) per pulse", Folder), 500, 0, 500, Folder);
    HBook1F(Hist->fStartTime, "tstart", Form("%s: start time, ns", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fEndTime, "tend", Form("%s: end time, ns", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fWidth, "wwidth", Form("%s: width, ns", Folder), 200, 0, 200, Folder);
    HBook2F(Hist->fXVsZ, "x_vs_z", Form("%s: X vs Z", Folder), 250, -5000, 25000, 200, -10000, 10000, Folder);
    HBook2F(Hist->fYVsZ, "y_vs_z", Form("%s: Y vs Z", Folder), 250, -5000, 25000, 200, 0, 8000, Folder);
    HBook1F(Hist->fCorrTime, "correctedtime", Form("%s: corrected time", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fCorrTimeProp, "time_prop", Form("%s: time at CRV", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fCorrTimeToF, "tof", Form("%s: time of flight from CRV", Folder), 100, 0, 200, Folder);
    HBook1F(Hist->fApproxTimeST, "apprx_t_st", Form("%s: time at ST center", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fApproxTimeCalo, "apprx_t_calo", Form("%s: time at Calo center", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fApproxTimeExtrap, "apprx_t_extrap", Form("%s: time at Extrapolation", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fApproxTimeSTToFront, "apprx_t_st_front", Form("%s: time at Trk front from ST center", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fApproxTimeCaloToFront, "apprx_t_calo_front", Form("%s: time at Trk front from Calo center", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fApproxTimeExtrapToFront, "apprx_t_extrap_front", Form("%s: time at Trk front from extrapolation", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fBarsOneEnd, "barsoneend", Form("%s: one ended bars", Folder), 20, 0, 20, Folder);
    HBook1F(Hist->fCrvPropdT, "crvpropdt  ", Form("%s: dT between CorrPropTime and StartTime", Folder), 200, -50, 50, Folder);
    HBook1F(Hist->fNSectors, "nsectors", Form("%s: Number of sectors in a CRV Cluster", Folder), 20, 0, 20, Folder);
    HBook1F(Hist->fNDiffLSectors, "ndifflsectors", Form("%s: Number of sectors in a CRV Cluster with different lengths", Folder), 20, 0, 20, Folder);
    HBook1F(Hist->fBarsTwoEnd, "barstwoend", Form("%s: two ended bars", Folder), 20, 0, 20, Folder);
    HBook1F(Hist->fStubSlope, "stub_slope", Form("%s: local stub slope", Folder), 200, -5, 5, Folder);
    HBook1F(Hist->fStubSlopeChi2, "stub_slope_chi2", Form("%s: stub slope chi2", Folder), 200, 0, 20, Folder);
    HBook1F(Hist->fStubSlopeDelta, "stub_slope_delta", Form("%s: delta local stub slope - MC stub slope; slope_local - slope_MC", Folder), 200, -5, 5, Folder);
    HBook1F(Hist->fStubQN, "stub_qn", Form("%s: stub qn: # of points in localXY", Folder), 10, 0, 10, Folder);
    HBook1F(Hist->fStubSlopeMCProduct, "stub_slope_prod", Form("%s: stub slope product: reco slope * MC slope", Folder), 200, -10, 10, Folder);
    HBook1F(Hist->fMCTime, "MC_time", Form("%s: MC cluster time;time (ns)", Folder), 200, 0., 2000., Folder);
    HBook1F(Hist->fMCdT[0], "MC_dt", Form("%s: Reco - MC cluster time;#Deltat (ns)", Folder), 200, -200., 200., Folder);
    HBook1F(Hist->fMCdT[1], "MC_dtCorr", Form("%s: Corrected Reco - MC cluster time;#Deltat (ns)", Folder), 200, -200., 200., Folder);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::FillEventHistograms(EventHist_t* Hist, EventPar_t* EvtPar, float Weight) {
    if(!Hist || !EvtPar)
      return;
    if(fDebugLevel > 0)
      fWatch->SetTime("TAna::FillEventHistograms");

    Hist->fInstLumi[0]->Fill(EvtPar->fInstLum, Weight);
    Hist->fInstLumi[1]->Fill(EvtPar->fInstLum, Weight);
    if(EvtPar->fPassedAprPath) {
      Hist->fInstLumiApr->Fill(EvtPar->fInstLum, Weight);
    }
    if(EvtPar->fPassedCprPath) {
      Hist->fInstLumiCpr->Fill(EvtPar->fInstLum, Weight);
    }
    if(EvtPar->fPassedAprPath || EvtPar->fPassedCprPath) {
      Hist->fInstLumiAprCpr->Fill(EvtPar->fInstLum, Weight);
    }
    Hist->fEventWeight[0]->Fill(Weight);
    Hist->fEventWeight[1]->Fill((Weight > 0.f) ? std::log10(Weight) : -1.e3);
    Hist->fNAprTracks->Fill(EvtPar->fNAprTracks, Weight);
    Hist->fNCprTracks->Fill(EvtPar->fNCprTracks, Weight);
    Hist->fNTracks->Fill(EvtPar->fNTracks, Weight);
    Hist->fNUeTracks->Fill(EvtPar->fNUeTracks, Weight);
    Hist->fNDmuTracks->Fill(EvtPar->fNDmuTracks, Weight);
    Hist->fNUmuTracks->Fill(EvtPar->fNUmuTracks, Weight);
    Hist->fNAprHelices->Fill(EvtPar->fNAprHelices, Weight);
    Hist->fNCprHelices->Fill(EvtPar->fNCprHelices, Weight);
    Hist->fNHelices->Fill(EvtPar->fNOfflineHelices, Weight);
    Hist->fNCRVClusters->Fill(EvtPar->fNCRVClusters, Weight);
    Hist->fNGoodCRVClusters->Fill(EvtPar->fNGoodCRVClusters, Weight);
    for(int bit = 0; bit < 30; ++bit) {
      if((EvtPar->fNonCRVVetoID & (1 << bit)) != 0)
        Hist->fNonCRVVetoID->Fill(bit, Weight);
    }
    Hist->fTrackerHits->Fill(EvtPar->fTrackerHits, Weight);
    Hist->fCaloHits->Fill(EvtPar->fCaloHits, Weight);

    if(fMakeTrigHists && fTriggerBlock) {
      double ps_wt = 1.;
      bool triggered = false;
      for(int bit : fTriggersPassed) {
        Hist->fTrigBits[0]->Fill(bit, Weight);
        Hist->fTrigPaths[0]->Fill(TriggerInfo::BitToName(bit).Data(), Weight);
        const double ps = TriggerInfo::BitToPrescale(bit);
        if(ps > 0.) {
          triggered = true;
          ps_wt *= (1. - 1. / ps);
          Hist->fTrigBits[1]->Fill(bit, Weight * (1. / ps));
          Hist->fTrigPaths[1]->Fill(TriggerInfo::BitToName(bit).Data(), Weight * (1. / ps));
        }
        // Check overlaps
        if(fMakeTrigHists > 1) {
          for(int bit_j : fTriggersPassed) {
            Hist->fTrigOverlap[0]->Fill(TriggerInfo::BitToName(bit).Data(), TriggerInfo::BitToName(bit_j).Data(), Weight);
          }
        }
      }
      if(triggered)
        Hist->fTriggered->Fill(1., Weight * (1. - ps_wt));
      if(ps_wt > 0.)
        Hist->fTriggered->Fill(0., Weight * ps_wt); // even if it passed, p(ps_wt) fraction of the time it doesn't
    }
    Hist->fNTriggerable->Fill(EvtPar->fNTriggerable, Weight);

    TString process = "Unknown";
    const int proc_code = (fPrimary) ? fPrimary->CreationCode() : -1;
    switch(proc_code) {
    case mu2e::ProcessCode::mu2eFlateMinus:
    case mu2e::ProcessCode::mu2eDIOLeadingLog:
    case mu2e::ProcessCode::mu2eMuonDecayAtRest: // FIXME: Ensure this isn't a frequent code
      process = "DIO";
      break;
    case mu2e::ProcessCode::mu2eExternalRMC:
    case mu2e::ProcessCode::mu2eInternalRMC:
    case mu2e::ProcessCode::mu2eGammaConversion:
    case mu2e::ProcessCode::mu2eFlatPhoton:
      process = "RMC";
      break;
    case mu2e::ProcessCode::mu2eExternalRPC:
    case mu2e::ProcessCode::mu2eInternalRPC:
      process = "RPC";
      break;
    case mu2e::ProcessCode::mu2eAntiproton:
      process = "Antiproton";
      break;
    default:
      process = "Unknown";
    }
    TString proc_type = (fPrimary) ? ProcessGroup(fPrimary->CreationCode()).Data() : "Unknown";
    if(proc_type == "DIO" && IsIPADIO(fPrimary))
      proc_type = "IPA DIO";
    Hist->fPrimaryCode->Fill((fPrimary) ? fPrimary->CreationCode() : 0, Weight);
    Hist->fPrimaryType->Fill(proc_type.Data(), Weight);
    Hist->fPrimaryGenE->Fill((fPrimary) ? fPrimary->fStartMom.E() : 0., Weight);

    Hist->fRMCEnergy->Fill(EvtPar->fRMCEnergy, Weight);
    if(fDebugLevel > 0)
      fWatch->StopTime("TAna::FillEventHistograms");
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::FillGenpHistograms(GenpHist_t* Hist, TGenParticle* Genp, float Weight) {
    if(!Hist || !Genp)
      return;
    int gen_id;
    float p, cos_th, z0, t0, r0, x0, y0;

    TLorentzVector mom, v;

    Genp->Momentum(mom);
    //  Genp->ProductionVertex(v);

    p = mom.P();
    cos_th = mom.CosTheta();

    x0 = Genp->Vx() + 3904.;
    y0 = Genp->Vy();

    z0 = Genp->Vz();
    t0 = Genp->T();
    r0 = sqrt(x0 * x0 + y0 * y0);
    gen_id = Genp->GetStatusCode();

    Hist->fPdgCode[0]->Fill(Genp->GetPdgCode(), Weight);
    Hist->fPdgCode[1]->Fill(Genp->GetPdgCode(), Weight);
    Hist->fGenID->Fill(gen_id, Weight);
    Hist->fZ0->Fill(z0, Weight);
    Hist->fT0->Fill(t0, Weight);
    Hist->fR0->Fill(r0, Weight);
    Hist->fP->Fill(p, Weight);
    Hist->fCosTh->Fill(cos_th, Weight);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp, float Weight) {
    if(!Hist || !Simp)
      return;

    Hist->fPdgCode->Fill(Simp->fPdgCode, Weight);
    Hist->fCreationCode->Fill(Simp->fCreationCode, Weight);
    Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd, Weight);
    Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront, Weight);
    Hist->fNStrawHits->Fill(Simp->fNStrawHits, Weight);
    Hist->fMomStart->Fill(Simp->fStartMom.P(), Weight);
    Hist->fTimeStart->Fill(Simp->fStartPos.T(), Weight);
    Hist->fTimeEnd->Fill(Simp->fEndPos.T(), Weight);

    const float x(Simp->fStartPos.X() + 3904.), y(Simp->fStartPos.Y()), z(Simp->fStartPos.Z());
    const float r(std::sqrt(x * x + y * y));
    Hist->fStartZ->Fill(z, Weight);
    Hist->fStartR->Fill(r, Weight);
    Hist->fStartXY->Fill(x, y, Weight);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::FillTrackHistograms(TrackHist_t* Hist, TrackPar_t* TrkPar, float Weight) {

    if(!Hist) {
      throw std::runtime_error(Form("TAnaModule::%s: Uninitialized Histogram set\n", __func__));
    }
    if(!TrkPar) {
      throw std::runtime_error(Form("TAnaModule::%s: Uninitialized TrkPar\n", __func__));
    }
    if(!TrkPar->fTrack) {
      throw std::runtime_error(Form("TAnaModule::%s: Uninitialized Track\n", __func__));
    }
    Hist->fP[0]->Fill((TrkPar->fTrack->fP) * (TrkPar->fTrack->fCharge), Weight);
    Hist->fP[1]->Fill(TrkPar->fTrack->fP, Weight);
    Hist->fObs->Fill(TrkPar->fObs, Weight);
    Hist->fPCorr->Fill(TrkPar->fTrack->fP + TrkPar->fApproxDpST, Weight);
    Hist->fPt->Fill(TrkPar->fTrack->fPt, Weight);
    Hist->fPCenter[0]->Fill((TrkPar->fTrack->fPTrackerMiddle) * (TrkPar->fTrack->fCharge), Weight);
    Hist->fPCenter[1]->Fill(TrkPar->fTrack->fPTrackerMiddle, Weight);
    Hist->fPTrkFront->Fill(TrkPar->fTrack->fPTrackerEntrance, Weight);
    Hist->fPExit->Fill((TrkPar->fTrack->fPTrackerExit) * (TrkPar->fTrack->fCharge), Weight);
    Hist->fPST[0]->Fill((TrkPar->fTrack->fPSTBack) * (TrkPar->fTrack->fCharge), Weight);
    Hist->fPST[1]->Fill(TrkPar->fTrack->fPSTBack, Weight);
    Hist->fPSTDiff->Fill(TrkPar->fTrack->fPSTBack - TrkPar->fTrack->fP, Weight);
    Hist->fPSTApproxDiff->Fill(TrkPar->fApproxDpST, Weight);
    Hist->fPExitDiff->Fill(TrkPar->fTrack->fP - TrkPar->fTrack->fPTrackerExit, Weight);
    Hist->fT0->Fill(TrkPar->fTrack->fT0, Weight);
    Hist->fT0Err->Fill(TrkPar->fTrack->fT0Err, Weight);
    Hist->fD0->Fill(TrkPar->fTrack->fD0, Weight);
    Hist->fDP->Fill((TrkPar->fTrack->fP) - (TrkPar->fTrack->fPFront), Weight);
    Hist->fDPvsP->Fill(TrkPar->fTrack->fPFront, (TrkPar->fTrack->fP) - (TrkPar->fTrack->fPFront), Weight);
    Hist->fDPvsNH->Fill(TrkPar->fTrack->NActive(), (TrkPar->fTrack->fP) - (TrkPar->fTrack->fPFront), Weight);
    Hist->fDPCorr->Fill((TrkPar->fTrack->fP + TrkPar->fApproxDpST) - (TrkPar->fTrack->fPStOut), Weight);
    Hist->fChi2NDof->Fill(TrkPar->fTrack->Chi2Dof(), Weight);
    Hist->fFitCons[0]->Fill(TrkPar->fTrack->fFitCons, Weight);
    Hist->fFitCons[1]->Fill(std::log10(std::max(1.e-10f, TrkPar->fTrack->fFitCons)), Weight);
    Hist->fFitMomErr->Fill(TrkPar->fTrack->fFitMomErr, Weight);
    Hist->fTanDip->Fill(TrkPar->fTrack->fTanDip, Weight);
    Hist->fCosTheta->Fill(TrkPar->CosTheta(), Weight);
    Hist->fRadius->Fill(TrkPar->fRadius, Weight);
    Hist->fRMax->Fill(TrkPar->fRMax, Weight);
    Hist->fNActive->Fill(TrkPar->fTrack->NActive(), Weight);
    Hist->fNActiveFrac->Fill(TrkPar->fTrack->NActive() * 1. / TrkPar->fTrack->NHits(), Weight);
    Hist->fTrkQual[0]->Fill(TrkPar->fTrack->fTrkQual, Weight);
    Hist->fTrkQual[1]->Fill(TrkPar->fTrkQual, Weight);
    Hist->fPID[0]->Fill(TrkPar->fPID, Weight);
    Hist->fPID[1]->Fill(TrkPar->fTrkPID, Weight);
    Hist->fCosmicID->Fill(TrkPar->fCosmicID, Weight);
    Hist->fClusterE->Fill(TrkPar->fTrack->fClusterE, Weight);
    Hist->fClusterDisk->Fill(TrkPar->fTrack->fDiskID, Weight);
    Hist->fDt->Fill(TrkPar->fTrack->fDt, Weight);
    Hist->fEp->Fill(TrkPar->fTrack->fEp, Weight);
    Hist->fTZSlope->Fill(TrkPar->fTZSlope, Weight);
    Hist->fTZSlopeSig->Fill(TrkPar->TZSlopeSig(), Weight);
    Hist->fTZSlopeRatio->Fill(TrkPar->TZSlopeRatio(), Weight);
    Hist->fBestAlg->Fill(TrkPar->fTrack->BestAlg(), Weight);
    Hist->fAlgMask->Fill(TrkPar->fTrack->AlgMask(), Weight);
    Hist->fSTBoundary->Fill((TrkPar->STBoundary()) ? 1. : 0., Weight);
    Hist->fSTInters->Fill(TrkPar->fTrack->NSTIntersections(), Weight);
    Hist->fIPAInters->Fill(TrkPar->fTrack->NIPAIntersections(), Weight);
    Hist->fOPAInters->Fill(TrkPar->fTrack->NOPAIntersections(), Weight);
    if(TrkPar->fIDWord[0] == 0) {
      Hist->fTrackID->Fill("Passed", Weight);
      Hist->fExlTrackID->Fill("Passed", Weight);
    } else {
      for(int bit = 0; bit < 30; ++bit) {
        if((TrkPar->fIDWord[0] & (1 << bit)) != 0) {
          TString name((kTrackIDNames.count(bit)) ? kTrackIDNames[bit] : Form("Unknown-%i", bit));
          Hist->fTrackID->Fill(name.Data(), Weight);
          if((TrkPar->fIDWord[0] & ~(1 << bit)) == 0)
            Hist->fExlTrackID->Fill(name.Data(), Weight);
        }
      }
    }

    auto crv_stub_par = TrkPar->fCRVStubPar;
    if(crv_stub_par) {
      auto crv_stub = crv_stub_par->fCluster;
      const float dt = TrkPar->fTrack->fT0 - crv_stub_par->fCorrTime;
      const float dt_crv = TrkPar->fTrack->fT0 - crv_stub_par->fTime;
      const float dt_extrap = TrkPar->fTrack->fT0 - crv_stub_par->fApproxTimeExtrapToFront;
      const float deltat_st = TrkPar->CRVSTDeltaT();
      const float deltat_calo = TrkPar->CRVCaloFrontDeltaT(false);     // TrkPar->CRVCaloDeltaT();
      const float deltat_calo_muon = TrkPar->CRVCaloFrontDeltaT(true); // TrkPar->CRVCaloDeltaT();
      // const float deltat_extrap = TrkPar->fTrack->fT0 - crv_stub_par->fApproxTimeExtrapToFront;

      // Minimum of the two hypotheses
      const float min_deltat = TrkPar->CRVMinDeltaT();

      // Use the Z extrapolation to determine the best time to use
      const float min_extrap = (crv_stub_par->fExtrapZ < 7000.)  ? deltat_st
                               : (crv_stub_par->fExtrapZ > 1.e5) ? deltat_calo // ignore only calo hypothesis for now
                                                                 : min_deltat;

      Hist->fCRVDeltaT->Fill(dt, Weight);
      Hist->fCRVDeltaTCRV->Fill(dt_crv, Weight);
      Hist->fCRVDeltaTST->Fill(deltat_st, Weight);
      Hist->fCRVDeltaTCalo[0]->Fill(deltat_calo, Weight);
      Hist->fCRVDeltaTCalo[1]->Fill(deltat_calo_muon, Weight);
      Hist->fCRVDeltaTExtrap->Fill(min_extrap /*deltat_extrap*/, Weight);
      Hist->fCRVMinDeltaT->Fill(min_deltat, Weight);
      Hist->fCRVExtrapZ->Fill(crv_stub_par->fExtrapZ, Weight);
      const float x = crv_stub->Position()->X();
      const float y = crv_stub->Position()->Y();
      const float z = crv_stub->Position()->Z();
      Hist->fCRVXZ->Fill(z, x, Weight);
      Hist->fCRVYZ->Fill(z, y, Weight);
      Hist->fCRVdTZ->Fill(z, dt, Weight);
      Hist->fCRVdTZCRV->Fill(z, dt_crv, Weight);
      Hist->fCRVdTZExtrap->Fill(z, dt_extrap, Weight);
    }

    auto us_trk = TrkPar->fUpstreamTrack;
    if(us_trk) {
      Hist->fUpstreamDt->Fill(TrkPar->fTrack->fT0 - us_trk->fT0, Weight);
      Hist->fUpstreamDp->Fill(TrkPar->fTrack->fP - us_trk->fP, Weight);
    }

    Hist->fMCPFront->Fill(TrkPar->fTrack->fPFront, Weight);
    Hist->fMCPStOut->Fill(TrkPar->fTrack->fPStOut, Weight);
    Hist->fMCPStDiff->Fill(TrkPar->fTrack->fPStOut - TrkPar->fTrack->fPFront, Weight);
    if(TrkPar->fTrack->fPStOut > 0. && TrkPar->fTrack->fPSTBack > 0.) { // only meaningful if all info is available
      const float mc_pstdiff = TrkPar->fTrack->fPStOut - TrkPar->fTrack->fPFront;
      const float reco_pstdiff = TrkPar->fTrack->fPSTBack - TrkPar->fTrack->fP;
      Hist->fMCPStDiffDiff->Fill(reco_pstdiff - mc_pstdiff, Weight);
    } else
      Hist->fMCPStDiffDiff->Fill(-999., Weight);
    if(TrkPar->fTrack->fPStOut > 0.) {
      const float mc_pstdiff = TrkPar->fTrack->fPStOut - TrkPar->fTrack->fPFront;
      const float reco_pstdiff = TrkPar->fApproxDpST;
      Hist->fMCApproxPStDiffDiff->Fill(reco_pstdiff - mc_pstdiff, Weight);
    } else
      Hist->fMCApproxPStDiffDiff->Fill(-999., Weight);
    Hist->fMCGenE->Fill(TrkPar->fGenE, Weight);
    Hist->fMCPGenEDiff->Fill(TrkPar->fTrack->P() - TrkPar->fGenE, Weight);
    Hist->fMCPSig->Fill((TrkPar->fTrack->fFitMomErr > 0.) ? (TrkPar->fTrack->fP - TrkPar->fTrack->fPFront) / TrkPar->fTrack->fFitMomErr : -999., Weight);
    Hist->fMCPdg[0]->Fill(TrkPar->fTrack->fPdgCode, Weight);
    Hist->fMCPdg[1]->Fill(std::abs(TrkPar->fTrack->fPdgCode), Weight);
    Hist->fMCStrawHits->Fill(TrkPar->fTrack->NMcStrawHits(), Weight);
    Hist->fMCGoodHits->Fill(TrkPar->fTrack->NGoodMcHits(), Weight);
    Hist->fMCTrajectory->Fill(TrkPar->fTrack->fMcDirection, Weight);
    Hist->fMCSimProc->Fill((TrkPar->fSimp) ? TrkPar->fSimp->CreationCode() : 0, Weight);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::FillHelixHistograms(HelixHist_t* Hist, HelixPar_t* HlxPar, float Weight) {
    if(!Hist)
      throw std::runtime_error(Form("TAnaModule::%s: Attempting to fill undefined histograms\n", __func__));
    if(!HlxPar)
      throw std::runtime_error(Form("TAnaModule::%s: Attempting to fill histograms with undefinded helix par\n", __func__));
    auto Helix = HlxPar->fHelix;

    const int nhits = Helix->NHits();
    const double clusterT = Helix->ClusterTime();
    const double clusterE = Helix->ClusterEnergy();

    const double radius = Helix->Radius();

    const double lambda = Helix->Lambda();
    const double tanDip = lambda / radius;
    const double pT = Helix->Pt();
    int genID = Helix->fSimpPDG1;
    if(abs(genID) == 2212)
      genID = genID * 25 * ((genID < 0) ? -1 : 1); // compress some particle IDs
    if(abs(genID) == 211)
      genID = genID * 26 * ((genID < 0) ? -1 : 1);
    const double gen_r = std::sqrt(std::pow(Helix->fSimpOrigin1.X() + 3904., 2) + std::pow(Helix->fSimpOrigin1.Y(), 2));
    const double gen_z = Helix->fSimpOrigin1.Z();
    const double gen_p = Helix->fSimpMom1.P();
    const double gen_pt = Helix->fSimpMom1.Pt();
    const double p = Helix->P();

    Hist->fHelicity->Fill(Helix->Helicity(), Weight);
    Hist->fNHits->Fill(nhits, Weight);
    Hist->fClusterTime->Fill(clusterT, Weight);
    Hist->fClusterEnergy->Fill(clusterE, Weight);

    Hist->fRadius->Fill(radius, Weight);
    Hist->fRMax->Fill(HlxPar->fRMax, Weight);
    Hist->fMom->Fill(p, Weight);
    Hist->fPt->Fill(pT, Weight);
    Hist->fLambda->Fill(lambda, Weight);
    Hist->fTanDip->Fill(tanDip, Weight);

    Hist->fGenMom->Fill(gen_p, Weight);
    Hist->fGenID->Fill(genID, Weight);
    Hist->fGenRZ->Fill(gen_z, gen_r, Weight);
    Hist->fDp->Fill(p - gen_p, Weight);
    Hist->fDpT->Fill(Helix->Pt() - gen_pt, Weight);

    Hist->fBestAlg->Fill(Helix->BestAlg(), Weight);
    Hist->fAlgMask->Fill(Helix->AlgMask(), Weight);
    Hist->fD0->Fill(Helix->D0(), Weight);
    Hist->fT0->Fill(Helix->T0(), Weight);
    Hist->fT0Err->Fill(Helix->T0Err(), Weight);
    Hist->fChi2XY->Fill(Helix->Chi2XY(), Weight);
    Hist->fChi2ZPhi->Fill(Helix->Chi2ZPhi(), Weight);
    Hist->fTZSlope->Fill(Helix->TZSlope(), Weight);
    Hist->fTZSlopeErr->Fill(Helix->TZSlopeError(), Weight);
    Hist->fTZSlopeSig->Fill((Helix->TZSlopeError() > 0.) ? Helix->TZSlope() / Helix->TZSlopeError() : -999., Weight); // sign the significance by the direction
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::FillClusterHistograms(ConvAna::ClusterHist_t* Hist, ConvAna::ClusterPar_t* Par, float Weight) {
    if(!Par)
      return;
    TStnCluster* Cluster = Par->fCluster;
    if(!Cluster)
      return;

    int row = Cluster->Ix1();
    int col = Cluster->Ix2();
    const float x = Cluster->fX;
    const float y = Cluster->fY;
    const float z = Cluster->fZ;
    const float r = sqrt(x * x + y * y);
    const float energy = Cluster->Energy();
    const float time = Cluster->Time();
    const float core_energy = Par->core_energy();

    if((row < 0) || (row > 9999))
      row = -9999;
    if((col < 0) || (col > 9999))
      col = -9999;

    Hist->fDiskID->Fill(Cluster->DiskID(), Weight);
    Hist->fEnergy->Fill(energy, Weight);
    Hist->fT0->Fill(time, Weight);
    Hist->fRow->Fill(row, Weight);
    Hist->fCol->Fill(col, Weight);
    Hist->fX->Fill(x, Weight);
    Hist->fY->Fill(y, Weight);
    Hist->fZ->Fill(z, Weight);
    Hist->fR->Fill(r, Weight);

    Hist->fYMean->Fill(Cluster->fYMean, Weight);
    Hist->fZMean->Fill(Cluster->fZMean, Weight);
    Hist->fSigY->Fill(Cluster->fSigY, Weight);
    Hist->fSigZ->Fill(Cluster->fSigZ, Weight);
    Hist->fSigR->Fill(Cluster->fSigR, Weight);
    Hist->fNCr0->Fill(Cluster->fNCrystals, Weight);
    Hist->fNCr1->Fill(Cluster->fNCr1, Weight);
    Hist->fFrE1->Fill(Cluster->fFrE1, Weight);
    Hist->fFrE2->Fill(Cluster->fFrE2, Weight);
    Hist->fSigE1->Fill(Cluster->fSigE1, Weight);
    Hist->fSigE2->Fill(Cluster->fSigE2, Weight);
    Hist->fTimeRMS->Fill(Cluster->fTimeRMS, Weight);
    Hist->fMaxR->Fill(Cluster->fMaxR, Weight);
    Hist->fE9OverE->Fill(Cluster->fE9 / energy, Weight);
    Hist->fE25OverE->Fill(Cluster->fE25 / energy, Weight);
    Hist->fRingEOverE->Fill(Cluster->RingE() / energy, Weight);
    Hist->fRingEOverE1->Fill(Cluster->RingE() / Cluster->E1(), Weight);
    Hist->fOutRingE->Fill(Cluster->fOutRingE, Weight);
    Hist->fOutRingEOverE->Fill(Cluster->fOutRingE / energy, Weight);
    Hist->fNCoreCrystals->Fill(Par->n_core_crystals(), Weight);
    Hist->fCoreEnergy->Fill(core_energy, Weight);
    Hist->fCoreEnergyFrac->Fill(core_energy / energy, Weight);

    // MC info
    const float mc_edep = Cluster->fMCEDep;
    const float mc_time = Cluster->fMCTime;
    const auto sim = Par->fSim;
    const float gen_energy = (sim) ? sim->fStartMom.E() : 0.;
    Hist->fMCSimEDep->Fill(Cluster->fMCSimEDep, Weight);
    Hist->fMCSimMomIn->Fill(Cluster->fMCSimMomIn, Weight);
    Hist->fMCSimPdg->Fill(Cluster->fMCSimPDG, Weight);
    Hist->fMCSimPdgName->Fill(NameFromPDG(Cluster->fMCSimPDG).Data(), Weight);
    Hist->fMCSimEStart->Fill(gen_energy, Weight);
    Hist->fMCEDep->Fill(mc_edep, Weight);
    Hist->fMCTime->Fill(mc_time, Weight);
    Hist->fMC_dE->Fill(energy - mc_edep, Weight);
    Hist->fMC_dt->Fill(time - mc_time, Weight);
    Hist->fMC_dGenE->Fill(energy - gen_energy, Weight);
    Hist->fMCGenE->Fill(gen_energy, Weight);
    Hist->fMCEvsE->Fill(mc_edep, energy, Weight);
    Hist->fMCGenEvsE->Fill(gen_energy, energy, Weight);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::FillCRVClusterHistograms(ConvAna::CRVClusterHist_t* Hist, ConvAna::CRVStubPar_t* CrvStubPar, float Weight) {
    if(!Hist)
      throw std::runtime_error(Form("TAnaModule::%s: Histogram set undefined!", __func__));
    if(!CrvStubPar)
      throw std::runtime_error(Form("TAnaModule::%s: CRV stub par undefined!", __func__));
    TCrvCoincidenceCluster* CrvCluster = CrvStubPar->fCluster;
    if(!CrvCluster)
      throw std::runtime_error(Form("TAnaModule::%s: Cluster undefined!", __func__));

    const float width = CrvCluster->EndTime() - CrvCluster->StartTime();
    const float x = CrvCluster->Position()->X();
    const float y = CrvCluster->Position()->Y();
    const float z = CrvCluster->Position()->Z();

    Hist->fSector->Fill(CrvStubPar->fSector, Weight);
    Hist->fFirstBar->Fill(CrvStubPar->fFirstBar, Weight);
    Hist->fNPulses->Fill(CrvCluster->NPulses(), Weight);
    Hist->fNPe->Fill(CrvCluster->NPe(), Weight);
    Hist->fNPePP->Fill(CrvStubPar->fNPePP, Weight);
    Hist->fStartTime->Fill(CrvCluster->StartTime(), Weight);
    Hist->fEndTime->Fill(CrvCluster->EndTime(), Weight);
    Hist->fWidth->Fill(width, Weight);
    Hist->fXVsZ->Fill(z, x, Weight);
    Hist->fYVsZ->Fill(z, y, Weight);
    Hist->fCorrTime->Fill(CrvStubPar->fCorrTime, Weight);
    Hist->fCorrTimeProp->Fill(CrvStubPar->fCorrTimeProp, Weight);
    Hist->fCorrTimeToF->Fill(CrvStubPar->fCorrTimeTof, Weight);
    Hist->fApproxTimeST->Fill(CrvStubPar->fApproxTimeST, Weight);
    Hist->fApproxTimeCalo->Fill(CrvStubPar->fApproxTimeCalo, Weight);
    Hist->fApproxTimeExtrap->Fill(CrvStubPar->fApproxTimeExtrap, Weight);
    Hist->fApproxTimeSTToFront->Fill(CrvStubPar->fApproxTimeSTToFront, Weight);
    Hist->fApproxTimeCaloToFront->Fill(CrvStubPar->fApproxTimeCaloToFront, Weight);
    Hist->fApproxTimeExtrapToFront->Fill(CrvStubPar->fApproxTimeExtrapToFront, Weight);
    Hist->fBarsOneEnd->Fill(CrvStubPar->fTotalBars - CrvStubPar[0].fTwoEndBars, Weight);
    Hist->fCrvPropdT->Fill(CrvStubPar->fCorrTimeProp - CrvCluster->StartTime(), Weight);
    Hist->fBarsTwoEnd->Fill(CrvStubPar->fTwoEndBars, Weight);
    Hist->fNSectors->Fill(CrvStubPar->fNSectors, Weight);
    Hist->fNDiffLSectors->Fill(CrvStubPar->fNDiffLSectors, Weight);
    Hist->fStubSlope->Fill(CrvStubPar->fStubDYDZ, Weight);
    Hist->fStubSlopeChi2->Fill(CrvStubPar->fStubSlopeChi2, Weight);
    Hist->fStubSlopeDelta->Fill(CrvStubPar->fStubDYDZ - CrvStubPar[0].fStubDYDZMC, Weight);
    Hist->fStubQN->Fill(CrvStubPar->fStubQN, Weight);
    Hist->fStubSlopeMCProduct->Fill(CrvStubPar->fStubSlopeMCProduct, Weight);
    Hist->fMCTime->Fill(CrvCluster->fMCAvgHitTime, Weight);
    Hist->fMCdT[0]->Fill(CrvCluster->StartTime() - CrvCluster->fMCAvgHitTime, Weight);
    Hist->fMCdT[1]->Fill(CrvStubPar->fTime - CrvCluster->fMCAvgHitTime, Weight);
  }

  void TAnaModule::InitEventInfo() { SetTriggersPassed(); }

  void TAnaModule::InitTrackPar(TStnTrack* Trk, ConvAna::TrackPar_t* TrkPar, TStnHelix* Hlx) {
    TrkPar->fCluster = nullptr;

    // clear track selection flags
    for(int id = 0; id < kTrackIDs; ++id)
      TrkPar->fIDWord[id] = 0;
    for(int ifit = 0; ifit < kMaxTrackFits; ++ifit)
      TrkPar->fAltHypotheses[ifit] = nullptr;
    TrkPar->fNAlt = 0;

    TrkPar->fCRVStubPar = nullptr;
    TrkPar->fUpstreamTrack = nullptr;
    TrkPar->fMCCRVStubPar = nullptr;

    // set pointer to TStnTrack
    TrkPar->fTrack = Trk;
    TrkPar->fHelix = Hlx;
    if(!Trk)
      return;

    // default to momentum as observable
    TrkPar->fObs = Trk->fP;

    // compute parameters not contained in TStnTrack
    TrkPar->fRadius = Trk->fPt / (mmTconversion * bz0);
    TrkPar->fRMax = -1 * Trk->fCharge * Trk->fD0 + 2 * TrkPar->fRadius;
    // pz/pt = tan(dip) --> tan(dip) = 0 = vertical
    //  length along IPA per intersection = c, a = thickness, and tan(dip) = b/a --> c = sqrt(a^2 +
    //  a^2*tan(dip)) = a*sqrt(1 + tan^2(dip)) TrkPar->fApproxDpST =
    //  0.10*Trk->NIPAIntersections()*std::sqrt(1. + std::pow(Trk->fTanDip, 2));
    TrkPar->fApproxDpST = 0.14 * Trk->NIPAIntersections(); // ignore the angle, assume just an average loss per intersection

    // Take the TZ slope from the helix for now
    if(Hlx) {
      TrkPar->fTZSlope = Hlx->fTZSlope;
      TrkPar->fTZSlopeErr = Hlx->fTZSlopeError;
    } else {
      TrkPar->fTZSlope = 0.f;
      TrkPar->fTZSlopeErr = 0.f;
    }

    // systematic shifts
    const float mom_scale = 1e-4f;
    TrkPar->fPUp = (1.f + mom_scale) * Trk->P();
    TrkPar->fPDown = (1.f - mom_scale) * Trk->P();

    TrkPar->fGenE = 0.f;
    TrkPar->fSimp = nullptr;
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::InitHelixPar(TStnHelix* Hlx, ConvAna::HelixPar_t* HlxPar) {

    // set pointer to TStnTrack
    HlxPar->fHelix = Hlx;
    if(!Hlx)
      return;

    // compute parameters not contained in TStnHelix
    HlxPar->fRMax = Hlx->D0() + 2. * Hlx->Radius();

    // truth-level info
    HlxPar->fIsMCDownstream = Hlx->fSimpMom1.Pz() > 0.; // sim particle with most hits is moving downstream
    HlxPar->fTZSigMC = (Hlx->fTZSlope / Hlx->fTZSlopeError) * ((HlxPar->fIsMCDownstream) ? 1.f : -1.f);
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::InitClusterPar(TStnCluster* Cluster, ConvAna::ClusterPar_t* ClusterPar) {
    if(!ClusterPar)
      return;
    ClusterPar->init(Cluster);
    ClusterPar->fSim = nullptr;
    if(!Cluster)
      return;
  }

  //-----------------------------------------------------------------------------
  void TAnaModule::InitCRVStubPar(TCrvClusterBlock* CrvClusterBlock, ConvAna::CRVStubPar_t* CrvStubPar, int maxStubs, TSimParticle* Simp) {
    if(!CrvClusterBlock || !CrvStubPar) {
      if(fDebugLevel > 0)
        printf("TAnaModule::%s: Block/Par not defined, exiting\n", __func__);
      return;
    }

    const int nstubs = std::min(CrvClusterBlock->NClusters(), maxStubs);

    if(fDebugLevel > 0)
      printf("---- N(CRV stubs) = %2i\n", nstubs);
    for(int is = 0; is < nstubs; is++) {
      if(fDebugLevel > 0)
        printf(" ---- new stub: index = %2i\n", is);
      TCrvCoincidenceCluster* crv_stub = CrvClusterBlock->Cluster(is);
      if(fDebugLevel > 2)
        printf(" Retrieved the CRV cluster\n");
      ConvAna::CRVStubPar_t& sp = CrvStubPar[is];
      sp.reset();
      if(!crv_stub) {
        printf(
            "TAnaModule::%s: WARNING: CRV cluster %i is null (NClusters = %i, TClonesArray entries "
            "= "
            "%i)\n",
            __func__, is, CrvClusterBlock->NClusters(), CrvClusterBlock->GetListOfClusters()->GetEntriesFast()
        );
        continue;
      }
      if(fDebugLevel > 2)
        printf(" Retrieved the CRVStubPar_t struct\n");
      sp.fCluster = crv_stub;

      int np = CrvClusterBlock->NClusterPulses(is);
      if(fDebugLevel > 2)
        printf(" The CRV cluster has %i pulses\n", np);
      sp.fStubDYDZ = crv_stub->Slope();
      // Time correction
      constexpr float crv_time_correction = -17.5f; // ~constant offset between start time and MC time
      // Positions
      const float x_cl(crv_stub->Position()->X()), y_cl(crv_stub->Position()->Y()), z_cl(crv_stub->Position()->Z());
      const float x_st(-3904.f), y_st(0.f), z_st(6271.f);     // stopping target exit
      const float x_cal(-3904.f), y_cal(0.f), z_cal(11820.f); // calo disk 0 front on solenoid axis
      // Evaluate time from stub dy/dz slope estimate to the solenoid axis and then to the tracker
      // front z = (y_1 - y_0) * dz/dy + z_0; y_1 = 0 = solenoid axis
      const float z_exr((sp.fStubDYDZ != 0.) ? z_cl + (y_st - y_cl) / sp.fStubDYDZ : (z_cl < z_st + 0.5 * (z_st + z_cal)) ? z_st : z_cal), x_exr(x_st), y_exr(y_st);

      // Distances
      const float dx_st(x_st - x_cl), dy_st(y_st - y_cl), dz_st(z_st - z_cl);
      const float dx_cal(x_cal - x_cl), dy_cal(y_cal - y_cl), dz_cal(z_cal - z_cl);
      const float dx_exr(x_exr - x_cl), dy_exr(y_exr - y_cl), dz_exr(z_exr - z_cl);

      // Times
      const float vlightinv = 1.f / 300.f; // light travels ~300 mm / ns
      const float approx_tof_st = vlightinv * std::sqrt(dx_st * dx_st + dy_st * dy_st + dz_st * dz_st);
      const float approx_tof_cal = vlightinv * std::sqrt(dx_cal * dx_cal + dy_cal * dy_cal + dz_cal * dz_cal);
      const float approx_tof_exr = vlightinv * std::sqrt(dx_exr * dx_exr + dy_exr * dy_exr + dz_exr * dz_exr);

      sp.fTime = crv_stub->StartTime() + crv_time_correction;
      sp.fZ = z_cl;
      sp.fExtrapZ = z_exr;
      sp.fApproxTimeST = sp.fTime + approx_tof_st;
      sp.fApproxTimeCalo = sp.fTime + approx_tof_cal;
      sp.fApproxTimeExtrap = sp.fTime + approx_tof_exr;
      sp.fApproxTimeSTToFront = sp.fApproxTimeST + 23.f; // account for rebounding etc.
      sp.fApproxTimeCaloToFront = sp.fApproxTimeCalo + 87.f;
      sp.fApproxTimeExtrapToFront = sp.fTime + approx_tof_exr + vlightinv * (z_exr - z_st) + 23.f;
      // FIXME: Using approximate time estimates for now
      sp.fCorrTime = (std::fabs(z_cl - z_st) < std::fabs(z_cl - z_cal)) ? sp.fApproxTimeSTToFront : sp.fApproxTimeCaloToFront;
      sp.fTotalBars = 0;
      sp.fNPePP = (crv_stub->NPulses() > 0) ? (float)crv_stub->NPe() / crv_stub->NPulses() : 0.f;
      if(fDebugLevel > 2)
        printf(" Finished processing stub %i\n", is);
    }
  }

  //_____________________________________________________________________________
  TStnTrack* TAnaModule::MatchUpstreamTrack(TStnTrack* track, TStnTrackBlock* tracks) {
    if(!track || !tracks)
      return nullptr;
    if(track->fMomentum.Pz() < 0.)
      return nullptr; // can't match an upstream leg to an upstream track

    const int ntrks = tracks->NTracks();
    TStnTrack* match = nullptr;
    const float tmin(50.), tmax(150.); // reflections are typically within this time window
    for(int itrk = 0; itrk < ntrks; ++itrk) {
      auto u_trk = tracks->Track(itrk);
      const double dt = track->fT0 - u_trk->fT0;
      if(dt > tmin && dt < tmax) { // reasonable time window
        const double current_dt = (match) ? track->fT0 - match->fT0 : dt + 1000.;
        if(current_dt > dt)
          match = u_trk; // take the closest in time
      }
    }
    return match;
  }

  //_____________________________________________________________________________
  ConvAna::CRVStubPar_t* TAnaModule::MatchCRVToTrack(TStnTrack* track, TCrvClusterBlock* clusters, ConvAna::CRVStubPar_t* stubPars, int maxStubs) {
    ConvAna::CRVStubPar_t* match(nullptr);
    if(!clusters || !stubPars || !track)
      return match;
    float match_t_cl(0.f), match_t_tk(0.f);

    const static float min_time_window(-250.), max_time_window(300.); // maximum delta t to consider
    const static int min_pe(0), min_pepp(0);                          // reduce impact from pileup clusters

    const int ncl = std::min(clusters->NClusters(), maxStubs);
    for(int icl = 0; icl < ncl; ++icl) {
      auto& par = stubPars[icl];
      // Defensive checks: verify cluster pointer is valid and properly initialized
      if(!par.fCluster)
        continue;

      // Additional safety check: ensure NPe() is called only on valid cluster
      const int npe = par.fCluster->NPe();
      if(npe < 0)
        continue; // sanity check for corrupted data

      const float t_cl = par.fCorrTime;
      const float t_tk = track->fT0;
      const float delta_t = t_tk - t_cl;
      // a reasonable match
      if(delta_t < max_time_window && delta_t > min_time_window && npe > min_pe && par.fNPePP > min_pepp) {
        // check if this match is better than the current one (if there's a current one)
        bool replace = !match;
        replace |= std::fabs(delta_t) < std::fabs(match_t_tk - match_t_cl);
        if(replace) {
          match = &par; // accept this one as the best match so far
          match_t_cl = t_cl;
          match_t_tk = t_tk;
        }
      }
    }
    return match;
  }

  //_____________________________________________________________________________
  void TAnaModule::Debug() {
    if(GetDebugBit(0)) {
      auto event = GetEvent();
      printf(">>> Event %5i/%5i/%6i:\n", event->fRunNumber, event->fSectionNumber, event->fEventNumber);
    }
  }

  //-----------------------------------------------------------------------------
  // Following SU2020 strategy, see su2020/ana/TAnaModule
  // Set cosmic veto bits based on the tracker and calorimeter information, without looking at the
  // CRV
  // FIXME: Updates needed:
  //        - Switch to actual upstream reco results to veto reflecting cosmics
  //        - Switch to actual muon reco results for muon veto
  //-----------------------------------------------------------------------------
  int TAnaModule::NonCRVCosmicVeto(CosmicVetoData_t* Data) {
    int VetoID(0);
    if(!Data)
      return VetoID;

    //-----------------------------------------------------------------------------
    // vetoing cosmics
    // 1. veto events with more than one track - upstream or downstream
    //    the bits are specified in EventPar_t.h
    //-----------------------------------------------------------------------------

    int ntrk_de = (Data->fTrackBlockDe) ? Data->fTrackBlockDe->NTracks() : 0;
    int trk_veto_de = 0;
    TrackPar_t trkpar_1, trkpar_2;

    for(int i1 = 0; i1 < ntrk_de; i1++) {
      TStnTrack* t1 = Data->fTrackBlockDe->Track(i1);
      InitTrackPar(t1, &trkpar_1);
      if(trkpar_1.fIDWord[1])
        continue;
      //-----------------------------------------------------------------------------
      // in principle, this is not needed. But, so far, can't rely on the upstream reco,
      // so call events with multiple tracks "cosmic candidates"
      //-----------------------------------------------------------------------------
      for(int i2 = 0; i2 < ntrk_de; i2++) {
        if(i2 == i1)
          continue;
        // Use a looser selection for the additional track veto selection: just chi^/dof and delta t
        TStnTrack* t2 = Data->fTrackBlockDe->Track(i2);
        InitTrackPar(t2, &trkpar_2);
        if(t2->Chi2Dof() > 3.)
          continue;
        float dt = t1->T0() - t2->T0();
        if(fabs(dt) < 30)
          continue;
        trk_veto_de += 1;
      }
    }

    int ntrk_ue = (Data->fTrackBlockUe) ? Data->fTrackBlockUe->NTracks() : 0;
    int trk_veto_ue = 0;
    for(int i1 = 0; i1 < ntrk_de; i1++) {
      TStnTrack* t1 = Data->fTrackBlockDe->Track(i1);
      // consider only good tracks
      InitTrackPar(t1, &trkpar_1);
      if(trkpar_1.fIDWord[1] != 0)
        continue;
      for(int i2 = 0; i2 < ntrk_ue; i2++) {
        TStnTrack* t2 = Data->fTrackBlockUe->Track(i2);
        InitTrackPar(t2, &trkpar_2);
        float dt = t1->T0() - t2->T0();
        if(fabs(dt) < 30) {
          // Ue track, close in time... could be a muon ... no muon reco so far
          if((t1->Charge() * t2->Charge() < 0) && (fabs(t1->P() - t2->P0()) < 5.)) {

            // feels like the same track - in reality, need to use the hit info

            if((t2->Chi2Dof() < t1->Chi2Dof()) && (t2->NActive() > t1->NActive())) {

              // upstream track looks better
              trk_veto_ue += 1;
            }
          }
        } else {
          // two tracks, an upstream one and a downstream one, far in time,
          // check quality of the upstram track
          if(t2->Chi2Dof() < 3.) {
            // the upstream track looks real
            trk_veto_ue += 1;
          }
        }
      }
    }

    if(trk_veto_de > 0)
      VetoID |= ConvAna::kNTrkDeVetoBit;
    if(trk_veto_ue > 0)
      VetoID |= ConvAna::kNTrkUeVetoBit;

    int ntc_ue = (Data->fTCFinderBlockUe) ? Data->fTCFinderBlockUe->NTimeClusters() : 0;
    int nhel_de = (Data->fHelixBlockDe) ? Data->fHelixBlockDe->NHelices() : 0;
    int nhel_ue = (Data->fHelixBlockUe) ? Data->fHelixBlockUe->NHelices() : 0;

    int tc_ue_veto = 0;
    int trk_ue_veto = 0;
    int helix_de_veto = 0;
    int helix_ue_veto = 0;

    for(int i = 0; i < ntrk_de; i++) {
      TStnTrack* trk = Data->fTrackBlockDe->Track(i);
      InitTrackPar(trk, &trkpar_1);

      float trk_t0 = trk->T0();
      //-----------------------------------------------------------------------------
      // 2. veto events with the upstream tracks not identical to downstream tracks
      //-----------------------------------------------------------------------------
      for(int i1 = 0; i1 < ntrk_ue; i1++) {
        TStnTrack* trk_ue = Data->fTrackBlockUe->Track(i1);
        float tu = trk_ue->T0();

        float dt = trk_t0 - tu;

        if((dt > 50) && (dt < 200)) {
          //-----------------------------------------------------------------------------
          // downstream and upstream tracks are separated by more than 50 ns,
          // the downstream track comes later
          //-----------------------------------------------------------------------------
          trk_ue_veto += 1;
        }
      }
      //-----------------------------------------------------------------------------
      // 3. assume an upstream leg has not been found - look for the corresponding time
      //    cluster and veto events with TCFinderUe clusters 50-200 ns earlier than
      //    the track in question
      //    count "offending" time clusters, but set just one bit
      //    consider only timeclusters with N(combo hits) > 25
      //-----------------------------------------------------------------------------
      for(int itc = 0; itc < ntc_ue; itc++) {
        TStnTimeCluster* tc = Data->fTCFinderBlockUe->TimeCluster(itc);
        if(tc->NComboHits() <= 25)
          continue;
        float dt = trk_t0 - tc->T0();
        if((dt > 50.) && (dt < 200)) {
          tc_ue_veto += 1;
        }
        if(fabs(dt) < 50.) {
          trkpar_1.fDNhitsUe = trk->NHits() - tc->NHits();
        }
      }

      //-----------------------------------------------------------------------------
      // 4. veto using the number of downstream and upstream helices
      //    helices to have Pt > 50 MeV and be far in time from the De track
      //-----------------------------------------------------------------------------
      for(int i2 = 0; i2 < nhel_de; i2++) {
        TStnHelix* hel = Data->fHelixBlockDe->Helix(i2);
        if(hel->NComboHits() <= 10)
          continue;
        if(hel->Pt() < 50)
          continue;
        float dt = trk_t0 - hel->T0();
        if((dt > 50.) && (dt < 200)) {
          helix_de_veto += 1;
        }
      }

      for(int i2 = 0; i2 < nhel_ue; i2++) {
        TStnHelix* hel = Data->fHelixBlockUe->Helix(i2);
        if(hel->NComboHits() <= 10)
          continue;
        if(hel->Pt() < 50)
          continue;
        float dt = trk_t0 - hel->T0();
        if((dt > 50.) && (dt < 200)) {
          helix_ue_veto += 1;
        }
      }
    }

    if(trk_ue_veto > 0)
      VetoID |= ConvAna::kTrkUeDtVetoBit;
    if(tc_ue_veto > 0)
      VetoID |= ConvAna::kTrkTcVetoBit;
    if(helix_de_veto > 0)
      VetoID |= ConvAna::kNHelDeVetoBit;
    if(helix_ue_veto > 0)
      VetoID |= ConvAna::kNHelUeVetoBit;
    //-----------------------------------------------------------------------------
    // 5. veto too energetic clusters in the event close to the track of interest
    //    15 ns is OK, 140 MeV may need to be tuned. For now, set the energy threshold
    //    just above max RPC energy
    //-----------------------------------------------------------------------------
    int ncalo_cl = (Data->fClusterBlock) ? Data->fClusterBlock->NClusters() : 0;
    for(int i = 0; i < ntrk_de; i++) {
      TStnTrack* trk = Data->fTrackBlockDe->Track(i);
      //    printf(" -- ttrack t0: %10.3f\n",trk->T0());
      for(int j = 0; j < ncalo_cl; j++) {
        TStnCluster* cl = Data->fClusterBlock->Cluster(j);
        float dt = trk->T0() - cl->Time();

        // printf(" -- cluster E, T0, dt: %10.3f %10.3f %10.3f \n",cl->Energy(), cl->Time(),dt);
        if(cl->Energy() > 140.) {
          if(fabs(dt) < 15.) {
            //-----------------------------------------------------------------------------
            // suspect misreconstructed upstream leg, so the cluster is in-time
            //-----------------------------------------------------------------------------
            VetoID |= ConvAna::kCaloInTimeVetoBit;
          } else if((dt > 50) && (dt <= 200)) {
            //-----------------------------------------------------------------------------
            // suspect missed upstream leg, and only the downstream leg reconstructed , so the
            // cluster is early
            //-----------------------------------------------------------------------------
            VetoID |= ConvAna::kCaloEarlyVetoBit;
          }
        }
      }
    }

    return VetoID;
  }

  //-----------------------------------------------------------------------------
  bool TAnaModule::IsSignal(int sim_code) {
    switch(sim_code) {
    case mu2e::ProcessCode::mu2eCeMinusEndpoint:
    case mu2e::ProcessCode::mu2eCePlusEndpoint:
    case mu2e::ProcessCode::mu2eCeMinusLeadingLog:
    case mu2e::ProcessCode::mu2eCePlusLeadingLog:
      return true;
    }
    return false;
  }

  //-----------------------------------------------------------------------------
  bool TAnaModule::IsBackground(int sim_code) {
    switch(sim_code) {
    case mu2e::ProcessCode::mu2eFlateMinus:
    case mu2e::ProcessCode::mu2eFlatePlus:
    case mu2e::ProcessCode::mu2eFlatPhoton:
    case mu2e::ProcessCode::mu2eAntiproton:
    case mu2e::ProcessCode::mu2eDIOLeadingLog:
    case mu2e::ProcessCode::mu2eMuonDecayAtRest: // FIXME: Ensure this isn't a frequent code
    case mu2e::ProcessCode::mu2eGammaConversion:
    case mu2e::ProcessCode::mu2eExternalRMC:
    case mu2e::ProcessCode::mu2eInternalRMC:
    case mu2e::ProcessCode::mu2eExternalRPC:
    case mu2e::ProcessCode::mu2eInternalRPC:
      return true;
    }
    // FIXME: Add cosmics check
    // if(sim_code == mu2e::ProcessCode::mu2ePrimary && abs(pdg) == 13) return true
    return false;
  }

  //-----------------------------------------------------------------------------
  bool TAnaModule::IsBeamProcess(int sim_code) {
    switch(sim_code) {
    case mu2e::ProcessCode::mu2eFlateMinus:
    case mu2e::ProcessCode::mu2eFlatePlus:
    case mu2e::ProcessCode::mu2eFlatPhoton:
    case mu2e::ProcessCode::mu2eAntiproton:
    case mu2e::ProcessCode::DIO:
    case mu2e::ProcessCode::mu2eMuonDecayAtRest:
    case mu2e::ProcessCode::mu2eDIOLeadingLog:
    case mu2e::ProcessCode::mu2eGammaConversion:
    case mu2e::ProcessCode::mu2eExternalRMC:
    case mu2e::ProcessCode::mu2eInternalRMC:
    case mu2e::ProcessCode::mu2eExternalRPC:
    case mu2e::ProcessCode::mu2eInternalRPC:
    case mu2e::ProcessCode::mu2eCeMinusEndpoint:
    case mu2e::ProcessCode::mu2eCePlusEndpoint:
    case mu2e::ProcessCode::mu2eCeMinusLeadingLog:
    case mu2e::ProcessCode::mu2eCePlusLeadingLog:
      return true;
    }
    return false;
  }

  //-----------------------------------------------------------------------------
  // Generated intensity distribution PDF
  double TAnaModule::BatchModeWeight(float lumi, int mode) {
    if(mode <= 0)
      return 1.;
    if(mode > 2)
      return 1.;
    //-----------------------------------------------------------------------------
    // Batch mode 1/2 log normal initialization
    //-----------------------------------------------------------------------------
    const static double mean_b1 = 1.6e7;
    const static double mean_b2 = 3.9e7;
    const static double sigma = 0.7147;
    const static double mub1 = log(mean_b1) - 0.5 * sigma * sigma;
    const static double mub2 = log(mean_b2) - 0.5 * sigma * sigma;
    const static double cut_off_norm_b1 = ROOT::Math::lognormal_cdf(1.2e8, mub1, sigma); // Due to max cutoff in generation
    const static double cut_off_norm_b2 = ROOT::Math::lognormal_cdf(1.2e8, mub2, sigma); // Due to max cutoff in generation
    if(mode == 1) {
      const double p1 = ROOT::Math::lognormal_pdf(lumi, mub1, sigma) / cut_off_norm_b1;
      return p1;
    }
    const double p2 = ROOT::Math::lognormal_pdf(lumi, mub2, sigma) / cut_off_norm_b2;
    return p2;
  }

  //-----------------------------------------------------------------------------
  // Re-weight a beam-related process for a given intensity mode
  double TAnaModule::BeamProcessWeight(float lumi, int mode) {
    if(mode <= 0)
      return 1.;
    if(mode > 2)
      return 1.;
    //-----------------------------------------------------------------------------
    // Batch mode 1/2 log normal initialization
    //-----------------------------------------------------------------------------
    const static double mean_b1 = 1.6e7;
    const static double mean_b2 = 3.9e7;
    const static double sigma = 0.7147;
    const static double mub1 = log(mean_b1) - 0.5 * sigma * sigma;
    const static double mub2 = log(mean_b2) - 0.5 * sigma * sigma;
    const static double xlognorm_norm_b1 = 6.293492e-8;                                  // evaluated in ROOT for cut_off = 1.2e8, above mub(1/2) and sigma for
                                                                                         // x*log-normal(x)
    const static double xlognorm_norm_b2 = 2.887949e-8;                                  // evaluated in ROOT
    const static double cut_off_norm_b1 = ROOT::Math::lognormal_cdf(1.2e8, mub1, sigma); // Due to max cutoff in generation
    const static double cut_off_norm_b2 = ROOT::Math::lognormal_cdf(1.2e8, mub2, sigma); // Due to max cutoff in generation
    if(mode == 1) {
      const double p1 = lumi * xlognorm_norm_b1 * cut_off_norm_b1;
      return p1;
    }
    const double p2 = lumi * xlognorm_norm_b2 * cut_off_norm_b2;
    return p2;
  }

  float TAnaModule::RMCWeight(const float gen_energy, const int spectrum, const float kmax) {
    float weight = 1.f;
    constexpr float ref_energy = 57.f;
    if(spectrum == kClosure) { // Closure approximation
      weight = TStntuple::RMC_ClosureAppxWeight(gen_energy, kmax);
      // Normalize the spectrum between 57 MeV - inf
      weight /= TStntuple::RMC_ClosureAppxIntegral(ref_energy, kmax, kmax);
    } else if(spectrum == kPlestid) { // Phase-space approximation
      weight = 0.f;
      // Fit results to the 1995 TRIUMF Al spectrum
      constexpr double br_0n         = 0.099; // defined for E_photon > 57 MeV
      constexpr double br_1n         = 0.901;
      constexpr double br_2n         = 0.000;
      constexpr double kmax_0n       = 101.866;
      constexpr double kmax_1n       =  95.449;
      constexpr double kmax_2n       =  84.395;
      constexpr double br_ref_energy = 57.; // Get the branching fractions for the full space
      const     double br_0n_full    = (br_0n <= 0.) ? 0. : br_0n / TStntuple::RMC_PlestidIntegral(br_ref_energy, kmax_0n, kmax_0n, 0);
      const     double br_1n_full    = (br_1n <= 0.) ? 0. : br_1n / TStntuple::RMC_PlestidIntegral(br_ref_energy, kmax_1n, kmax_1n, 1);
      const     double br_2n_full    = (br_2n <= 0.) ? 0. : br_2n / TStntuple::RMC_PlestidIntegral(br_ref_energy, kmax_2n, kmax_2n, 2);
      weight += br_0n_full * TStntuple::RMC_PlestidWeight(gen_energy, kmax_0n, 0);
      weight += br_1n_full * TStntuple::RMC_PlestidWeight(gen_energy, kmax_1n, 1);
      weight += br_2n_full * TStntuple::RMC_PlestidWeight(gen_energy, kmax_2n, 2);

      // Normalize the spectrum between 57 MeV - inf
      double integral = 0.;
      integral += br_0n_full * TStntuple::RMC_PlestidIntegral(ref_energy, kmax_0n, kmax_0n, 0);
      integral += br_1n_full * TStntuple::RMC_PlestidIntegral(ref_energy, kmax_1n, kmax_1n, 1);
      integral += br_2n_full * TStntuple::RMC_PlestidIntegral(ref_energy, kmax_2n, kmax_2n, 2);
      weight /= integral;
    } else {
      throw std::runtime_error("Unknown RMC spectrum option!\n");
    }
    return weight;
  }


  //_____________________________________________________________________________
  int TAnaModule::EndJob() {
    printf("----- end job: ---- %s\n", GetName());

    return 0;
  }
} // namespace ConvAna
