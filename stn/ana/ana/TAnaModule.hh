///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __ConvAna_ana_TAnaModule_hh__
#define __ConvAna_ana_TAnaModule_hh__

// standard includes
#include <climits>

// ROOT includes
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

// CLHEP includes
#include "CLHEP/Units/PhysicalConstants.h"

// mu2e includes
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/Mu2eUtilities/inc/StopWatch.hh"

// Stntuple includes
#include "Stntuple/alg/RMCSpectra.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TCrvCoincidenceCluster.hh"
#include "Stntuple/obj/TGenParticle.hh"
#include "Stntuple/obj/TSimParticle.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnTrack.hh"
#include "Stntuple/obj/TStnTriggerBlock.hh"

// local includes
#include "ConvAna/ana/CRVHist_t.hh"
#include "ConvAna/ana/CRVStubPar_t.hh"
#include "ConvAna/ana/ClusterHist_t.hh"
#include "ConvAna/ana/ClusterPar_t.hh"
#include "ConvAna/ana/CosmicVetoData_t.hh"
#include "ConvAna/ana/EventHist_t.hh"
#include "ConvAna/ana/EventPar_t.hh"
#include "ConvAna/ana/GenpHist_t.hh"
#include "ConvAna/ana/HelixHist_t.hh"
#include "ConvAna/ana/HelixPar_t.hh"
#include "ConvAna/ana/SimpHist_t.hh"
#include "ConvAna/ana/TrackHist_t.hh"
#include "ConvAna/ana/TrackPar_t.hh"
#include "ConvAna/ana/TriggerInfo.hh"

namespace ConvAna {
  class TAnaModule : public TStnModule {
  public:
    //-----------------------------------------------------------------------------
    //  define constexpr
    //-----------------------------------------------------------------------------
    static constexpr float mmTconversion = CLHEP::c_light / 1000.0;
    static constexpr float bz0 = 1.0;

    // track selection bits
    enum { kP = 0, kRMax = 1, kTrkQual = 2, kT0 = 3, kFitCon = 4, kClusterE = 5, kD0 = 6, kTDip = 7, kT0Loose = 8, kUpstream = 10, kPID = 11, kFitHyp = 12, kCosmicID = 13, kCRV = 15, kMC = 20 };
    std::map<int, TString> kTrackIDNames = {{kP, "P"},           {kRMax, "R(max)"},           {kTrkQual, "TrkQual"},   {kT0, "T_{0}"}, {kFitCon, "p(#chi^2)"}, {kClusterE, "E_{CL}"},    {kD0, "D_{0}"},
                                            {kTDip, "tan(dip)"}, {kT0Loose, "T_{0} (loose)"}, {kUpstream, "Upstream"}, {kPID, "PID"},  {kFitHyp, "Fit hyp."},  {kCosmicID, "Cosmic ID"}, {kCRV, "CRV"},
                                            {kMC, "MC"}};

    // PDG -> name
    static TString NameFromPDG(const int pdg) {
      const int abs_pdg = std::abs(pdg);
      if(abs_pdg == 11)
        return "Electron";
      if(abs_pdg == 13)
        return "Muon";
      if(abs_pdg == 22)
        return "Photon";
      if(pdg == -2212)
        return "Antiproton";
      if(abs_pdg == 2212)
        return "Proton";
      if(abs_pdg == 2112)
        return "Neutron";
      if(abs_pdg == 211)
        return "Pion";
      if(abs_pdg == 111)
        return "Pi-0";
      if(abs_pdg > 10000)
        return "Ion";
      return "Unknown";
    }

    //-----------------------------------------------------------------------------
    //  data members
    //-----------------------------------------------------------------------------
  public:
    // normalization info
    struct NormInfo_t {
      Long64_t nseen_ = 0;
    };

    TStnTriggerBlock* fTriggerBlock = nullptr;
    TString fTriggerBlockName;

    ConvAna::EventPar_t fEvtPar;
    TTree* fNormTree;
    NormInfo_t fNormInfo;
    TStntuple* fStntuple; // STNTUPLE singleton, for algorithm access
    int fDebugLevel = 0;
    TSimParticle* fPrimary = nullptr; // Physics primary

    //-----------------------------------------------------------------------------
    //  Define the RMC spectrum information to be used
    //-----------------------------------------------------------------------------
    enum { kClosure, kClosureFlat, kClosureExp, kClosureTransition, kPlestid, kLAST };
    int fSpectrum;
    int fIntSpectrum;
    /** Spectrum parameters **/
    // Closure: Not used
    // Closure+Flat: 1 = Br(flat > kmax) / Br(RMC) 2 = E(max) of flat
    // Closure+Exp: 1 = intersection point of exp and Closure 2 = slope of exp
    // Closure+Transition: 1 = Br(transition) / Br(RMC) 2 = Transition energy +- 0.1 MeV
    float fSpectrumParam[2];
    float fKMax;
    float fKinematicLimit;
    RMCSpectra* fRMCSpectra; // information on RMC spectra

    std::vector<int> fTriggersPassed; // List of trigger bits that passed
    std::map<int, int> fTrigBitToBin;
    int fMakeTrigHists = 0;

    mu2e::StopWatch* fWatch; // for timing

    //-----------------------------------------------------------------------------
    //  functions
    //-----------------------------------------------------------------------------
  public:
    TAnaModule(const char* name = "ConvAna_Ana", const char* title = "Ana");
    ~TAnaModule();

    //-----------------------------------------------------------------------------
    // overloaded methods of TStnModule
    //-----------------------------------------------------------------------------
    virtual int BeginJob();
    virtual int BeginRun();
    virtual int EndJob();

    //-----------------------------------------------------------------------------
    // other methods
    //-----------------------------------------------------------------------------
    void BookNormTree();
    void BookEventHistograms(ConvAna::EventHist_t* Hist, const char* Folder);
    void BookGenpHistograms(ConvAna::GenpHist_t* Hist, const char* Folder);
    void BookSimpHistograms(ConvAna::SimpHist_t* Hist, const char* Folder);
    void BookTrackHistograms(ConvAna::TrackHist_t* Hist, const char* Folder);
    void BookHelixHistograms(ConvAna::HelixHist_t* Hist, const char* Folder);
    void BookClusterHistograms(ConvAna::ClusterHist_t* Hist, const char* Folder);
    void BookCRVClusterHistograms(ConvAna::CRVClusterHist_t* Hist, const char* Folder);

    void FillEventHistograms(ConvAna::EventHist_t* Hist, ConvAna::EventPar_t* EvtPar, float Weight = 1.f);
    void FillGenpHistograms(ConvAna::GenpHist_t* Hist, TGenParticle* Genp, float Weight = 1.f);
    void FillSimpHistograms(ConvAna::SimpHist_t* Hist, TSimParticle* Simp, float Weight = 1.f);
    void FillTrackHistograms(ConvAna::TrackHist_t* Hist, ConvAna::TrackPar_t* TrkPar, float Weight = 1.f);
    void FillHelixHistograms(ConvAna::HelixHist_t* Hist, ConvAna::HelixPar_t* HlxPar, float Weight = 1.f);
    void FillClusterHistograms(ConvAna::ClusterHist_t* Hist, ConvAna::ClusterPar_t* Par, float Weight = 1.f);
    void FillCRVClusterHistograms(ConvAna::CRVClusterHist_t* Hist, ConvAna::CRVStubPar_t* CrvPar, float Weight = 1.f);

    void InitEventInfo();
    float RMCWeight(const float gen_energy, const int spectrum, const float kmax = 90.1f);

    void InitTrackPar(TStnTrack* Trk, ConvAna::TrackPar_t* TrkPar, TStnHelix* Hlx = nullptr);
    void InitHelixPar(TStnHelix* Hlx, ConvAna::HelixPar_t* HlxPar);
    void InitClusterPar(TStnCluster* Cluster, ConvAna::ClusterPar_t* Par);
    void InitCRVStubPar(TCrvClusterBlock* CrvClusterBlock, ConvAna::CRVStubPar_t* CrvStubPar, int maxStubs = INT_MAX, TSimParticle* Simp = nullptr);
    CRVStubPar_t* MatchCRVToTrack(TStnTrack* track, TCrvClusterBlock* clusters, CRVStubPar_t* stubPars, int maxStubs = INT_MAX);
    TStnTrack* MatchUpstreamTrack(TStnTrack* track, TStnTrackBlock* tracks);
    bool ResolveAmbiguity(TStnTrack* t_1, TStnTrack* t_2) {
      if(!t_2)
        return true;
      if(!t_1)
        return false;
      // if only 1 has a cluster, prefer that one
      if(t_1->fClusterE > 0. && t_2->fClusterE <= 0.)
        return true;
      if(t_2->fClusterE > 0. && t_1->fClusterE <= 0.)
        return false;
      // use the fit quality to decide
      return t_1->fFitCons >= t_2->fFitCons;
    }

    static TString ProcessGroup(int code) {
      switch(code) {
      case mu2e::ProcessCode::mu2eAntiproton:
        return "Antiproton";
      case mu2e::ProcessCode::mu2eFlateMinus:
      case mu2e::ProcessCode::DIO:
      case mu2e::ProcessCode::mu2eMuonDecayAtRest:
      case mu2e::ProcessCode::mu2eDIOLeadingLog:
        return "DIO";
      case mu2e::ProcessCode::mu2eFlatPhoton:
      case mu2e::ProcessCode::mu2eGammaConversion:
      case mu2e::ProcessCode::mu2eExternalRMC:
        return "Ext-RMC";
      case mu2e::ProcessCode::mu2eFlatePlus:
      case mu2e::ProcessCode::mu2eInternalRMC:
        return "Int-RMC";
      case mu2e::ProcessCode::mu2eExternalRPC:
        return "Ext-RPC";
      case mu2e::ProcessCode::mu2eInternalRPC:
        return "Int-RPC";
      case mu2e::ProcessCode::mu2eCeMinusEndpoint:
      case mu2e::ProcessCode::mu2eCeMinusLeadingLog:
        return "CE-";
      case mu2e::ProcessCode::mu2eCePlusEndpoint:
      case mu2e::ProcessCode::mu2eCePlusLeadingLog:
        return "CE+";
      case mu2e::ProcessCode::mu2ePrimary:
        return "Primary"; // This should typically be cosmics
      case 301:
        return "IPA DIO"; // Created codes
      }
      return Form("Unknown-%i", code);
    }

    static const bool IsIPADIO(TSimParticle* sim) {
      if(!sim)
        return false;
      if(sim->fPdgCode != 11)
        return false;
      const float x(sim->fStartPos.X() + 3904.), y(sim->fStartPos.Y()), z(sim->fStartPos.Z());
      const float r(std::sqrt(x * x + y * y));
      return std::fabs(r - 300.) < 10. && z > 6800. && z < 8000.;
    }

    bool IsSignal(int sim_code);
    bool IsBackground(int sim_code);
    bool IsBeamProcess(int sim_code);
    double BatchModeWeight(float lumi, int mode);
    double BeamProcessWeight(float lumi, int mode);
    double PBarWeight(float z, float time, float r) {
      // FIXME: Initial corrections
      const float z_wt = (z < 5520.) ? 0.689113 : 57.4314;
      const float t_wt = (time < 850.) ? 5.49869 * std::exp(-0.5 * std::pow((time - 636.651) / 98.3452, 2)) : std::max(0., 0.906741 - 0.00422318 * time);
      const float r_wt = (r < 25.f) ? std::exp(6.5611 - 0.279425 * r) : (r < 65.f) ? 0.472027 + 0.00398991 * r : (r < 75.f) ? std::exp(-15.0391 + 0.223124 * r) : 0.f;
      return z_wt * t_wt * r_wt;
    }

    int NonCRVCosmicVeto(CosmicVetoData_t* Data);

    // initialize the list of triggers that passed
    void SetTriggersPassed() {
      fTriggersPassed.clear();
      if(!fTriggerBlock)
        return;
      for(int bit = 100; bit < 1000; ++bit) {
        if(fTriggerBlock->PathPassed(bit))
          fTriggersPassed.push_back(bit);
      }
    }

    // Set trigger overlap to N(X && Y) / N(X)
    void NormalizeTriggerOverlap(EventHist_t* Hist) {
      if(!Hist)
        return;
      TH2* h = Hist->fTrigOverlap[0];
      if(!h)
        return;
      for(int xbin = 1; xbin <= h->GetNbinsX(); ++xbin) {
        const char* xbin_label = h->GetXaxis()->GetBinLabel(xbin);
        const int ybin_norm = h->GetYaxis()->FindBin(xbin_label);
        const double norm = h->GetBinContent(xbin, ybin_norm);
        if(norm <= 0.)
          continue;
        for(int ybin = 1; ybin <= h->GetNbinsY(); ++ybin) {
          const double value = h->GetBinContent(xbin, ybin);
          if(value <= 0.)
            continue;
          h->SetBinContent(xbin, ybin, value / norm);
          h->SetBinError(xbin, ybin, std::sqrt(value) / norm);
        }
      }
    }

    // Print a cut-flow result
    void PrintCutFlow(TrackHist_t* Hist, const char* name = nullptr) {
      if(!Hist) return;
      TH1* h_i = Hist->fTrackID   ; // inclusive
      TH1* h_e = Hist->fExlTrackID; // exclusive
      TH1* h_p = Hist->fP[0]; // momentum for normalization
      if(!h_i || !h_e || !h_p) return;
      const double norm = h_p->Integral(0, h_p->GetNbinsX()+1);
      if(norm <= 0.) return;
      printf("[TAnaModule::%s] Cut-Flow information", __func__);
      if(name) printf(" for set %s", name);
      printf(":\n");
      printf("%15s: %10s %10s %10s %10s\n",
             "Selection", "Inclusive", "Exclusive",
             "Frac(I)", "Frac(E)");
      for(int icut = -1; icut <= kMC; ++icut) {
        if(icut >= 0 && !kTrackIDNames.contains(icut)) continue; // cut isn't defined
        TString cut_name = (icut == -1) ? "Passed" : kTrackIDNames.at(icut);
        const int bin = h_i->GetXaxis()->FindBin(cut_name.Data());
        double n_i = h_i->GetBinContent(bin);
        double n_e = h_e->GetBinContent(bin);
        if(icut >= 0) { // re-write as accepted counts
          n_i = norm - n_i;
          n_e = norm - n_e;
        }
        printf("%15s: %10.4g %10.4g %10.4g %10.4g\n",
               cut_name.Data(), n_i, n_e, n_i/norm, n_e/norm);
      }
    }

    void Debug();

    ClassDef(TAnaModule, 0)
  };
} // namespace ConvAna
#endif
