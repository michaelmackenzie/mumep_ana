#ifndef __mumep_ana_ana_TrackPar_t__
#define __mumep_ana_ana_TrackPar_t__

#include "TString.h"

#include <map>

#include "mumep_ana/ana/CRVStubPar_t.hh"
#include "Stntuple/obj/TSimParticle.hh"
#include "Stntuple/obj/TStnCluster.hh"
#include "Stntuple/obj/TStnHelix.hh"
#include "Stntuple/obj/TStnTrack.hh"

namespace mumep_ana {

  enum { kTrackIDs = 20, kMaxTrackFits = 10 };

  struct TrackPar_t {

    TStnTrack* fTrack = nullptr;

    float fRMax = 0.f;
    float fRadius = 0.f;
    int fIDWord[kTrackIDs]; // different track selection options
    int fDNhitsUe = 0;
    float fApproxDpST = 0.; // appoximated energy loss estimate from ST exit -> Tracker front based
                            // on IPA intersections

    // tracker hits slope
    float fTZSlope = 0.f;
    float fTZSlopeErr = 0.f;

    // MVAs
    float fTrkQual = -999.f;
    float fPID = -999.f;
    float fTrkPID = -999.f;
    float fCosmicID = -999.f;
    float fOfflinePID = -999.f;

    // Observable
    float fObs = 0.f;

    // Systematic shifts
    float fPUp = 0.f;
    float fPDown = 0.f;

    TStnHelix* fHelix = nullptr;
    TStnCluster* fCluster = nullptr;
    CRVStubPar_t* fCRVStubPar = nullptr; // matched CRV stub
    TStnTrack* fUpstreamTrack = nullptr; // matched upstream leg
    int fNAlt = 0;
    TStnTrack* fAltHypotheses[kMaxTrackFits];

    // MC truth
    float fGenE = 0.f; // energy at generation
    TSimParticle* fSimp = nullptr;
    CRVStubPar_t* fMCCRVStubPar = nullptr; // CRV cluster associated with a true cosmic signal

    //--------------------------------------------
    // functions

    // Reco fit hypothesis
    int FitDirection() const {
      if(!fTrack)
        return 0;
      return (fTrack->fMomentum.Pz() < 0.) ? -1 : 1;
    }

    // Extrapolated intersection with the ST boundary
    int STBoundary() const {
      if(!fTrack)
        return 0;
      if(fTrack->fPSTBack > 0.1)
        return 1;
      if(fTrack->fPSTFront > 0.1)
        return 1;
      return 0;
    }

    // Cos(theta) at the tracker front
    float CosTheta() const {
      if(!fTrack)
        return -999.;
      if(fTrack->P() <= 0.)
        return -2.;
      float cos_theta = fTrack->fMomentum.Pz() / fTrack->P();
      if(cos_theta > 1.f || cos_theta < -1.f) {
        printf("TrackPar::%s: Cos theta out of bounds, pz = %.3f, p = %.3f\n", __func__, fTrack->fMomentum.Pz(), fTrack->P());
        cos_theta = std::min(1.f, std::max(-1.f, cos_theta));
      }
      return cos_theta;
    }

    // Helix-based TZ slope functions
    float TZSlopeSig() const {
      if(fTZSlopeErr > 0.)
        return fTZSlope / fTZSlopeErr;
      return 0.;
    }
    float TZSlopeRatio() const {
      if(!fTrack)
        return 0.;
      const float expected_vel = 300.f * velocity(fTrack->fP, 0.511) * CosTheta(); // FIXME: Using electron mass, mm / ns
      const float ratio = fTZSlope * expected_vel;                                 // dt / dz * (c * dz / dt)
      return ratio;
    }

    // CRV info
    float CRVSTDeltaT() const {
      if(!fTrack || !fCRVStubPar)
        return 1.e5; // no CRV cluster or track
      const float deltat_st = fTrack->fT0 - fCRVStubPar->fApproxTimeSTToFront;
      return deltat_st;
    }
    float CRVCaloDeltaT() const {
      if(!fTrack || !fCRVStubPar)
        return 1.e5; // no CRV cluster or track
      const float deltat_calo = fTrack->fT0 - fCRVStubPar->fApproxTimeCaloToFront;
      return deltat_calo;
    }
    float CRVMinDeltaT() const {
      if(!fTrack || !fCRVStubPar)
        return 1.e5; // no CRV cluster or track
      const float deltat_st = CRVSTDeltaT();
      const float deltat_calo = CRVCaloDeltaT();

      // Minimum of the two hypotheses
      const float min_deltat = (std::fabs(deltat_st) < std::fabs(deltat_calo)) ? deltat_st : deltat_calo;
      return min_deltat;
    }
    float CRVCaloFrontDeltaT(bool muon = false) const { // evaluate the time from calo, rebound, then to tracker front
      if(!fTrack || !fCRVStubPar || fTrack->fP <= 0.)
        return 1.e5;
      const float v = 300.f * velocity(fTrack->fP, (muon) ? 105.66 : 0.511) * 0.5;
      const float z_cal(11820.f), z_st(6271.f), z_trk(8540.f), delta_z(z_cal - z_st);
      const float dt = delta_z / v;
      const float time_rebound = (z_trk - z_st) / v + ((muon) ? 50.f : 35.f); // rough rebound time
      const float crv_time = fCRVStubPar->fApproxTimeCalo + dt + time_rebound;
      return fTrack->fT0 - crv_time;
    }

    // Static functions
    static float velocity(double p, double m) {
      if(p <= 0. && m <= 0.)
        return 0.;
      return p / std::sqrt(p * p + m * m);
    }
  };
} // namespace mumep_ana
#endif
