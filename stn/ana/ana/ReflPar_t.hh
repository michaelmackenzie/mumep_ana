// Information about a reflection candidate
#ifndef __mumep_ana_ana_ReflPar_t__
#define __mumep_ana_ana_ReflPar_t__

#include "Stntuple/obj/TSimParticle.hh"
#include "Stntuple/obj/TStnHelix.hh"
#include "Stntuple/obj/TStnTrack.hh"

namespace mumep_ana {

  struct ReflPar_t {

    TStnTrack* fUpstreamTrk = nullptr;
    TStnTrack* fDownstreamTrk = nullptr;
    TStnHelix* fUpstreamHlx = nullptr;
    TStnHelix* fDownstreamHlx = nullptr;
    TStnTrack* fUpstreamTrigTrk = nullptr;
    TStnTrack* fDownstreamTrigTrk = nullptr;

    void reset() {
      fUpstreamTrk = nullptr;
      fDownstreamTrk = nullptr;
      fUpstreamHlx = nullptr;
      fDownstreamHlx = nullptr;
      fUpstreamTrigTrk = nullptr;
      fDownstreamTrigTrk = nullptr;
    }

    void set(TStnTrack* Trk_u, TStnTrack* Trk_d) {
      reset();
      fUpstreamTrk = Trk_u;
      fDownstreamTrk = Trk_d;
    }

    float deltaT0() { return (fUpstreamTrk && fDownstreamTrk) ? fDownstreamTrk->fT0 - fUpstreamTrk->fT0 : -1.e10; }
    float deltaTFront() { return (fUpstreamTrk && fDownstreamTrk) ? fDownstreamTrk->fTFront - fUpstreamTrk->fTFront : -1.e10; }
    float deltaP() { return (fUpstreamTrk && fDownstreamTrk) ? fDownstreamTrk->P() - fUpstreamTrk->P() : -1.e10; }
  };
} // namespace mumep_ana
#endif
