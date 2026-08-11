#ifndef __mumep_ana_ana_PBarPar_t__
#define __mumep_ana_ana_PBarPar_t__

#include "Stntuple/obj/TSimParticle.hh"
#include "Stntuple/obj/TStnHelix.hh"
#include "Stntuple/obj/TStnTrack.hh"

namespace mumep_ana {

  struct PBarPar_t {

    TStnTrack* fLeadTrk = nullptr;
    TStnTrack* fTrailTrk = nullptr;
    TStnHelix* fLeadHlx = nullptr;
    TStnHelix* fTrailHlx = nullptr;

    void reset() {
      fLeadTrk = nullptr;
      fTrailTrk = nullptr;
      fLeadHlx = nullptr;
      fTrailHlx = nullptr;
    }
  };
} // namespace mumep_ana
#endif
