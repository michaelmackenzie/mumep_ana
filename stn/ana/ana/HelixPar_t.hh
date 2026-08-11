#ifndef __mumep_ana_ana_HelixPar_t__
#define __mumep_ana_ana_HelixPar_t__

#include "Stntuple/obj/TStnHelix.hh"

namespace mumep_ana {

  struct HelixPar_t {

    TStnHelix* fHelix;

    float fRMax;

    // truth-level info
    bool fIsMCDownstream;
    float fTZSigMC; // TZ slope significance, signed by the MC true particle direction
  };
} // namespace mumep_ana
#endif
