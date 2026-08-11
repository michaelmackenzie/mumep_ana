// Structure to facilitate passing multiple data blocks for CRV veto checks
#ifndef __mumep_ana_ana_CosmicVetoData_t_hh
#define __mumep_ana_ana_CosmicVetoData_t_hh

#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"

namespace mumep_ana {

  enum {
    kNTrkDeVetoBit = 0x0001,     // N(De tracks) > 1
    kNTrkUeVetoBit = 0x0002,     // N(Ue tracks) > 1
    kTrkUeDtVetoBit = 0x0004,    // Dt(De - Ue) > 50 ns
    kTrkTcVetoBit = 0x0008,      // Dt(De - Ue) > 50 ns
    kNHelDeVetoBit = 0x0010,     // N(De helices) > 1
    kNHelUeVetoBit = 0x0020,     // N(Ue helices) > 1
    kCaloInTimeVetoBit = 0x0040, // too energetic cluster in the calorimeter
    kCaloEarlyVetoBit = 0x0080,  // too energetic cluster in the calorimeter
    kCrvStubVetoBit = 0x0100     // CRV stub close in time
  };

  struct CosmicVetoData_t {
    TStnTimeClusterBlock* fTCFinderBlockUe = nullptr; // don't immediately need the De one
    TStnHelixBlock* fHelixBlockDe = nullptr;
    TStnHelixBlock* fHelixBlockUe = nullptr;
    TStnTrackBlock* fTrackBlockDe = nullptr;
    TStnTrackBlock* fTrackBlockUe = nullptr;
    TStnClusterBlock* fClusterBlock = nullptr;
  };

}; // namespace mumep_ana

#endif
