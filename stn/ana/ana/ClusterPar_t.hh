#ifndef __mumep_ana_ana_ClusterPar_t__
#define __mumep_ana_ana_ClusterPar_t__

#include "Stntuple/obj/TSimParticle.hh"
#include "Stntuple/obj/TStnCluster.hh"

namespace mumep_ana {

  struct ClusterPar_t {

    TStnCluster* fCluster;
    TSimParticle* fSim;
    int core_crystals[100];

    float gen_energy;

    ClusterPar_t() { init(); }

    void init(TStnCluster* cluster = nullptr) {
      fCluster = cluster;
      fSim = nullptr;
      gen_energy = 0.f;
      for(int i = 0; i < 100; ++i)
        core_crystals[i] = -1;
      set_core_crystals();
    }

    void set_core_crystals() {
      if(!fCluster)
        return;
      const int ncr = fCluster->NCrystals();
      if(ncr == 0)
        return;

      // Starting from the core, add crystals connected to the core until none are added
      core_crystals[0] = 0;
      int n_added = 1;
      bool added = true;
      while(added) {
        added = false;
        // Check each crystal
        for(int icr = 1; icr < ncr; ++icr) {
          if(n_added >= 100)
            break; // prevent overflow of core_crystals[100]
          const float x0 = fCluster->CrystalX(icr);
          const float y0 = fCluster->CrystalY(icr);
          bool already_added = false;
          for(int index = 0; index < n_added; ++index) {
            const int nbr = core_crystals[index];
            if(nbr == icr) {
              already_added = true;
              break; // this crystal is already added
            }
          }
          if(already_added)
            continue;
          for(int index = 0; index < n_added; ++index) {
            const int nbr = core_crystals[index];
            const float x1 = fCluster->CrystalX(nbr);
            const float y1 = fCluster->CrystalY(nbr);
            const double r = std::sqrt(std::pow(x0 - x1, 2) + std::pow(y0 - y1, 2));
            if(r < 50.) { // neighbor
              added = true;
              core_crystals[n_added] = icr;
              ++n_added;
              break;
            }
          }
        }
        if(n_added >= 100)
          break; // stop if array is full
      }
    }

    float core_energy() {
      if(!fCluster)
        return 0.f;
      if(core_crystals[0] < 0)
        set_core_crystals();
      float energy = 0.f;
      int index = 0;
      while(index < 100 && core_crystals[index] >= 0) {
        energy += fCluster->CrystalE(core_crystals[index]);
        ++index;
      }
      return energy;
    }

    int n_core_crystals() {
      if(!fCluster)
        return 0;
      if(core_crystals[0] < 0)
        set_core_crystals();
      int index = 0;
      while(index < 100 && core_crystals[index] >= 0) {
        ++index;
      }
      return index;
    }
  };
} // namespace mumep_ana
#endif
