#ifndef __mumep_ana_ana_SysHist_t_hh
#define __mumep_ana_ana_SysHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace mumep_ana {

  enum { kMaxSystematics = 100 };
  struct SysHist_t {
    TH1F* fObs[kMaxSystematics];
    TH1F* fDeltaObs[kMaxSystematics];
    TH1F* fWeight[kMaxSystematics];
    TH1F* fDeltaWeight[kMaxSystematics];

    SysHist_t() {
      for(int i = 0; i < kMaxSystematics; ++i) {
        fObs[i] = nullptr;
        fDeltaObs[i] = nullptr;
        fWeight[i] = nullptr;
        fDeltaWeight[i] = nullptr;
      }
    }
  };
} // namespace mumep_ana
#endif
