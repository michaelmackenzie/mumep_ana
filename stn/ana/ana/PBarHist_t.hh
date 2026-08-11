#ifndef __mumep_ana_ana_PBarHist_t_hh
#define __mumep_ana_ana_PBarHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace mumep_ana {

  struct PBarHist_t {
    TH1F* fPLead;
    TH1F* fPTrail;
    TH1F* fPTotal;
    TH1F* fDeltaT;
    TH1F* fDeltaP;
    TH1F* fHelSum;
  };

} // namespace mumep_ana
#endif
