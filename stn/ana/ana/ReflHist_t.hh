#ifndef __mumep_ana_ana_ReflHist_t_hh
#define __mumep_ana_ana_ReflHist_t_hh

#include "TH1.h"
#include "TH2.h"

namespace mumep_ana {

  struct ReflHist_t {
    TH1F* fPUpstream;
    TH1F* fPDownstream;
    TH1F* fDeltaT0;
    TH1F* fDeltaTFront;
    TH1F* fDeltaP;
    TH1F* fEpDownstream;
    TH1F* fPdgUpstream;
    TH1F* fPdgDownstream;
  };

} // namespace mumep_ana
#endif
