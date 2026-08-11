#ifndef __mumep_ana_ana_SimpHist_t_hh
#define __mumep_ana_ana_SimpHist_t_hh

#include "TH1.h"

namespace mumep_ana {

  struct SimpHist_t {
    TH1F* fPdgCode;
    TH1F* fCreationCode;
    TH1F* fMomTargetEnd;
    TH1F* fMomTrackerFront;
    TH1F* fMomStart;
    TH1F* fTimeStart;
    TH1F* fTimeEnd;
    TH1F* fStartZ;
    TH1F* fStartR;
    TH2F* fStartXY;
    TH1F* fNStrawHits;
  };
} // namespace mumep_ana
#endif
