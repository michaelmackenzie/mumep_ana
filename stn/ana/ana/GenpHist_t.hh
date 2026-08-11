#ifndef __mumep_ana_ana_GenpHist_t_hh
#define __mumep_ana_ana_GenpHist_t_hh

#include "TH1.h"

namespace mumep_ana {

  struct GenpHist_t {
    TH1F* fPdgCode[2]; // same distribution in different scale
    TH1F* fGenID;      //
    TH1F* fZ0;         //
    TH1F* fT0;         //
    TH1F* fR0;         //
    TH1F* fP;          //
    TH1F* fCosTh;      //
  };
} // namespace mumep_ana
#endif
