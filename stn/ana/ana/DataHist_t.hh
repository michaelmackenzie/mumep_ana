#ifndef __mumep_ana_ana_DataHist_t_hh
#define __mumep_ana_ana_DataHist_t_hh

#include "TH1.h"

namespace mumep_ana {

  struct EventHist_t {
    TH1F* fInstLumi;
    TH1F* fInstLumiApr;    // only filled if apr triggered event
    TH1F* fInstLumiCpr;    // only filled if cpr triggered event
    TH1F* fInstLumiAprCpr; // filled if either apr or cpr triggered event
    TH1F* fEventWeight[2];
    TH1F* fNAprHelices;
    TH1F* fNCprHelices;
    TH1F* fNHelices;
    TH1F* fNMatchedHelices;
    TH1F* fNAprTracks;
    TH1F* fNCprTracks;
    TH1F* fNTracks;
  };

} // namespace mumep_ana
#endif
