#ifndef ROOSTITCHEDPDF_H
#define ROOSTITCHEDPDF_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooStitchedPdf : public RooAbsPdf {
public:
  RooStitchedPdf();
  RooStitchedPdf(const char *name, const char *title,
                 RooAbsReal& _obs,
                 RooAbsPdf&  _pdf1,
                 RooAbsPdf&  _pdf2,
                 RooAbsReal& _cutVal);
  
  RooStitchedPdf(const RooStitchedPdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const override { return new RooStitchedPdf(*this, newname); }
  inline virtual ~RooStitchedPdf() { }

protected:
  RooRealProxy obs;
  RooRealProxy pdf1;
  RooRealProxy pdf2;
  RooRealProxy cutVal;

  Double_t evaluate() const override;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override;
private:
  ClassDefOverride(RooStitchedPdf, 1) // Piecewise continuous stitched PDF
};

#endif
