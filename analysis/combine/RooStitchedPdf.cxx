#include "RooStitchedPdf.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include <cmath>

// Constructors
RooStitchedPdf::RooStitchedPdf() {}

RooStitchedPdf::RooStitchedPdf(const char *name, const char *title,
                               RooAbsReal& _obs,
                               RooAbsPdf&  _pdf1,
                               RooAbsPdf&  _pdf2,
                               RooAbsReal& _cutVal) :
  RooAbsPdf(name, title),
  obs("obs", "Observable", this, _obs),
  pdf1("pdf1", "First PDF (Left Subrange)", this, _pdf1),
  pdf2("pdf2", "Second PDF (Right Subrange)", this, _pdf2),
  cutVal("cutVal", "Stitching Boundary Cut Point", this, _cutVal)
{
}

RooStitchedPdf::RooStitchedPdf(const RooStitchedPdf& other, const char* name) :
  RooAbsPdf(other, name),
  obs("obs", this, other.obs),
  pdf1("pdf1", this, other.pdf1),
  pdf2("pdf2", this, other.pdf2),
  cutVal("cutVal", this, other.cutVal)
{
}

Int_t RooStitchedPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const {
  // Check if we are integrating over our primary observable
  if (matchArgs(allVars, analVars, obs)) return 1;
  return 0;
}

Double_t RooStitchedPdf::analyticalIntegral(Int_t code, const char* rangeName) const {
  if (code == 1) {
    double min = obs.min(rangeName);
    double max = obs.max(rangeName);
    double cut = cutVal;
    double integral = 0.0;

    // Get direct pointers to the underlying PDFs from the proxies
    const RooAbsPdf* pdf1_ptr = dynamic_cast<const RooAbsPdf*>(&pdf1.arg());
    const RooAbsPdf* pdf2_ptr = dynamic_cast<const RooAbsPdf*>(&pdf2.arg());

    if (!pdf1_ptr || !pdf2_ptr) {
      return 0.0; // Safety fallback if proxies are invalid
    }

    // Create a RooArgSet of the observable to pass down to sub-PDFs
    RooArgSet obsSet(obs.arg());

    // Define temporary range names to isolate the piecewise integration
    TString rName1 = Form("%s_subRange1", GetName());
    TString rName2 = Form("%s_subRange2", GetName());

    // 1. Integrate pdf1 below the cut boundary
    if (min < cut) {
      double upper_poly = std::min(max, cut);
      // Cast away constness safely on the observable to define a local sub-range
      (const_cast<RooRealVar*>(dynamic_cast<const RooRealVar*>(&obs.arg()))) -> setRange(rName1, min, upper_poly);

      int subCode1 = pdf1_ptr->getAnalyticalIntegral(obsSet, obsSet, rName1);
      // Pass the address (&obsSet) to match the signature
      integral += pdf1_ptr->analyticalIntegralWN(subCode1, &obsSet, rName1);
    }

    // 2. Integrate pdf2 above the cut boundary
    if (max > cut) {
      double lower_exp = std::max(min, cut);
      (const_cast<RooRealVar*>(dynamic_cast<const RooRealVar*>(&obs.arg()))) -> setRange(rName2, lower_exp, max);

      int subCode2 = pdf2_ptr->getAnalyticalIntegral(obsSet, obsSet, rName2);
      // Pass the address (&obsSet) to match the signature
      integral += pdf2_ptr->analyticalIntegralWN(subCode2, &obsSet, rName2);
    }

    return integral;
  }
  return 0;
}

Double_t RooStitchedPdf::evaluate() const
{
  Double_t x = obs;
  Double_t c = cutVal;

  // Region 1: Left of the cut point
  if (x <= c) {
    return pdf1;
  }

  // Region 2: Right of the cut point
  // 1. Fetch pointers to the underlying base arguments
  RooRealVar* obsArg = const_cast<RooRealVar*>(static_cast<const RooRealVar*>(&obs.arg()));
  RooAbsPdf* pdf1Arg = const_cast<RooAbsPdf*>(static_cast<const RooAbsPdf*>(&pdf1.arg()));
  RooAbsPdf* pdf2Arg = const_cast<RooAbsPdf*>(static_cast<const RooAbsPdf*>(&pdf2.arg()));

  // 2. Save current state
  Double_t x_saved = x;

  // 3. Force evaluation at the boundary point 'c'
  obsArg->setVal(c);
  // Setting dirty flags prevents RooFit from serving up a stale cached value
  pdf1Arg->setValueDirty();
  pdf2Arg->setValueDirty();

  Double_t f_c = pdf1Arg->getVal();
  Double_t g_c = pdf2Arg->getVal();

  // 4. Restore state and force-re-evaluate pdf2 at the original 'x' coordinate
  obsArg->setVal(x_saved);
  pdf1Arg->setValueDirty();
  pdf2Arg->setValueDirty();

  Double_t g_x = pdf2Arg->getVal();

  obsArg->setVal(x);

  // 5. If denominator is problematic, fall back safely
  if (g_c <= 0.0) return g_x;
  if (f_c <= 0.0) return g_x;

  // Return smoothly matched scaling factor multiplied by the clean dynamic pdf value
  Double_t scale = f_c / g_c;
  return scale * g_x;
}
