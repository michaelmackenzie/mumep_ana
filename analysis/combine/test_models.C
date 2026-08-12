// Simple script to test model retrieval

#include "background_model.C"
#include "signal_model.C"

int test_models() {
  auto cosmic_model = read_model("cosmic", "mumem", 20);
  auto dio_model    = read_model("dio"   , "mumem", 20);
  auto signal_model = read_model("signal", "mumem", 20);
  if(!cosmic_model.pdf_ || !dio_model.pdf_ || !signal_model.pdf_ || !signal_model.obs_) return 1;

  auto cosmic_pdf = cosmic_model.pdf_;
  auto dio_pdf    = dio_model   .pdf_;
  auto signal_pdf = signal_model.pdf_;

  RooRealVar* obs = cosmic_model.obs_;

  auto frame = obs->frame();
  cosmic_pdf->plotOn(frame, RooFit::Normalization(cosmic_model.rate_), RooFit::LineColor(cosmic_model.color_));
  dio_pdf   ->plotOn(frame, RooFit::Normalization(dio_model   .rate_), RooFit::LineColor(dio_model   .color_));
  signal_pdf->plotOn(frame, RooFit::Normalization(signal_model.rate_), RooFit::LineColor(signal_model.color_));

  TCanvas* c = new TCanvas();
  frame->Draw();
  return 0;
}
