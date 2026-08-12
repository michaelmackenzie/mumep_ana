// Shared helpers for exporting fitted PDFs/histograms into workspaces
#ifndef __CONVANA_ANALYSIS_FITWORKSPACEUTILS__
#define __CONVANA_ANALYSIS_FITWORKSPACEUTILS__

#include "../tools/utilities.C"

TH1* convolve_with_resolution(TH1* theory, TH1* resolution) {
  TH1* result = (TH1*) theory->Clone("result");
  result->Reset();
  const double estep = 0.01;
  const double e_low  = result->GetXaxis()->GetBinLowEdge(1);
  const double e_high = result->GetXaxis()->GetBinUpEdge(result->GetNbinsX());
  const double de_low  = resolution->GetXaxis()->GetBinLowEdge(1);
  const double de_high = resolution->GetXaxis()->GetBinUpEdge(resolution->GetNbinsX());
  for(double energy = e_low; energy <= e_high; energy += estep) {
    const double p_theory = estep*theory->Interpolate(energy);
    for(double de = de_low; de <= de_high; de += estep) {
      const double p_reco = estep*resolution->Interpolate(de);
      result->Fill(energy + de, p_theory*p_reco);
    }
  }
  if(result->Integral() > 0.) result->Scale(1./result->Integral());
  return result;
}

TH1* convolve_with_resolution(TH1* theory, TF1* resolution) {
  TH1* result = (TH1*) theory->Clone("result");
  result->Reset();
  const double estep = 0.01;
  const double e_low  = result->GetXaxis()->GetBinLowEdge(1);
  const double e_high = result->GetXaxis()->GetBinUpEdge(result->GetNbinsX());
  const double de_low  = resolution->GetXmin();
  const double de_high = resolution->GetXmax();
  for(double energy = e_low; energy <= e_high; energy += estep) {
    const double p_theory = estep*theory->Interpolate(energy);
    for(double de = de_low; de <= de_high; de += estep) {
      const double p_reco = estep*resolution->Eval(de);
      result->Fill(energy + de, p_theory*p_reco);
    }
  }
  if(result->Integral() > 0.) result->Scale(1./result->Integral());
  return result;
}

int save_fit_plot_pair(TCanvas* c,
                       RooPlot* frame,
                       const TString out_base,
                       const double log_ymin = -1.,
                       const double log_ymax = -1.) {
  if(!c || !frame) return 1;

  c->SaveAs(Form("%s.png", out_base.Data()));

  if(log_ymin >= 0. && log_ymax > log_ymin) {
    frame->GetYaxis()->SetRangeUser(log_ymin, log_ymax);
  }

  TPad* pad1 = (TPad*) c->GetPrimitive("pad1");
  if(pad1) pad1->SetLogy();
  c->SaveAs(Form("%s_log.png", out_base.Data()));

  delete frame;
  delete c;
  return 0;
}

int enforce_uniform_if_sparse(TH1* h,
                              RooRealVar& obs,
                              RooAbsPdf*& pdf,
                              const TString pdf_name,
                              const TString pdf_title,
                              const double min_filled_fraction = 0.75) {
  if(!h || !pdf) return -1;

  int filled = 0;
  const int nbins = h->GetNbinsX();
  for(int bin = 1; bin <= nbins; ++bin) {
    if(h->GetBinContent(bin) > 0.) ++filled;
  }

  if(filled < min_filled_fraction*nbins) {
    cout << "Only " << filled << " / " << nbins << " bins have content --> using a uniform model!\n";
    delete pdf;
    pdf = new RooUniform(pdf_name, pdf_title, obs);
  }
  return filled;
}

int save_fit_workspace(TString process,
                       int selection,
                       TString tag,
                       TString component,
                       TString component_title,
                       RooAbsPdf*& pdf,
                       RooRealVar& obs,
                       RooRealVar& norm,
                       const bool hist_pdfs,
                       const TString out_suffix = "") {
  const char* hist_name = Form("%s_%i_%s_hist", process.Data(), selection, component.Data());
  auto h_fit = pdf->createHistogram(hist_name, obs);
  h_fit->SetName(hist_name);

  if(hist_pdfs) {
    pdf->SetName("tmp_pdf");
    RooDataHist* fit_data = new RooDataHist(Form("%s_%i_%s_data_hist", process.Data(), selection, component.Data()),
                                            Form("%s data hist", component_title.Data()),
                                            obs, h_fit);
    pdf = new RooHistPdf(Form("%s_%i_%s_pdf", process.Data(), selection, component.Data()),
                         Form("%s PDF", component_title.Data()),
                         obs, *fit_data);
  }

  const char* fitdir = Form("workspaces/%s%s", process.Data(), (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", fitdir, fitdir));
  TFile* fout = new TFile(Form("%s/%s_fit_%i%s.root", fitdir, component.Data(), selection, out_suffix.Data()), "RECREATE");
  RooWorkspace ws("workspace", "workspace");
  ws.import(*pdf);
  ws.import(norm);
  ws.Write();
  h_fit->Write();
  fout->Close();

  print_pdf(pdf);
  return 0;
}

int save_fit_workspace_with_hist(TString process,
                                 int selection,
                                 TString tag,
                                 TString component,
                                 RooAbsPdf* pdf,
                                 RooRealVar& norm,
                                 TH1* hist,
                                 const TString out_suffix = "") {
  if(!pdf || !hist) return 1;

  const char* fitdir = Form("workspaces/%s%s", process.Data(), (tag == "") ? "" : ("_"+tag).Data());
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", fitdir, fitdir));
  TFile* fout = new TFile(Form("%s/%s_fit_%i%s.root", fitdir, component.Data(), selection, out_suffix.Data()), "RECREATE");
  RooWorkspace ws("workspace", "workspace");
  ws.import(*pdf);
  ws.import(norm);
  ws.Write();
  hist->SetName(Form("%s_%i_%s_hist", process.Data(), selection, component.Data()));
  hist->Write();
  fout->Close();

  print_pdf(pdf);
  return 0;
}

#endif
