// Make a copy of a single-channel shape card and its workspace, as written by build_model, with the
// observable restricted to a narrower momentum range. The PDFs are kept as fit over the full
// range: each process's card rate is scaled by the fraction of its PDF inside the new range, and
// Combine renormalizes the PDFs over the narrowed observable. This gives the model of a narrower
// fit range without refitting the components (see sensitivity_vs_fit_range.py).
//
// A freely floating process, i.e. the data-driven background envelope, takes the same treatment:
// Combine multiplies its card rate (1) by its floating <pdf>_norm, so scaling the rate by the
// fraction of the envelope's current function inside the range starts the narrowed model at the
// envelope's yield in that range, while <pdf>_norm still floats freely. Combine renormalizes every
// envelope function over the narrowed observable when it profiles the function index.
//
// Limitations: the range edges are moved to the nearest bin edges of the observable; the shape
// systematic normalization effects, MC statistical uncertainties, and lnN priors are kept as
// evaluated over the full range; the envelope functions keep the shapes fit to the full-range
// sidebands; 2D (p, t0) models are not supported.
#ifndef __CONVANA_ANALYSIS_NARROWMODEL__
#define __CONVANA_ANALYSIS_NARROWMODEL__

#include <fstream>
#include <set>
#include <sstream>

namespace narrow_model_utils {
  std::vector<TString> tokens(const TString& line) {
    std::vector<TString> res;
    std::istringstream stream(line.Data());
    std::string token;
    while(stream >> token) res.push_back(token.c_str());
    return res;
  }

  TString join(const std::vector<TString>& items, const int width = 10) {
    if(items.empty()) return "";
    TString res = Form("%-15s", items[0].Data());
    for(size_t i = 1; i < items.size(); ++i) res += Form(" %-*s", width, items[i].Data());
    return res;
  }

  // nearest bin edge of the binning to x
  double snap_to_edge(const RooAbsBinning& binning, const double x) {
    double best = binning.lowBound();
    for(int i = 0; i <= binning.numBins(); ++i) {
      const double edge = (i < binning.numBins()) ? binning.binLow(i) : binning.highBound();
      if(std::fabs(edge - x) < std::fabs(best - x)) best = edge;
    }
    return best;
  }
}

// card     : input shape card of one channel
// xlow/high: requested momentum range, moved to the nearest bin edges
// out_card : output card; the narrowed workspace is written next to it (.txt --> .root)
int narrow_model(TString card, double xlow, double xhigh, TString out_card) {
  using namespace narrow_model_utils;

  // Read the card
  std::ifstream infile(card.Data());
  if(!infile.is_open()) {
    cout << __func__ << ": Unable to open card " << card.Data() << endl;
    return 1;
  }
  std::vector<TString> lines;
  std::string line_str;
  while(std::getline(infile, line_str)) lines.push_back(line_str.c_str());
  infile.close();

  // Find the shapes, process, and rate lines
  int i_shapes(-1), i_data(-1), i_obs(-1), i_proc(-1), i_rate(-1);
  for(size_t i = 0; i < lines.size(); ++i) {
    const auto t = tokens(lines[i]);
    if(t.empty()) continue;
    if(t[0] == "imax" && t.size() > 1 && t[1] != "1") {
      cout << __func__ << ": Only single-channel cards are supported, not " << card.Data() << endl;
      return 2;
    }
    if(t[0] == "shapes" && t.size() >= 5 && t[1] == "*") i_shapes = i;
    if(t[0] == "shapes" && t.size() >= 5 && t[1] == "data_obs") i_data = i;
    if(t[0] == "observation") i_obs = i;
    if(t[0] == "process" && i_proc < 0) i_proc = i; // the first process line has the names
    if(t[0] == "rate" && i_rate < 0) i_rate = i; // not the "rate group = ..." line
  }
  if(i_shapes < 0 || i_data < 0 || i_obs < 0 || i_proc < 0 || i_rate < 0) {
    cout << __func__ << ": Card " << card.Data() << " is not a shape card written by build_model" << endl;
    return 2;
  }
  const auto shapes = tokens(lines[i_shapes]); // shapes * <obs> <file> workspace:<prefix>$PROCESS_pdf [...]
  const TString obs_name = shapes[2];
  TString ws_file = shapes[3];
  TString pdf_pattern = shapes[4];
  pdf_pattern.ReplaceAll("workspace:", "");
  const auto processes = tokens(lines[i_proc]);
  auto rates = tokens(lines[i_rate]);
  if(processes.size() != rates.size()) {
    cout << __func__ << ": The process and rate lines of " << card.Data() << " differ in length" << endl;
    return 2;
  }

  // The workspace path is relative to the working directory, the card, or the directory above the
  // card (build_model writes datacards/<card> referring to workspaces/<file>)
  if(!ws_file.BeginsWith("/") && gSystem->AccessPathName(ws_file)) {
    const TString card_dir = gSystem->GetDirName(card);
    for(const TString& alt : {card_dir + "/" + ws_file, TString(gSystem->GetDirName(card_dir)) + "/" + ws_file}) {
      if(!gSystem->AccessPathName(alt)) { ws_file = alt; break; }
    }
  }
  TFile* f = TFile::Open(ws_file, "READ");
  if(!f || f->IsZombie()) {
    cout << __func__ << ": Unable to open workspace file " << ws_file.Data() << endl;
    return 3;
  }
  RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
  RooRealVar* obs = (ws) ? ws->var(obs_name) : nullptr;
  RooAbsData* data = (ws) ? ws->data("data_obs") : nullptr;
  if(!ws || !obs || !data) {
    cout << __func__ << ": Workspace, observable " << obs_name.Data() << ", or data_obs missing in " << ws_file.Data() << endl;
    return 3;
  }
  if(data->get()->getSize() != 1) {
    cout << __func__ << ": 2D models are not supported" << endl;
    return 3;
  }

  // Snap the requested range to the observable's bin edges
  const RooAbsBinning& binning = obs->getBinning();
  const double full_min(obs->getMin()), full_max(obs->getMax());
  const double bin_width = (full_max - full_min) / obs->getBins();
  const double lo = snap_to_edge(binning, std::max(xlow , full_min));
  const double hi = snap_to_edge(binning, std::min(xhigh, full_max));
  if(hi - lo < 0.5*bin_width) {
    cout << __func__ << ": Empty range [" << xlow << ", " << xhigh << "] in the observable range ["
         << full_min << ", " << full_max << "]" << endl;
    return 4;
  }
  const int nbins = std::lround((hi - lo) / bin_width);

  // Fraction of each process's PDF inside the new range
  obs->setRange("narrow", lo, hi);
  std::vector<double> fractions;
  for(size_t i = 1; i < processes.size(); ++i) {
    TString pdf_name = pdf_pattern;
    pdf_name.ReplaceAll("$PROCESS", processes[i]);
    RooAbsPdf* pdf = ws->pdf(pdf_name);
    if(!pdf) {
      cout << __func__ << ": PDF " << pdf_name.Data() << " not found in " << ws_file.Data() << endl;
      return 5;
    }
    std::unique_ptr<RooAbsReal> integral(pdf->createIntegral(*obs, RooFit::NormSet(*obs), RooFit::Range("narrow")));
    fractions.push_back(integral->getVal());
  }

  // Data restricted to the new range, binned like the observable
  std::unique_ptr<TH1> h_full(data->createHistogram("narrow_model_data_full", *obs));
  h_full->SetDirectory(nullptr); // owned here, not by the input file
  obs->setRange(lo, hi);
  obs->setBins(nbins);
  TH1D h_data("narrow_model_data", "", nbins, lo, hi);
  h_data.SetDirectory(nullptr);
  for(int bin = 1; bin <= nbins; ++bin) {
    const int full_bin = h_full->FindBin(h_data.GetBinCenter(bin));
    h_data.SetBinContent(bin, h_full->GetBinContent(full_bin));
    h_data.SetBinError  (bin, h_full->GetBinError  (full_bin));
  }
  RooDataHist data_narrow("data_obs_narrow", "Data in the narrowed range", *obs, &h_data);
  ws->import(data_narrow);

  // Write the narrowed workspace next to the output card
  TString out_ws = out_card;
  out_ws.ReplaceAll(".txt", ".root");
  if(out_ws == out_card) out_ws += ".root";
  const TString out_dir = gSystem->GetDirName(out_card);
  if(out_dir != "" && out_dir != ".") gSystem->mkdir(out_dir, true);
  ws->writeToFile(out_ws, true);
  // the card refers to the workspace by its absolute path, so it works from any directory
  const TString out_ws_abs = (out_ws.BeginsWith("/")) ? out_ws : TString(gSystem->WorkingDirectory()) + "/" + out_ws;

  // Write the narrowed card: the same card with the new workspace, observation, and rates
  for(size_t i = 1; i < rates.size(); ++i) rates[i] = Form("%.6g", rates[i].Atof() * fractions[i-1]);

  // Combine drops processes without any yield in the range, and with them the nuisances that only
  // affect those processes, so such nuisances are removed from the nuisance groups as well
  std::set<TString> dropped;
  for(size_t i = 0; i < lines.size(); ++i) {
    const auto t = tokens(lines[i]);
    if(t.size() != processes.size() + 1 || (t[1] != "lnN" && t[1] != "shape")) continue;
    bool used = false;
    for(size_t j = 2; j < t.size(); ++j) used |= t[j] != "-" && rates[j-1].Atof() > 0.;
    if(!used) dropped.insert(t[0]);
  }
  std::ofstream outfile(out_card.Data());
  if(!outfile.is_open()) {
    cout << __func__ << ": Unable to write card " << out_card.Data() << endl;
    return 6;
  }
  for(size_t i = 0; i < lines.size(); ++i) {
    if((int) i == i_shapes) {
      auto t = shapes;
      t[3] = out_ws_abs;
      TString line; for(auto& s : t) line += s + " ";
      outfile << line.Strip(TString::kTrailing).Data() << endl;
    } else if((int) i == i_data) {
      outfile << Form("shapes data_obs %s %s workspace:data_obs_narrow", obs_name.Data(), out_ws_abs.Data()) << endl;
    } else if((int) i == i_obs) {
      outfile << Form("observation %.0f", data_narrow.sumEntries()) << endl;
    } else if((int) i == i_rate) {
      outfile << join(rates).Data() << endl;
    } else if(lines[i].BeginsWith("#") && lines[i].Contains("normalization floats freely")) {
      // "# <process> normalization floats freely (starting yield: <yield>)": the yield in the range
      const auto t = tokens(lines[i]);
      double yield = -1.;
      for(size_t j = 1; t.size() > 1 && j < processes.size(); ++j) {
        if(processes[j] == t[1]) yield = rates[j].Atof() * TString(lines[i](lines[i].Last(':')+1, lines[i].Length())).Atof();
      }
      if(yield >= 0.) outfile << Form("# %s normalization floats freely (starting yield in the range: %.4f)", t[1].Data(), yield) << endl;
      else            outfile << lines[i].Data() << endl;
    } else if(tokens(lines[i]).size() > 2 && tokens(lines[i])[1] == "group") {
      const auto t = tokens(lines[i]); // <name> group = <nuisances>
      TString line = t[0] + " group =";
      int nkept = 0;
      for(size_t j = 3; j < t.size(); ++j) {
        if(dropped.count(t[j])) continue;
        line += " " + t[j];
        ++nkept;
      }
      if(nkept > 0) outfile << line.Data() << endl;
    } else {
      outfile << lines[i].Data() << endl;
    }
    if(i == 1) outfile << Form("# Fit range narrowed to %.2f < p < %.2f MeV/c from %s\n", lo, hi, card.Data());
  }
  outfile.close();

  printf("narrow_model: [%.2f, %.2f] MeV/c, %i bins, observation %.0f\n", lo, hi, nbins, data_narrow.sumEntries());
  for(size_t i = 1; i < processes.size(); ++i)
    printf("narrow_model: %-12s fraction = %.4f rate = %s\n", processes[i].Data(), fractions[i-1], rates[i].Data());
  // the input file stays open: closing it deletes the observable data_narrow still refers to
  return 0;
}

#endif
