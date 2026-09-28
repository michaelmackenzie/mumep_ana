#ifndef __CONVANA_TOOLS_WRITEDATACARD__
#define __CONVANA_TOOLS_WRITEDATACARD__

#include <functional>

struct card_info_t {
  TString name_ = "";
  double  rate_ = 0.;
  int     selection_ = 0;
  bool    floating_ = false; // normalization floats freely via a <pdf>_norm variable (e.g. the envelope)
  card_info_t(TString name, double rate, int selection, bool floating = false) :
    name_(name), rate_(rate), selection_(selection), floating_(floating) {}
};

//---------------------------------------------------------------------------------------------------------------------------
// Rate (lnN) uncertainties: the one definition used by the shape and counting cards
struct rate_sys_t {
  TString name_;
  double  kappa_; // lnN kappa, 1 + relative uncertainty
  std::function<bool(const TString&)> applies_; // process name --> whether the uncertainty applies
  rate_sys_t(TString name, double kappa, std::function<bool(const TString&)> applies) :
    name_(name), kappa_(kappa), applies_(applies) {}
};

// Process families are named <base> or <base>_<variant>, e.g. rmc_ext_0n
bool is_process(const TString& name, const TString& base) {
  return name == base || name.BeginsWith(base + "_");
}

std::vector<rate_sys_t> rate_systematics() {
  std::vector<rate_sys_t> sys;
  sys.push_back(rate_sys_t("lumi", 1.1  , [](const TString& p) { return !is_process(p, "cosmic"); })); // beam intensity, incl. the signal
  sys.push_back(rate_sys_t("csmN", 1.2  , [](const TString& p) { return  is_process(p, "cosmic"); }));
  sys.push_back(rate_sys_t("dioN", 1.025, [](const TString& p) { return  is_process(p, "dio"   ); }));
  sys.push_back(rate_sys_t("rpcN", 1.27 , [](const TString& p) { return  is_process(p, "rpc"   ); }));
  sys.push_back(rate_sys_t("pbrN", 2.0  , [](const TString& p) { return  is_process(p, "pbar"  ); }));
  sys.push_back(rate_sys_t("rmcN", 1.079, [](const TString& p) { return  is_process(p, "rmc"   ); })); // 1.40 +- 0.11 (TRIUMF, 1999)
  sys.push_back(rate_sys_t("intN", 1.045, [](const TString& p) { return  p.Contains("_int"     ); })); // internal conversion: 0.0069 +- 0.00031
  return sys;
}

//---------------------------------------------------------------------------------------------------------------------------
// Effect of one shape systematic on one process
struct shape_sys_t {
  bool   has_up_    = false;
  bool   has_down_  = false;
  double norm_nom_  = 0.;
  double norm_up_   = 0.;
  double norm_down_ = 0.;
  bool   complete    () const { return has_up_ && has_down_; }
  double kappa_up    () const { return (norm_nom_ > 0.) ? norm_up_   / norm_nom_ : 1.; }
  double kappa_down  () const { return (norm_nom_ > 0.) ? norm_down_ / norm_nom_ : 1.; }
  bool   changes_norm() const { return std::fabs(kappa_up() - 1.) > 1.e-3 || std::fabs(kappa_down() - 1.) > 1.e-3; }
};
typedef std::map<TString, std::map<TString, shape_sys_t>> shape_sys_map_t; // systematic --> process --> effect

//---------------------------------------------------------------------------------------------------------------------------
// Card writing helpers shared by the shape and counting cards
namespace datacard {
  bool is_signal(const card_info_t& info, const TString& signal_name) {
    return info.name_ == signal_name || info.name_ == "signal";
  }

  int column_width(const std::vector<card_info_t>& infos, const TString& bin_name) {
    int width = std::max(10, bin_name.Length());
    for(auto& info : infos) width = std::max(width, info.name_.Length());
    return width;
  }

  TString separator(const std::vector<card_info_t>& infos, const int width) {
    return TString(std::string(15 + infos.size()*(width + 1), '-'));
  }

  TString row(const TString& label, const std::vector<TString>& entries, const int width) {
    TString line = Form("%-15s", label.Data());
    for(auto& entry : entries) line += Form(" %-*s", width, entry.Data());
    return line;
  }

  bool open(std::ofstream& outfile, const TString& outname) {
    const TString dir = gSystem->GetDirName(outname);
    if(dir != "" && dir != ".") gSystem->mkdir(dir, true);
    outfile.open(outname.Data());
    if(!outfile.is_open()) {
      cout << "datacard::" << __func__ << ": Unable to open " << outname.Data() << " for writing\n";
      return false;
    }
    return true;
  }

  // bin/process/rate block: the signal is process 0, backgrounds are numbered 1..N
  void write_processes(std::ofstream& outfile, const std::vector<card_info_t>& infos, const TString& signal_name,
                       const TString& bin_name, const int width, const TString& sep) {
    std::vector<TString> bins, names, ids, rates;
    int nbkg = 0;
    for(auto& info : infos) {
      bins .push_back(bin_name);
      names.push_back(info.name_);
      ids  .push_back(Form("%i", (is_signal(info, signal_name)) ? 0 : ++nbkg));
      // Combine multiplies a parametric shape by its <pdf>_norm variable, so a freely floating
      // process takes a unit rate here and carries its yield in the workspace
      rates.push_back(Form("%.4f", (info.floating_) ? 1. : info.rate_));
    }
    outfile << sep.Data() << std::endl;
    outfile << row("bin"    , bins , width).Data() << std::endl;
    outfile << row("process", names, width).Data() << std::endl;
    outfile << row("process", ids  , width).Data() << std::endl;
    outfile << row("rate"   , rates, width).Data() << std::endl << std::endl;
  }

  // lnN rate uncertainties. A data-driven, freely floating process takes none, and an
  // uncertainty that affects no process in the card is left out.
  void write_rate_systematics(std::ofstream& outfile, const std::vector<card_info_t>& infos, const int width, const TString& sep) {
    outfile << sep.Data() << std::endl;
    for(auto& sys : rate_systematics()) {
      std::vector<TString> entries;
      bool used = false;
      for(auto& info : infos) {
        const bool applies = !info.floating_ && sys.applies_(info.name_);
        used |= applies;
        entries.push_back((applies) ? TString(Form("%.3f", sys.kappa_)) : TString("-"));
      }
      if(used) outfile << row(Form("%-10s %-4s", sys.name_.Data(), "lnN"), entries, width).Data() << std::endl;
    }
    outfile << sep.Data() << std::endl;
  }

  void write_footer(std::ofstream& outfile) {
    // yield scale factor, useful for scanning livetimes
    outfile << "yieldScale rateParam * * 1." << std::endl;
    outfile << "nuisance edit freeze yieldScale" << std::endl;
  }
}

//---------------------------------------------------------------------------------------------------------------------------
// Shape card for the workspace written by build_model. The up/down PDFs of a shape systematic are
// expected as <nominal pdf>_<sys>Up/Down. Combine ignores the yield of PDF shape variations and does
// not allow a shape and an lnN line of the same name, so a normalization effect has to be carried
// by a <nominal pdf>_norm term in the workspace (build_model adds it). use_es enables the
// energy-scale nuisance of the function models.
int write_datacard(TString signal_name, std::vector<card_info_t> infos, TString file_in, const shape_sys_map_t& shape_sys,
                   TString outname = "", std::vector<TString> extra_lines = {}, const bool use_es = true) {

  if(infos.empty()) {
    cout << __func__ << ": No process information was given\n";
    return -1;
  }

  // Determine the outfile name
  if(outname == "") {
    outname = file_in;
    outname.ReplaceAll(".root", ".txt");
    if(outname.Contains("/")) outname = outname(outname.Last('/')+1, outname.Sizeof());
    outname.ReplaceAll("workspace", "combine");
  }
  outname = "datacards/" + outname;

  // Open the input file
  TFile* f = TFile::Open(file_in.Data(), "READ");
  if(!f) return 1;
  RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
  if(!ws) {
    cout << "Workspace not found in file " << file_in.Data() << endl;
    return 2;
  }
  RooRealVar* ref_br = (RooRealVar*) ws->var("ref_signal_br");
  if(!ref_br) {
    cout << "Reference signal branching fraction not found in file " << file_in.Data() << endl;
    return 3;
  }
  RooRealVar* npot = (RooRealVar*) ws->var("npot");
  if(!npot) {
    cout << "Reference N(POT) not found in file " << file_in.Data() << endl;
  }
  RooRealVar* livetime = (RooRealVar*) ws->var("livetime");
  if(!livetime) {
    cout << "Reference livetime not found in file " << file_in.Data() << endl;
  }
  RooRealVar* nmuons = (RooRealVar*) ws->var("nmuons");
  if(!nmuons) {
    cout << "Reference N(muons) not found in file " << file_in.Data() << endl;
  }
  RooRealVar* sig_eff = (RooRealVar*) ws->var("signal_eff");
  if(!sig_eff) {
    cout << "Reference signal efficiency not found in file " << file_in.Data() << endl;
  }
  RooAbsData* data = ws->data("data_obs");
  if(!data) {
    cout << "No data found in file " << file_in.Data() << endl;
    return 4;
  }

  const int selection = infos[0].selection_; //assume fixed for all categories
  RooRealVar* obs = (RooRealVar*) ws->var(Form("obs_%i", selection));
  if(!obs) {
    cout << "Observable is not defined!\n";
    return 5;
  }
  const char* obs_name = obs->GetName();

  // Every process needs its PDF, or the card columns and the workspace disagree
  auto pdf_name = [&](const TString& proc) { return TString(Form("%s_%i_%s_pdf", signal_name.Data(), selection, proc.Data())); };
  for(auto& info : infos) {
    if(!ws->pdf(pdf_name(info.name_))) {
      cout << __func__ << ": PDF " << pdf_name(info.name_).Data() << " not found in file " << file_in.Data() << endl;
      return 6;
    }
  }

  //Make the combine card
  std::ofstream outfile;
  if(!datacard::open(outfile, outname)) return 7;
  const int width = datacard::column_width(infos, obs_name);
  const TString filler = datacard::separator(infos, width);
  outfile << "# -*- mode:tcl; eval: (whitespace-mode 0) -*-\n# Auto-generated Combine data card\n";
  outfile << Form("# R_mue used for signal: %.3e\n", ref_br->getVal());
  if(npot) outfile << Form("# N(POT): %.3e\n", npot->getVal());
  if(livetime) outfile << Form("# Livetime: %.3e\n", livetime->getVal());
  if(nmuons) outfile << Form("# N(muons): %.3e\n", nmuons->getVal());
  if(sig_eff) outfile << Form("# Signal efficiency: %.3e\n", sig_eff->getVal());
  for(auto& info : infos) {
    // A floating process carries its yield in the workspace, so the card rate below is 1
    if(info.floating_) outfile << Form("# %s normalization floats freely (starting yield: %.4f)\n",
                                       info.name_.Data(), info.rate_);
  }
  outfile << filler.Data() << std::endl;
  outfile << "\nimax 1 #number of bins\njmax * #number of processes\nkmax * #number of systematics\n\n";
  outfile << filler.Data() << std::endl;

  outfile << Form("shapes * %s %s workspace:%s_%i_$PROCESS_pdf workspace:%s_%i_$PROCESS_pdf_$SYSTEMATIC\n", obs_name, file_in.Data(),
                  signal_name.Data(), selection, signal_name.Data(), selection);
  outfile << Form("shapes data_obs %s %s workspace:data_obs\n\n", obs_name, file_in.Data());
  outfile << filler.Data() << std::endl;
  outfile << "observation " << Form("%.0f", data->sumEntries()) << std::endl << std::endl;

  datacard::write_processes(outfile, infos, signal_name, obs_name, width, filler);

  // rate uncertainties
  datacard::write_rate_systematics(outfile, infos, width, filler);

  // shape uncertainties: Combine needs both the up and the down PDF of every process a line names
  std::vector<TString> shape_lines;
  std::vector<TString> shape_comments;
  for(auto& sys : shape_sys) {
    std::vector<TString> entries;
    bool used = false;
    for(auto& info : infos) {
      bool applies = false;
      auto effect = sys.second.find(info.name_);
      if(effect != sys.second.end() && !info.floating_) {
        const TString base = pdf_name(info.name_);
        const bool found_up   = ws->pdf(base + "_" + sys.first + "Up"  );
        const bool found_down = ws->pdf(base + "_" + sys.first + "Down");
        applies = effect->second.complete() && found_up && found_down;
        if(!applies) {
          cout << __func__ << ": Dropping shape systematic " << sys.first.Data() << " for " << info.name_.Data()
               << ", its " << ((found_up) ? "down" : (found_down) ? "up" : "up and down") << " PDF is missing\n";
        } else if(effect->second.changes_norm()) {
          shape_comments.push_back(Form("# %s changes the %s yield by %+.2f%% / %+.2f%% (up / down) via %s_norm\n",
                                        sys.first.Data(), info.name_.Data(),
                                        100.*(effect->second.kappa_up() - 1.), 100.*(effect->second.kappa_down() - 1.),
                                        base.Data()));
        }
      }
      used |= applies;
      entries.push_back((applies) ? "1" : "-");
    }
    if(used) shape_lines.push_back(datacard::row(Form("%-10s %-4s", sys.first.Data(), "shape"), entries, width));
  }
  if(!shape_lines.empty()) {
    outfile << std::endl << filler.Data() << std::endl;
    for(auto& comment : shape_comments) outfile << comment.Data();
    for(auto& line : shape_lines) outfile << line.Data() << std::endl;
    outfile << filler.Data() << std::endl;
  }

  // constrained params: the energy-scale nuisance is the <process>_<selection>_es variable the
  // signal/DIO function models shift with (signal_model.C, background_model.C). Combine makes it
  // float when a param line names it; histogram models have no such variable, and a param line
  // would then only add a disconnected nuisance, so it is written only when the model uses it.
  // Without the param line the variable stays constant at its nominal value.
  const TString es_name = Form("%s_%i_es", signal_name.Data(), selection);
  RooRealVar* es_var = ws->var(es_name.Data());
  bool es_used = false;
  if(use_es && es_var) {
    for(auto& info : infos) {
      if(ws->pdf(pdf_name(info.name_))->dependsOn(*es_var)) { es_used = true; break; }
    }
  }
  if(es_used) {
    outfile << es_name.Data() << " param 0.0  1.0\n";
    outfile << filler.Data() << std::endl;
  }

  // Additional card lines, e.g. the discrete index Combine profiles over for an envelope
  if(!extra_lines.empty()) {
    outfile << std::endl;
    for(auto& line : extra_lines) outfile << line.Data() << std::endl;
    outfile << filler.Data() << std::endl;
  }

  datacard::write_footer(outfile);
  outfile.close();
  f->Close();

  return 0;
}

//---------------------------------------------------------------------------------------------------------------------------
// Counting card, written to outname as given
int write_counting_datacard(TString signal_name, std::vector<card_info_t> infos,
                            TString outname, int nobs, double npot = -1., double livetime = -1., double nmuons = -1.,
                            double ref_br = -1., double signal_eff = -1.,
                            double xmin = 1., double xmax = -1.) {

  if(infos.empty()) {
    cout << __func__ << ": No process information was given\n";
    return -1;
  }

  const int selection = infos[0].selection_; //assume fixed for all categories
  const TString obs_name = Form("obs_%i", selection);

  //Make the combine card
  std::ofstream outfile;
  if(!datacard::open(outfile, outname)) return 7;
  const int width = datacard::column_width(infos, obs_name);
  const TString filler = datacard::separator(infos, width);
  outfile << "# -*- mode:tcl; eval: (whitespace-mode 0) -*-\n# Auto-generated Combine data card\n";
  if(ref_br > 0.) outfile << Form("# R_mue used for signal: %.3e\n", ref_br);
  if(npot > 0.) outfile << Form("# N(POT): %.3e\n", npot);
  if(livetime > 0.) outfile << Form("# Livetime: %.3e\n", livetime);
  if(nmuons > 0.) outfile << Form("# N(muons): %.3e\n", nmuons);
  if(signal_eff > 0.) outfile << Form("# Signal efficiency: %.3e\n", signal_eff);
  if(xmin < xmax) outfile << Form("# Selection: %.2f < p < %.2f MeV/c\n", xmin, xmax);
  outfile << filler.Data() << std::endl;
  outfile << "\nimax 1 #number of bins\njmax * #number of processes\nkmax * #number of systematics\n\n";
  outfile << filler.Data() << std::endl;

  outfile << "observation " << Form("%i", nobs) << std::endl << std::endl;

  // A counting card has no workspace to carry a floating yield, so every rate is used as given
  for(auto& info : infos) info.floating_ = false;
  datacard::write_processes(outfile, infos, signal_name, obs_name, width, filler);

  // rate uncertainties
  datacard::write_rate_systematics(outfile, infos, width, filler);

  datacard::write_footer(outfile);
  outfile.close();

  return 0;
}

#endif
