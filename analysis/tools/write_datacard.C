#ifndef __CONVANA_TOOLS_WRITEDATACARD__
#define __CONVANA_TOOLS_WRITEDATACARD__

#include <functional>

struct card_info_t {
  TString name_ = "";
  double  rate_ = 0.;
  int     selection_ = 0;
  bool    floating_ = false; // normalization floats freely via a <pdf>_norm variable (e.g. the envelope)
  double  mc_rel_unc_ = 0.;  // relative MC statistical uncertainty of the rate, 0 if not evaluated
  card_info_t(TString name, double rate, int selection, bool floating = false, double mc_rel_unc = 0.) :
    name_(name), rate_(rate), selection_(selection), floating_(floating), mc_rel_unc_(mc_rel_unc) {}
};

//---------------------------------------------------------------------------------------------------------------------------
// Rate (lnN) uncertainties: the one definition used by the shape and counting cards. These are
// named without the selection, so they are correlated across all channels of a merged card.
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
  // internal conversion fraction: RPC and RMC both use BR(internal RPC) / BR(RPC) = 0.0069 +- 0.00031
  // (physics.C), so one nuisance covers both
  sys.push_back(rate_sys_t("intN", 1.045, [](const TString& p) { return  p.Contains("_int"     ); }));
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
  // Smallest relative MC statistical uncertainty written to a card
  const double min_mc_rel_unc = 1.e-3;

  // One systematic line: its name, type (lnN, shape), and the entry for each process ("-" if none)
  struct sys_line_t {
    TString name_;
    TString type_;
    std::vector<TString> entries_;
  };

  // Column layout: the label column holds "<name> <type>", the process columns follow
  struct layout_t {
    int label_  = 15;
    int column_ = 10;
  };

  // Rate as written to the card: Combine multiplies a parametric shape by its <pdf>_norm variable,
  // so a freely floating process takes a unit rate and carries its yield in the workspace
  TString card_rate(const card_info_t& info) {
    return Form("%.4f", (info.floating_) ? 1. : info.rate_);
  }

  // Whether a process enters the model: Combine and combineCards.py drop zero-rate processes, and
  // with them every nuisance that only affects them, so such a process takes no uncertainties here
  // (a nuisance group naming a dropped nuisance fails the merged card). The rate is checked as
  // written, since a small positive rate is written, and read by Combine, as 0.
  bool contributes(const card_info_t& info) {
    return card_rate(info).Atof() > 0.;
  }

  bool is_signal(const card_info_t& info, const TString& signal_name) {
    return info.name_ == signal_name || info.name_ == "signal";
  }

  layout_t make_layout(const std::vector<card_info_t>& infos, const TString& bin_name,
                       const std::vector<sys_line_t>& lines) {
    layout_t layout;
    layout.column_ = std::max(10, bin_name.Length());
    for(auto& info : infos) layout.column_ = std::max(layout.column_, info.name_.Length());
    for(auto& line : lines) layout.label_ = std::max(layout.label_, line.name_.Length() + 1 + 5);
    return layout;
  }

  TString separator(const std::vector<card_info_t>& infos, const layout_t& layout) {
    return TString(std::string(layout.label_ + infos.size()*(layout.column_ + 1), '-'));
  }

  TString row(const TString& label, const std::vector<TString>& entries, const layout_t& layout) {
    TString line = Form("%-*s", layout.label_, label.Data());
    for(auto& entry : entries) line += Form(" %-*s", layout.column_, entry.Data());
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
                       const TString& bin_name, const layout_t& layout, const TString& sep) {
    std::vector<TString> bins, names, ids, rates;
    int nbkg = 0;
    for(auto& info : infos) {
      bins .push_back(bin_name);
      names.push_back(info.name_);
      ids  .push_back(Form("%i", (is_signal(info, signal_name)) ? 0 : ++nbkg));
      rates.push_back(card_rate(info));
    }
    outfile << sep.Data() << std::endl;
    outfile << row("bin"    , bins , layout).Data() << std::endl;
    outfile << row("process", names, layout).Data() << std::endl;
    outfile << row("process", ids  , layout).Data() << std::endl;
    outfile << row("rate"   , rates, layout).Data() << std::endl << std::endl;
  }

  // lnN rate uncertainties. A data-driven, freely floating process takes none, nor does a zero-rate
  // process, and an uncertainty that affects no process in the card is left out.
  std::vector<sys_line_t> rate_lines(const std::vector<card_info_t>& infos) {
    std::vector<sys_line_t> lines;
    for(auto& sys : rate_systematics()) {
      sys_line_t line{sys.name_, "lnN", {}};
      bool used = false;
      for(auto& info : infos) {
        const bool applies = !info.floating_ && contributes(info) && sys.applies_(info.name_);
        used |= applies;
        line.entries_.push_back((applies) ? TString(Form("%.3f", sys.kappa_)) : TString("-"));
      }
      if(used) lines.push_back(line);
    }
    return lines;
  }

  // MC statistical uncertainty of each process's rate, from the finite size of its simulated
  // sample. Each process in each channel has its own nuisance, <signal>_<selection>_<process>_mcstat,
  // so they are uncorrelated across processes and channels. Combine's autoMCStats does not apply
  // to the parametric (RooHistPdf) shapes, so this is the only MC statistics term in the model.
  std::vector<sys_line_t> mc_stat_lines(const std::vector<card_info_t>& infos, const TString& signal_name) {
    std::vector<sys_line_t> lines;
    for(size_t index = 0; index < infos.size(); ++index) {
      const auto& info = infos[index];
      if(info.floating_ || info.rate_ <= 0. || info.mc_rel_unc_ < min_mc_rel_unc) continue;
      sys_line_t line{Form("%s_%i_%s_mcstat", signal_name.Data(), info.selection_, info.name_.Data()), "lnN",
                      std::vector<TString>(infos.size(), "-")};
      line.entries_[index] = Form("%.3f", 1. + info.mc_rel_unc_);
      lines.push_back(line);
    }
    return lines;
  }

  // Normalization effect of the shape systematics as asymmetric lnN (kappa_down/kappa_up), for
  // the counting card, which has no PDFs to carry them
  std::vector<sys_line_t> shape_norm_lines(const std::vector<card_info_t>& infos, const shape_sys_map_t& shape_sys) {
    std::vector<sys_line_t> lines;
    for(auto& sys : shape_sys) {
      sys_line_t line{sys.first, "lnN", {}};
      bool used = false;
      for(auto& info : infos) {
        auto effect = sys.second.find(info.name_);
        const bool applies = effect != sys.second.end() && !info.floating_ && contributes(info)
                          && effect->second.complete() && effect->second.changes_norm();
        used |= applies;
        line.entries_.push_back((applies) ? TString(Form("%.3f/%.3f", effect->second.kappa_down(), effect->second.kappa_up()))
                                          : TString("-"));
      }
      if(used) lines.push_back(line);
    }
    return lines;
  }

  void write_lines(std::ofstream& outfile, const std::vector<sys_line_t>& lines, const layout_t& layout) {
    for(auto& line : lines) {
      outfile << row(Form("%s %s", line.name_.Data(), line.type_.Data()), line.entries_, layout).Data() << std::endl;
    }
  }

  // Nuisance groups, e.g. for impact breakdowns or freezing a class of uncertainties
  // (--freezeNuisanceGroups). combineCards.py merges groups of the same name across channels.
  void write_group(std::ofstream& outfile, const TString& group, const std::vector<TString>& names) {
    if(names.empty()) return;
    outfile << group.Data() << " group =";
    for(auto& name : names) outfile << " " << name.Data();
    outfile << std::endl;
  }

  std::vector<TString> line_names(const std::vector<sys_line_t>& lines) {
    std::vector<TString> names;
    for(auto& line : lines) names.push_back(line.name_);
    return names;
  }

  // Name of the scale factor of one process's yield
  TString scale_name(const TString& process) {
    return "scale_" + process;
  }

  void write_footer(std::ofstream& outfile, const std::vector<card_info_t>& infos) {
    // yield scale factor, useful for scanning livetimes
    outfile << "yieldScale rateParam * * 1." << std::endl;
    outfile << "nuisance edit freeze yieldScale" << std::endl;
    // per-process yield scale factors, e.g. combine --setParameters scale_dio=2. They are frozen, so they only rescale
    // the nominal yields. A zero-rate process is dropped by Combine, so it takes none.
    for(auto& info : infos) {
      if(!contributes(info)) continue;
      outfile << Form("%s rateParam * %s 1.", scale_name(info.name_).Data(), info.name_.Data()) << std::endl;
      outfile << Form("nuisance edit freeze %s", scale_name(info.name_).Data()) << std::endl;
    }
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
  // a control region without signal has no signal efficiency
  bool has_signal = false;
  for(auto& info : infos) has_signal |= datacard::is_signal(info, signal_name);
  RooRealVar* sig_eff = (RooRealVar*) ws->var("signal_eff");
  if(!sig_eff && has_signal) {
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

  // Systematic lines
  const auto rate_lines    = datacard::rate_lines(infos);
  const auto mc_stat_lines = datacard::mc_stat_lines(infos, signal_name);

  // shape uncertainties: Combine needs both the up and the down PDF of every process a line names
  std::vector<datacard::sys_line_t> shape_lines;
  std::vector<TString> shape_comments;
  for(auto& sys : shape_sys) {
    datacard::sys_line_t line{sys.first, "shape", {}};
    bool used = false;
    for(auto& info : infos) {
      bool applies = false;
      auto effect = sys.second.find(info.name_);
      if(effect != sys.second.end() && !info.floating_ && datacard::contributes(info)) {
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
      line.entries_.push_back((applies) ? "1" : "-");
    }
    if(used) shape_lines.push_back(line);
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

  //Make the combine card
  std::ofstream outfile;
  if(!datacard::open(outfile, outname)) return 7;
  std::vector<datacard::sys_line_t> all_lines = rate_lines;
  all_lines.insert(all_lines.end(), mc_stat_lines.begin(), mc_stat_lines.end());
  all_lines.insert(all_lines.end(), shape_lines.begin(), shape_lines.end());
  const auto layout = datacard::make_layout(infos, obs_name, all_lines);
  const TString filler = datacard::separator(infos, layout);
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

  datacard::write_processes(outfile, infos, signal_name, obs_name, layout, filler);

  // rate and MC statistical uncertainties
  outfile << filler.Data() << std::endl;
  datacard::write_lines(outfile, rate_lines, layout);
  datacard::write_lines(outfile, mc_stat_lines, layout);
  outfile << filler.Data() << std::endl;

  if(!shape_lines.empty()) {
    outfile << std::endl << filler.Data() << std::endl;
    for(auto& comment : shape_comments) outfile << comment.Data();
    datacard::write_lines(outfile, shape_lines, layout);
    outfile << filler.Data() << std::endl;
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

  // Nuisance groups
  std::vector<TString> shape_names = datacard::line_names(shape_lines);
  if(es_used) shape_names.push_back(es_name);
  datacard::write_group(outfile, "rate"  , datacard::line_names(rate_lines));
  datacard::write_group(outfile, "mcstat", datacard::line_names(mc_stat_lines));
  datacard::write_group(outfile, "shape" , shape_names);

  datacard::write_footer(outfile, infos);
  outfile.close();
  f->Close();

  return 0;
}

//---------------------------------------------------------------------------------------------------------------------------
// Counting card, written to outname as given. The shape systematics enter through their
// normalization effect on each process, evaluated over the full fit range.
int write_counting_datacard(TString signal_name, std::vector<card_info_t> infos,
                            TString outname, int nobs, double npot = -1., double livetime = -1., double nmuons = -1.,
                            double ref_br = -1., double signal_eff = -1.,
                            double xmin = 1., double xmax = -1.,
                            const shape_sys_map_t& shape_sys = {}) {

  if(infos.empty()) {
    cout << __func__ << ": No process information was given\n";
    return -1;
  }

  const int selection = infos[0].selection_; //assume fixed for all categories
  const TString obs_name = Form("obs_%i", selection);

  // A counting card has no workspace to carry a floating yield, so every rate is used as given
  for(auto& info : infos) info.floating_ = false;

  // Systematic lines
  const auto rate_lines       = datacard::rate_lines(infos);
  const auto mc_stat_lines    = datacard::mc_stat_lines(infos, signal_name);
  const auto shape_norm_lines = datacard::shape_norm_lines(infos, shape_sys);

  //Make the combine card
  std::ofstream outfile;
  if(!datacard::open(outfile, outname)) return 7;
  std::vector<datacard::sys_line_t> all_lines = rate_lines;
  all_lines.insert(all_lines.end(), mc_stat_lines.begin(), mc_stat_lines.end());
  all_lines.insert(all_lines.end(), shape_norm_lines.begin(), shape_norm_lines.end());
  const auto layout = datacard::make_layout(infos, obs_name, all_lines);
  const TString filler = datacard::separator(infos, layout);
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

  datacard::write_processes(outfile, infos, signal_name, obs_name, layout, filler);

  // rate, MC statistical, and shape normalization uncertainties
  outfile << filler.Data() << std::endl;
  datacard::write_lines(outfile, rate_lines, layout);
  datacard::write_lines(outfile, mc_stat_lines, layout);
  datacard::write_lines(outfile, shape_norm_lines, layout);
  outfile << filler.Data() << std::endl;

  // Nuisance groups
  datacard::write_group(outfile, "rate"  , datacard::line_names(rate_lines));
  datacard::write_group(outfile, "mcstat", datacard::line_names(mc_stat_lines));
  datacard::write_group(outfile, "shape" , datacard::line_names(shape_norm_lines));

  datacard::write_footer(outfile, infos);
  outfile.close();

  return 0;
}

#endif
