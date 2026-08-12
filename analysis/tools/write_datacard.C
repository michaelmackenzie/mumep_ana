#ifndef __CONVANA_TOOLS_WRITEDATACARD__
#define __CONVANA_TOOLS_WRITEDATACARD__

struct card_info_t {
  TString name_ = "";
  double  rate_ = 0.;
  int     selection_ = 0;
  card_info_t(TString name, double rate, int selection) : name_(name), rate_(rate), selection_(selection) {}
};

int write_datacard(TString signal_name, std::vector<card_info_t> infos, TString file_in, map<TString, map<TString, bool>> sys_map,
                   TString outname = "") {

  if(infos.empty()) {
    cout << __func__ << ": Not process information was given\n";
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
  RooDataHist* data = (RooDataHist*) ws->data("data_obs");
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

  //Make the combine card
  gSystem->Exec("[ ! -d datacards ] && mkdir datacards");
  std::ofstream outfile;
  outfile.open(outname.Data());
  TString filler = std::string((infos.size()+1)*10 + 15, '-');
  outfile << "# -*- mode:tcl; eval: (whitespace-mode 0) -*-\n# Auto-generated Combine data card\n";
  outfile << Form("# R_mue used for signal: %.3e\n", ref_br->getVal());
  if(npot) outfile << Form("# N(POT): %.3e\n", npot->getVal());
  if(livetime) outfile << Form("# Livetime: %.3e\n", livetime->getVal());
  if(nmuons) outfile << Form("# N(muons): %.3e\n", nmuons->getVal());
  if(sig_eff) outfile << Form("# Signal efficiency: %.3e\n", sig_eff->getVal());
  outfile << filler.Data() << std::endl;
  outfile << "\nimax 1 #number of bins\njmax * #number of processes\nkmax * #number of systematics\n\n";
  outfile << filler.Data() << std::endl;

  outfile << Form("shapes * %s %s workspace:%s_%i_$PROCESS_pdf workspace:%s_%i_$PROCESS_pdf_$SYSTEMATIC\n", obs_name, file_in.Data(),
                  signal_name.Data(), selection, signal_name.Data(), selection);
  outfile << Form("shapes data_obs %s %s workspace:data_obs\n\n", obs_name, file_in.Data());
  outfile << filler.Data() << std::endl;
  // outfile << "bin " << obs_name << endl;
  outfile << "observation " << Form("%.0f", data->sumEntries()) << std::endl << std::endl;
  // outfile << "observation 0\n\n";

  TString bins   = Form("%-15s", "bin"    );
  TString proc_n = Form("%-15s", "process");
  TString proc_i = Form("%-15s", "process");
  TString rates  = Form("%-15s", "rate"   );
  std::vector<TString> systematics;
  int ncats = 0;
  for(size_t index = 0; index < infos.size(); ++index) {
    auto& info = infos[index];
    // Check for the input information in the workspace
    const char* pdf_name = Form("%s_%i_%s_pdf", signal_name.Data(), selection, info.name_.Data());
    auto pdf = ws->pdf(pdf_name);
    if(!pdf) {
      cout << "PDF " << pdf_name << " not found in file " << file_in.Data() << endl;
      continue;
    }

    const bool is_signal = info.name_ == signal_name || info.name_ == "signal";
    if(is_signal) {
      ++ncats;
    }
    const bool is_cosmic = info.name_.Contains("cosmic");
    const bool is_dio    = info.name_.Contains("dio");
    const bool is_rpc    = info.name_.Contains("rpc");
    const bool is_pbar   = info.name_.Contains("pbar");
    const bool is_rmc    = info.name_.Contains("rmc");
    const int category = (is_signal) ? 0 : ncats;
    bins += Form(" %-10s", obs_name);
    proc_n += Form(" %-10s", info.name_.Data());
    proc_i += Form(" %-10i", category);
    rates  += Form(" %-10.4f", info.rate_);

    if(index == 0) {
      systematics.push_back(Form("%-10s %-4s", "lumi", "lnN"));
      systematics.push_back(Form("%-10s %-4s", "csmN", "lnN"));
      systematics.push_back(Form("%-10s %-4s", "dioN", "lnN"));
      systematics.push_back(Form("%-10s %-4s", "rpcN", "lnN"));
      systematics.push_back(Form("%-10s %-4s", "pbrN", "lnN"));
      // systematics.push_back(Form("%-10s %-4s", "rmcN", "lnN"));
    }
    if(!is_cosmic) systematics[0] += Form(" %-10.3f", 1.1);
    else           systematics[0] += Form(" %-10s", "-");
    if(is_cosmic ) systematics[1] += Form(" %-10.3f", 1.2);
    else           systematics[1] += Form(" %-10s", "-");
    if(is_dio    ) systematics[2] += Form(" %-10.3f", 1.025);
    else           systematics[2] += Form(" %-10s", "-");
    if(is_rpc    ) systematics[3] += Form(" %-10.3f", 1.27);
    else           systematics[3] += Form(" %-10s", "-");
    if(is_pbar   ) systematics[4] += Form(" %-10.3f", 2.0);
    else           systematics[4] += Form(" %-10s", "-");
    // // Systematics
    // outfile << "lumi   lnN     1.1        1.1         -         1.1        1.1        1.1\n";
    // outfile << "sigN   lnN     1.04        -          -          -          -          -\n";
    // outfile << "dioN   lnN      -        1.025        -          -          -          -\n";
    // outfile << "csmN   lnN      -          -         1.2         -          -          -\n";
    // outfile << "pbrN   lnN      -          -          -         2.0         -          -\n";
    // outfile << "rpcN   lnN      -          -          -          -        1.093      1.093\n";
    // outfile << "rpiN   lnN      -          -          -          -          -        1.045\n";
    // outfile << "pion   lnN      -          -          -          -         1.27       1.27\n";
    // outfile << "MomScale shape  1          1          -          -          -          -\n";
  }

  outfile << filler.Data() << std::endl;
  outfile << bins.Data() << std::endl;
  outfile << proc_n.Data() << std::endl;
  outfile << proc_i.Data() << std::endl;
  outfile << rates.Data() << std::endl << std::endl;

  // rate uncertainties
  outfile << filler.Data() << std::endl;
  for(auto sys : systematics) outfile << sys.Data() << std::endl;
  outfile << filler.Data() << std::endl;

  // shape uncertainties
  if(sys_map.size() > 0) outfile << std::endl << filler.Data() << std::endl;
  for(auto sys : sys_map) {
    TString line = Form("%-10s shape", sys.first.Data());
    for(size_t index = 0; index < infos.size(); ++index) {
      auto& info = infos[index];
      if(sys.second[info.name_]) line += Form(" %-10i", 1);
      else                       line += Form(" %-10s", "-");
    }
    outfile << line.Data() << std::endl;
  }
  if(sys_map.size() > 0) outfile << filler.Data() << std::endl;

  // outfile << "\n* autoMCStats 0\n"; //MC statistics uncertainty

  // yield scale factor, useful for scanning livetimes
  outfile << "yieldScale rateParam * * 1." << std::endl;
  outfile << "nuisance edit freeze yieldScale" << std::endl;

  outfile.close();

  return 0;
}

int write_counting_datacard(TString signal_name, std::vector<card_info_t> infos,
                            TString outname, int nobs, double npot = -1., double livetime = -1., double nmuons = -1.,
                            double ref_br = -1., double signal_eff = -1.,
                            double xmin = 1., double xmax = -1.) {

  if(infos.empty()) {
    cout << __func__ << ": Not process information was given\n";
    return -1;
  }

  const int selection = infos[0].selection_; //assume fixed for all categories
  const char* obs_name = Form("obs_%i", selection);

  //Make the combine card
  gSystem->Exec("[ ! -d datacards ] && mkdir datacards");
  std::ofstream outfile;
  outfile.open(outname.Data());
  TString filler = std::string((infos.size()+1)*10 + 15, '-');
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

  TString bins   = Form("%-15s", "bin"    );
  TString proc_n = Form("%-15s", "process");
  TString proc_i = Form("%-15s", "process");
  TString rates  = Form("%-15s", "rate"   );
  std::vector<TString> systematics;
  int ncats = 0;
  for(size_t index = 0; index < infos.size(); ++index) {
    auto& info = infos[index];

    const bool is_signal = info.name_ == signal_name || info.name_ == "signal";
    if(!is_signal) {
      ++ncats;
    }
    const bool is_cosmic = info.name_.Contains("cosmic");
    const bool is_dio    = info.name_.Contains("dio");
    const bool is_rpc    = info.name_.Contains("rpc");
    const bool is_pbar   = info.name_.Contains("pbar");
    const bool is_rmc    = info.name_.Contains("rmc");
    const int category = (is_signal) ? 0 : ncats;
    bins += Form(" %-10s", obs_name);
    proc_n += Form(" %-10s", info.name_.Data());
    proc_i += Form(" %-10i", category);
    rates  += Form(" %-10.4f", info.rate_);

    if(index == 0) {
      systematics.push_back(Form("%-10s %-4s", "lumi", "lnN"));
      systematics.push_back(Form("%-10s %-4s", "csmN", "lnN"));
      systematics.push_back(Form("%-10s %-4s", "dioN", "lnN"));
      systematics.push_back(Form("%-10s %-4s", "rpcN", "lnN"));
      systematics.push_back(Form("%-10s %-4s", "pbrN", "lnN"));
    }
    if(!is_cosmic) systematics[0] += Form(" %-10.3f", 1.1);
    else           systematics[0] += Form(" %-10s", "-");
    if(is_cosmic ) systematics[1] += Form(" %-10.3f", 1.2);
    else           systematics[1] += Form(" %-10s", "-");
    if(is_dio    ) systematics[2] += Form(" %-10.3f", 1.025);
    else           systematics[2] += Form(" %-10s", "-");
    if(is_rpc    ) systematics[3] += Form(" %-10.3f", 1.27);
    else           systematics[3] += Form(" %-10s", "-");
    if(is_pbar   ) systematics[4] += Form(" %-10.3f", 2.0);
    else           systematics[4] += Form(" %-10s", "-");
  }

  outfile << filler.Data() << std::endl;
  outfile << bins.Data() << std::endl;
  outfile << proc_n.Data() << std::endl;
  outfile << proc_i.Data() << std::endl;
  outfile << rates.Data() << std::endl << std::endl;

  // rate uncertainties
  outfile << filler.Data() << std::endl;
  for(auto sys : systematics) outfile << sys.Data() << std::endl;
  outfile << filler.Data() << std::endl;

  // yield scale factor, useful for scanning livetimes
  outfile << "yieldScale rateParam * * 1." << std::endl;
  outfile << "nuisance edit freeze yieldScale" << std::endl;

  outfile.close();

  return 0;
}

#endif
