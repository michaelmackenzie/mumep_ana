// Create dataset histogram files with a modified selection, for fast selection sensitivity testing.
//
// For each dataset, the selection tree Ana/Hist/trs_<set>/tree is read and the model inputs (obs, t0, and
// a few track variables) are re-filled for tracks passing an additional cut, using the binning of the
// original Ana/Hist/trk_<set> histograms. The output files have the same names and layout as the inputs,
// so the model building can read them directly. Each tree is read once for all cut variations.
//
// sets     : selection sets to re-fill, e.g. "40" or "40 1040 2040" (sets must have a trs_<set> tree).
//            A token "IN>OUT" (e.g. "81>1081") reads the trk_IN histograms (binning) and the trs_IN tree and
//            writes the histograms to trk_OUT; a plain "81" means "81>81". Output sets are not copied by copy_all.
// cuts     : TTree cut string (e.g. "trk_qual > 0.9 && trk_pid > 0.95"), or a file with one "<tag> : <cut>"
//            variation per line ('#' comments allowed); each variation is written to <out_dir>/<tag>/
// out_dir  : output directory (default: selection_hists/ in the working directory)
// datasets : dataset keys to process, e.g. "mumep cosmic" (default: all non-data EventNtuple datasets).
//            A token that is not a dataset key but has a histogram file ConvAna.<func>.<token>.m<mode>.<type> in
//            the histogram path (CONVANA_HIST_PATH) is processed with the token as the dataset name, e.g.
//            "cele0b0s5r0101"; other unknown tokens are reported as issues.
// copy_all : copy the model inputs (obs, t0, and the DIO tail convolution inputs MC_GenE and dP) of all other
//            trk_<set> directories, so other sets (e.g. shape or control region sets) are still available
//            to the model building
// set_cuts : extra cuts per output set, ';'-separated "<out_set>:<cut>" entries (split on the first ':').
//            For an output set with an entry the cut applied is "(<variation cut>) && (<set cut>)".
//
// The cut is applied on top of the set's selection, so it can only tighten it. Use the outputs with:
//   ./build_model.py --hist-path <out_dir>[/<tag>] ... (or ./full_loop.py)
//
// Examples:
//   root.exe -q -b 'make_selection_hists.C("40", "trk_qual > 0.9", "selection_hists/qual_0p9")'
//   Split set 81 (no CRV veto) into a CRV-vetoed set 81 and a CRV-tagged set 1081:
//   root.exe -q -b 'make_selection_hists.C("81 81>1081", "trk_cluster > 0.", "selection_hists/crv", "cosmic", true,
//     "81: !(trk_crv_raw_dt > 0. && trk_crv_raw_dt < 150.) ; 1081: (trk_crv_raw_dt > 0. && trk_crv_raw_dt < 150.)")'

#include "../defaults.C"
#include "../datasets.C"

// Histograms re-filled from the selection tree: histogram name --> tree expression.
// The binning is taken from the histogram of the same name in the input trk_<set> directory.
const std::vector<std::pair<TString, TString>> selection_hist_vars_ = {
  {"obs"     , "trk_p"      }, // momentum observable used by the model
  {"t0"      , "trk_t0"     }, // time observable used by the model
  {"d0"      , "trk_d0"     },
  {"tanDip"  , "trk_tandip" },
  {"rMax"    , "trk_rmax"   },
  {"trkQual" , "trk_qual"   },
  {"pid"     , "trk_pid"    },
  {"clusterE", "trk_cluster"},
};
// Histograms copied for sets that are not re-filled
// (obs, t0, and the DIO tail convolution inputs MC_GenE and dP read by build_model.py)
const std::vector<TString> copied_hist_vars_ = {"obs", "t0", "MC_GenE", "dP"};

struct SelectionSet_t {
  int in_;  // set read (histogram binning and selection tree)
  int out_; // set written
};

struct SelectionVariation_t {
  TString tag_;
  TString cut_;
  TString out_dir_;
};

// Parse the cut argument into a list of variations
std::vector<SelectionVariation_t> parse_selection_variations(TString cuts, TString out_dir) {
  std::vector<SelectionVariation_t> variations;
  if(cuts != "" && !gSystem->AccessPathName(cuts)) {
    std::ifstream in(cuts.Data());
    std::string line;
    while(std::getline(in, line)) {
      TString l(line.c_str());
      if(l.Contains("#")) l = l(0, l.Index("#"));
      l = l.Strip(TString::kBoth);
      if(l == "") continue;
      const int colon = l.Index(":");
      if(colon <= 0) {
        cout << __func__ << ": Skipping line without \"<tag> : <cut>\" format: " << line << endl;
        continue;
      }
      TString tag = TString(l(0, colon)).Strip(TString::kBoth);
      TString cut = TString(l(colon+1, l.Length())).Strip(TString::kBoth);
      variations.push_back({tag, cut, out_dir + "/" + tag});
    }
  } else {
    variations.push_back({"", cuts, out_dir});
  }
  return variations;
}

// Parse the extra cuts per output set: ';'-separated "<out_set>:<cut>" entries, split on the first ':' only
std::map<int, TString> parse_set_cuts(TString set_cuts) {
  std::map<int, TString> cuts;
  std::unique_ptr<TObjArray> entries(set_cuts.Tokenize(";"));
  for(auto e : *entries) {
    TString entry = ((TObjString*) e)->GetString().Strip(TString::kBoth);
    const int colon = entry.Index(":");
    if(colon <= 0) {
      cout << __func__ << ": Skipping entry without \"<out_set>:<cut>\" format: " << entry.Data() << endl;
      continue;
    }
    const int set = TString(entry(0, colon)).Atoi();
    cuts[set] = TString(entry(colon+1, entry.Length())).Strip(TString::kBoth);
  }
  return cuts;
}

// Create a (nested) directory, returning the deepest directory (TDirectory::mkdir returns the top-most one)
TDirectory* make_dir(TFile* f, const TString& path) {
  f->mkdir(path, "", true);
  return f->GetDirectory(path);
}

// Copy the requested histograms of a directory into the output directory
void copy_hists(TDirectory* in, TDirectory* out, const std::vector<TString>& names) {
  for(const TString& name : names) {
    TH1* h = (TH1*) in->Get(name);
    if(!h) continue;
    out->cd();
    h->Clone(name)->Write();
  }
}

int make_selection_hists(TString sets = "40", TString cuts = "", TString out_dir = "selection_hists",
                         TString datasets = "", const bool copy_all = true, TString set_cuts = "") {
  use_evtana_ = true;
  set_evtana_defaults();
  init_dataset_info();

  //Parse the inputs
  std::vector<SelectionSet_t> set_list;
  {
    TString s = sets; s.ReplaceAll(",", " ");
    std::unique_ptr<TObjArray> tokens(s.Tokenize(" "));
    for(auto t : *tokens) {
      const TString token = ((TObjString*) t)->GetString();
      const int arrow = token.Index(">");
      if(arrow < 0) set_list.push_back({token.Atoi(), token.Atoi()});
      else          set_list.push_back({TString(token(0, arrow)).Atoi(), TString(token(arrow+1, token.Length())).Atoi()});
    }
  }
  if(set_list.empty()) {cout << __func__ << ": No sets given\n"; return 1;}
  const std::map<int, TString> set_cut_map = parse_set_cuts(set_cuts);
  int status = 0;
  for(auto& sc : set_cut_map) {
    const bool is_out = std::any_of(set_list.begin(), set_list.end(), [&](const SelectionSet_t& s) {return s.out_ == sc.first;});
    if(!is_out) {cout << __func__ << ": Extra cut given for set " << sc.first << ", which is not an output set\n"; ++status;}
  }
  // Effective cut of a variation for an output set
  auto effective_cut = [&](const TString& cut, const int out_set) {
    const TString base = (cut == "") ? "1" : cut;
    auto it = set_cut_map.find(out_set);
    return (it == set_cut_map.end()) ? base : TString("(" + base + ") && (" + it->second + ")");
  };

  const std::vector<SelectionVariation_t> variations = parse_selection_variations(cuts, out_dir);
  if(variations.empty()) {cout << __func__ << ": No cut variations found\n"; return 1;}

  //Dataset files to process, skipping aliases that point to the same file
  std::vector<std::pair<TString, TString>> files; // (dataset key, dataset name)
  {
    std::vector<TString> keys;
    if(datasets != "") {
      std::unique_ptr<TObjArray> tokens(datasets.Tokenize(" ,"));
      for(auto t : *tokens) keys.push_back(((TObjString*) t)->GetString());
    } else {
      for(auto& d : datasets_) if(!d.first.BeginsWith("data")) keys.push_back(d.first);
    }
    std::set<TString> seen;
    for(const TString& key : keys) {
      TString name = key; // a raw dataset name if not a dataset key
      if(datasets_.count(key)) {
        name = datasets_[key].name_;
        if(name == "") continue;
      } else if(gSystem->AccessPathName(Form("%sConvAna.%s.%s.m%i.%s", hist_path_, hist_func_, key.Data(), hist_mode_, file_type_.Data()))) {
        cout << __func__ << ": Warning: " << key.Data() << " is not a dataset key and has no histogram file in " << hist_path_ << endl;
        ++status;
        continue;
      }
      if(seen.count(name)) continue;
      seen.insert(name);
      files.push_back({key, name});
    }
  }

  cout << __func__ << ": Re-filling sets " << sets.Data() << " for " << files.size() << " dataset files\n";
  for(auto& v : variations) {
    gSystem->mkdir(v.out_dir_, true);
    cout << "  " << ((v.tag_ == "") ? "(nominal)" : v.tag_.Data()) << ": cut = \"" << v.cut_.Data()
         << "\" --> " << v.out_dir_.Data() << endl;
    for(auto& s : set_list) {
      cout << "    set " << s.in_ << " --> " << s.out_ << ": effective cut = \"" << effective_cut(v.cut_, s.out_).Data() << "\"\n";
    }
  }

  for(auto& file_info : files) {
    const TString in_name = Form("%sConvAna.%s.%s.m%i.%s", hist_path_, hist_func_, file_info.second.Data(), hist_mode_, file_type_.Data());
    TFile* fin = TFile::Open(in_name, "READ");
    if(!fin || fin->IsZombie()) {
      cout << __func__ << ": Unable to open " << in_name.Data() << endl;
      ++status;
      continue;
    }
    TDirectory* hist_dir = (TDirectory*) fin->Get(Form("%sHist", dir_path_.Data()));
    TTree* t_norm = (TTree*) fin->Get(Form("%sdata/Norm", dir_path_.Data()));
    if(!hist_dir || !t_norm) {
      cout << __func__ << ": Histogram directory or Norm tree not found in " << in_name.Data() << endl;
      fin->Close();
      ++status;
      continue;
    }

    //Create the output files with the normalization tree
    const TString base_name = gSystem->BaseName(in_name);
    std::vector<TFile*> fouts;
    for(auto& v : variations) {
      TFile* fout = new TFile(Form("%s/%s", v.out_dir_.Data(), base_name.Data()), "RECREATE");
      TDirectory* data_dir = make_dir(fout, Form("%sdata", dir_path_.Data()));
      data_dir->cd();
      t_norm->CloneTree(-1, "fast")->Write();
      make_dir(fout, Form("%sHist", dir_path_.Data()));
      fouts.push_back(fout);
    }

    //Copy the model inputs of the sets that are not re-filled
    if(copy_all) {
      for(auto key : *hist_dir->GetListOfKeys()) {
        const TString name = key->GetName();
        if(!name.BeginsWith("trk_")) continue;
        const int set = TString(name(4, name.Length())).Atoi();
        if(std::any_of(set_list.begin(), set_list.end(), [&](const SelectionSet_t& s) {return s.out_ == set;})) continue;
        TDirectory* in_set = (TDirectory*) hist_dir->Get(name);
        for(TFile* fout : fouts) {
          TDirectory* out_set = make_dir(fout, Form("%sHist/%s", dir_path_.Data(), name.Data()));
          copy_hists(in_set, out_set, copied_hist_vars_);
        }
      }
    }

    //Re-fill the requested sets from the selection trees
    for(const SelectionSet_t& set_info : set_list) {
      const int set = set_info.in_;
      const int out_set_id = set_info.out_;
      TDirectory* in_set = (TDirectory*) hist_dir->Get(Form("trk_%i", set));
      TTree* tree = (TTree*) hist_dir->Get(Form("trs_%i/tree", set));
      if(!in_set || !tree) {
        cout << __func__ << ": Set " << set << " histograms or selection tree not found in " << in_name.Data() << endl;
        ++status;
        continue;
      }

      //Output histograms per variation, with the input binning
      std::vector<std::vector<TH1*>> hists(variations.size());
      std::vector<TTreeFormula*> var_forms, cut_forms;
      for(auto& var : selection_hist_vars_) {
        TH1* h_in = (TH1*) in_set->Get(var.first);
        if(!h_in) {
          cout << __func__ << ": Histogram " << var.first.Data() << " not found for set " << set << " in " << in_name.Data() << endl;
          var_forms.push_back(nullptr);
          for(auto& h : hists) h.push_back(nullptr);
          continue;
        }
        var_forms.push_back(new TTreeFormula(Form("var_%s", var.first.Data()), var.second, tree));
        for(unsigned iv = 0; iv < variations.size(); ++iv) {
          TH1* h = (TH1*) h_in->Clone(var.first);
          h->SetDirectory(nullptr);
          h->Reset();
          hists[iv].push_back(h);
        }
      }
      for(unsigned iv = 0; iv < variations.size(); ++iv) {
        const TString cut = effective_cut(variations[iv].cut_, out_set_id);
        TTreeFormula* form = new TTreeFormula(Form("cut_%u", iv), cut, tree);
        if(form->GetNdim() == 0) {
          cout << __func__ << ": Invalid cut \"" << cut.Data() << "\"\n";
          fin->Close();
          return 1;
        }
        cut_forms.push_back(form);
      }
      float weight = 1.f;
      tree->SetBranchAddress("weight", &weight);

      //Loop through the selected tracks
      const Long64_t nentries = tree->GetEntries();
      std::vector<double> npass(variations.size(), 0.);
      for(Long64_t entry = 0; entry < nentries; ++entry) {
        tree->GetEntry(entry);
        for(unsigned iv = 0; iv < variations.size(); ++iv) {
          cut_forms[iv]->GetNdata();
          if(cut_forms[iv]->EvalInstance() == 0.) continue;
          npass[iv] += weight;
          for(unsigned ih = 0; ih < var_forms.size(); ++ih) {
            if(!var_forms[ih]) continue;
            var_forms[ih]->GetNdata();
            hists[iv][ih]->Fill(var_forms[ih]->EvalInstance(), weight);
          }
        }
      }

      //Write the outputs
      TH1* h_ref = (TH1*) in_set->Get("obs");
      const double nref = (h_ref) ? h_ref->Integral(0, h_ref->GetNbinsX()+1) : 0.;
      for(unsigned iv = 0; iv < variations.size(); ++iv) {
        TDirectory* out_set = make_dir(fouts[iv], Form("%sHist/trk_%i", dir_path_.Data(), out_set_id));
        out_set->cd();
        for(TH1* h : hists[iv]) if(h) {h->Write(); delete h;}
        printf("  %-12s set %4i %-12s: %10.1f / %10.1f selected (%.4f)\n", file_info.first.Data(), out_set_id,
               (variations[iv].tag_ == "") ? "" : variations[iv].tag_.Data(), npass[iv], nref, (nref > 0.) ? npass[iv]/nref : 0.);
      }
      for(auto f : var_forms) delete f;
      for(auto f : cut_forms) delete f;
      tree->ResetBranchAddresses();
    }

    for(TFile* fout : fouts) fout->Close();
    fin->Close();
  }

  if(status) cout << __func__ << ": Finished with " << status << " issues\n";
  return status;
}
