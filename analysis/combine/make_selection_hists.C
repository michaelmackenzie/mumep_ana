// Create dataset histogram files with a modified selection, for fast selection sensitivity testing.
//
// For each dataset, the selection tree Ana/Hist/trs_<set>/tree is read and the model inputs (obs, t0, and
// a few track variables) are re-filled for tracks passing an additional cut, using the binning of the
// original Ana/Hist/trk_<set> histograms. The output files have the same names and layout as the inputs,
// so the model building can read them directly. Each tree is read once for all cut variations.
//
// sets     : selection sets to re-fill, e.g. "40" or "40 1040 2040" (sets must have a trs_<set> tree)
// cuts     : TTree cut string (e.g. "trk_qual > 0.9 && trk_pid > 0.95"), or a file with one "<tag> : <cut>"
//            variation per line ('#' comments allowed); each variation is written to <out_dir>/<tag>/
// out_dir  : output directory (default: selection_hists/ in the working directory)
// datasets : dataset keys to process, e.g. "mumep cosmic" (default: all non-data EventNtuple datasets)
// copy_all : copy obs and t0 of all other trk_<set> directories, so other sets (e.g. shape or control
//            region sets) are still available to the model building
//
// The cut is applied on top of the set's selection, so it can only tighten it. Use the outputs with:
//   ./build_model.py --hist-path <out_dir>[/<tag>] ... (or ./full_loop.py)
//
// Example:
//   root.exe -q -b 'make_selection_hists.C("40", "trk_qual > 0.9", "selection_hists/qual_0p9")'

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
const std::vector<TString> copied_hist_vars_ = {"obs", "t0"};

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
                         TString datasets = "", const bool copy_all = true) {
  use_evtana_ = true;
  set_evtana_defaults();
  init_dataset_info();

  //Parse the inputs
  std::vector<int> set_list;
  {
    TString s = sets; s.ReplaceAll(",", " ");
    std::unique_ptr<TObjArray> tokens(s.Tokenize(" "));
    for(auto t : *tokens) set_list.push_back(((TObjString*) t)->GetString().Atoi());
  }
  if(set_list.empty()) {cout << __func__ << ": No sets given\n"; return 1;}

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
      auto info = get_dataset_info(key);
      if(info.name_ == "" || seen.count(info.name_)) continue;
      seen.insert(info.name_);
      files.push_back({key, info.name_});
    }
  }

  cout << __func__ << ": Re-filling sets " << sets.Data() << " for " << files.size() << " dataset files\n";
  for(auto& v : variations) {
    gSystem->mkdir(v.out_dir_, true);
    cout << "  " << ((v.tag_ == "") ? "(nominal)" : v.tag_.Data()) << ": cut = \"" << v.cut_.Data()
         << "\" --> " << v.out_dir_.Data() << endl;
  }

  int status = 0;
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
        if(std::find(set_list.begin(), set_list.end(), set) != set_list.end()) continue;
        TDirectory* in_set = (TDirectory*) hist_dir->Get(name);
        for(TFile* fout : fouts) {
          TDirectory* out_set = make_dir(fout, Form("%sHist/%s", dir_path_.Data(), name.Data()));
          copy_hists(in_set, out_set, copied_hist_vars_);
        }
      }
    }

    //Re-fill the requested sets from the selection trees
    for(const int set : set_list) {
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
        const TString cut = (variations[iv].cut_ == "") ? "1" : variations[iv].cut_;
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
        TDirectory* out_set = make_dir(fouts[iv], Form("%sHist/trk_%i", dir_path_.Data(), set));
        out_set->cd();
        for(TH1* h : hists[iv]) if(h) {h->Write(); delete h;}
        printf("  %-12s set %4i %-12s: %10.1f / %10.1f selected (%.4f)\n", file_info.first.Data(), set,
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
