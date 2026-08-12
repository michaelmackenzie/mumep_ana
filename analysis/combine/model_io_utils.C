#ifndef __CONVANA_ANALYSIS_MODEL_IO_UTILS__
#define __CONVANA_ANALYSIS_MODEL_IO_UTILS__

#include "../defaults.C"
#include "../datasets.C"

Long64_t read_norm_tree_entries(TFile* f, const TString& process) {
  if(!f) return 0;
  TTree* t_norm = (TTree*) f->Get(Form("%sdata/Norm", dir_path_.Data()));
  if(!t_norm) {
    cout << __func__ << ": Normalization tree for process " << process.Data() << " not found\n";
    return 0;
  }

  Long64_t nseen(0), ntotal(0);
  t_norm->SetBranchAddress("nseen", &nseen);
  for(Long64_t entry = 0; entry < t_norm->GetEntries(); ++entry) {
    t_norm->GetEntry(entry);
    ntotal += nseen;
  }
  return ntotal;
}

void apply_expected_event_scaling(TH1* h,
                                  TFile* f,
                                  const TString& process,
                                  const Long64_t nexpect) {
  if(!h || !f || nexpect <= 0) return;

  const Long64_t ntotal = read_norm_tree_entries(f, process);
  if(ntotal <= 0) {
    cout << __func__ << ": No normalization entries found for process " << process.Data() << endl;
    return;
  }

  if(nexpect != ntotal) {
    const double ratio = nexpect * 1. / ntotal;
    cout << __func__ << ": See " << ntotal << " events but expect " << nexpect
         << " for process " << process.Data() << " --> scaling by " << ratio << endl;
    h->Scale(ratio);
  }
}

TH1* load_component_hist_from_dataset(const TString& dataset_key,
                                      const int selection,
                                      const TString& name,
                                      const int isys = -1,
                                      TString var_name = var_,
                                      const bool apply_norm_scaling = true) {
  auto info = get_dataset_info(dataset_key);
  if(info.name_ == "") {
    cout << __func__ << ": No dataset info found for key " << dataset_key.Data() << endl;
    return nullptr;
  }

  TFile* f = TFile::Open(Form("%sConvAna.%s.%s.m%i.%s",
                              hist_path_, hist_func_, info.name_.Data(), hist_mode_, file_type_.Data()),
                         "READ");
  if(!f) {
    cout << __func__ << ": Unable to open histogram file for " << dataset_key.Data() << endl;
    return nullptr;
  }

  const TString hist_name = (isys < 0)
    ? Form("%sHist/trk_%i/%s", dir_path_.Data(), selection, var_name.Data())
    : Form("%sHist/sys_%i/%s_%i", dir_path_.Data(), selection, var_name.Data(), isys);

  TH1* h = (TH1*) f->Get(hist_name.Data());
  if(!h) {
    cout << __func__ << ": Input histogram for selection " << selection
         << " not found in file " << f->GetName() << ": Hist name = " << hist_name.Data() << endl;
    f->Close();
    return nullptr;
  }

  h = (TH1*) h->Clone(Form("%s_%i_%s%s",
                           dataset_key.Data(), selection, name.Data(),
                           (isys < 0) ? "" : Form("_sys_%i", isys)));
  h->SetDirectory(0);

  const int rebin = bin_width_ / h->GetBinWidth(1) + 1.e-3;
  if(rebin > 1) h->Rebin(rebin);

  if(apply_norm_scaling) {
    if(info.ndigi_ <= 0) {
      cout << __func__ << ": No estimate for expected events for process " << dataset_key.Data() << endl;
    } else {
      apply_expected_event_scaling(h, f, dataset_key, info.ndigi_);
    }
  }

  f->Close();
  return h;
}

TH1* load_component_hist_from_sets(const TString& dataset_key,
                                   const std::vector<int>& selections,
                                   const TString& name,
                                   const int isys = -1,
                                   TString var_name = var_,
                                   const bool apply_norm_scaling = true) {
  if(selections.empty()) return nullptr;

  TH1* h_sum = nullptr;
  int loaded = 0;
  for(const int set : selections) {
    TH1* h = load_component_hist_from_dataset(dataset_key,
                          set,
                          Form("%s_set_%i", name.Data(), set),
                          isys,
                          var_name,
                          apply_norm_scaling);
    if(!h) continue;
    if(!h_sum) {
      h_sum = (TH1*) h->Clone(Form("%s_merged", name.Data()));
      h_sum->SetDirectory(0);
    } else {
      h_sum->Add(h);
    }
    ++loaded;
  }

  if(!h_sum) {
    cout << __func__ << ": Failed to load any set for " << dataset_key.Data() << endl;
    return nullptr;
  }

  if(loaded != (int) selections.size()) {
    cout << __func__ << ": Loaded " << loaded << " / " << selections.size()
         << " sets for " << dataset_key.Data() << endl;
  }
  return h_sum;
}

#endif
