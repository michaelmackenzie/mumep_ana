// Study how well the RMC e- background can be constrained in e+ data
#include "mumep_ana/analysis/datasets.C"

double get_norm(DatasetInfo_t& info, TFile* f) {
  TTree* t_norm = (TTree*) f->Get("Ana/data/Norm");
  const char* name = f->GetName();
  double norm = info.norm();
  if(!t_norm) {
    cout << "Normalization tree for input " << name << " not found!\n";
  } else {
    Long64_t nseen(0), ntotal(0);
    t_norm->SetBranchAddress("nseen", &nseen);
    for(Long64_t entry = 0; entry < t_norm->GetEntries(); ++entry) {
      t_norm->GetEntry(entry);
      ntotal += nseen;
    }
    if(ntotal == 0) cout << __func__ << ": No normalization contained in the normalization tree for file " << name << endl;
    else {
      const Long64_t nexpect = info.ndigi_;
      if(nexpect == 0) cout << __func__ << ": No estimate for the number of expected events for file " << name << endl;
      else {
        if(nexpect != ntotal) {
          const double ratio = nexpect * 1. / ntotal;
          cout << __func__ << ": See " << ntotal << " events but expect " << nexpect << " for file " << name
               << " --> scaling by " << ratio << endl;
          norm *= ratio;
        }
      }
    }
  }
  return norm;
}

int study(TString dir = "/exp/mu2e/data/projects/run1a/mumep_ana/histograms/",
          const int base_set = 75) {

  //------------------------------------------
  // Retrieve the inputs
  //------------------------------------------

  // Configure the dataset info setup
  use_evtana_ = true;
  set_evtana_defaults();
  init_dataset_info();

  // Get the dataset infos
  auto rmc_ext_info = get_dataset_info("rmc_ext_0n");
  auto rmc_int_info = get_dataset_info("rmc_int_0n");

  // Get the files
  TString rmc_ext_name = "ConvAna.cnv_ana." + rmc_ext_info.name_ + ".m1.root";
  TString rmc_int_name = "ConvAna.cnv_ana." + rmc_int_info.name_ + ".m1.root";
  TFile* f_ext = TFile::Open(dir + rmc_ext_name, "READ");
  TFile* f_int = TFile::Open(dir + rmc_int_name, "READ");
  if(!f_ext || !f_int) return 1;

  // Get the histograms
  TH1* h_ext_m = (TH1*) f_ext->Get(Form("Ana/Hist/trk_%i/p_2", base_set));
  TH1* h_ext_p = (TH1*) f_ext->Get(Form("Ana/Hist/trk_%i/p_2", base_set + 4000));
  TH1* h_int_m = (TH1*) f_int->Get(Form("Ana/Hist/trk_%i/p_2", base_set));
  TH1* h_int_p = (TH1*) f_int->Get(Form("Ana/Hist/trk_%i/p_2", base_set + 4000));

  if(!h_ext_m || !h_ext_p || !h_int_m || !h_int_p) {
    cout << "Unable to retrieve all histograms!\n";
    return 1;
  }

  // Detach the histograms from the files
  h_ext_m->SetDirectory(0);
  h_ext_p->SetDirectory(0);
  h_int_m->SetDirectory(0);
  h_int_p->SetDirectory(0);

  //------------------------------------------
  // Normalize inputs
  //------------------------------------------

  const double npot = npot_run1a_;
  const double norm_ext = npot*get_norm(rmc_ext_info, f_ext);
  const double norm_int = npot*get_norm(rmc_int_info, f_int);

  h_ext_m->Scale(norm_ext);
  h_ext_p->Scale(norm_ext);
  h_int_m->Scale(norm_int);
  h_int_p->Scale(norm_int);

  h_ext_m->Rebin(5);
  h_ext_p->Rebin(5);
  h_int_m->Rebin(5);
  h_int_p->Rebin(5);

  // Close the files
  f_ext->Close();
  f_int->Close();

  //------------------------------------------
  // Plot the inputs
  //------------------------------------------

  gSystem->Exec("mkdir -p figures");
  gStyle->SetOptStat(0);
  TCanvas c;

  h_ext_m->SetLineColor(kBlue);
  h_int_m->SetLineColor(kBlue);
  h_ext_p->SetLineColor(kRed );
  h_int_p->SetLineColor(kRed );
  h_int_m->SetLineStyle(kDashed);
  h_int_p->SetLineStyle(kDashed);
  h_ext_m->SetLineWidth(2);
  h_int_m->SetLineWidth(2);
  h_ext_p->SetLineWidth(2);
  h_int_p->SetLineWidth(2);

  h_ext_m->Draw("hist");
  h_ext_p->Draw("hist same");
  h_int_m->Draw("hist same");
  h_int_p->Draw("hist same");

  h_ext_m->GetXaxis()->SetRangeUser(95., 105.);
  h_ext_m->SetTitle(Form("RMC tracks by charge;p (MeV/c);Events / %.2g MeV/c",
                         h_ext_m->GetBinWidth(1)));

  TLegend leg(0.3, 0.7, 0.89, 0.89);
  leg.SetNColumns(2);
  leg.SetLineWidth(0);
  leg.SetLineColor(0);
  leg.AddEntry(h_ext_m, "Exteral electrons");
  leg.AddEntry(h_ext_p, "Exteral positrons");
  leg.AddEntry(h_int_m, "Internal electrons");
  leg.AddEntry(h_int_p, "Internal positrons");
  leg.Draw();

  c.SaveAs("figures/mom.png");
  c.SetLogy();
  c.SaveAs("figures/mom_log.png");

  return 0;
}
