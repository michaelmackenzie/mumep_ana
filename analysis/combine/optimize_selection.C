// Optimize the event selection

#include "../defaults.C"
#include "../physics.C"
#include "../datasets.C"

int debug_level_ = 0;

const bool use_asimov_  = false; // use Asimov significance instead of S/sqrt(B)
const bool use_combine_ = false;


double getCombine(double S, double B) {
  // Print a temporary card
  std::ofstream outfile;
  outfile.open("tmp_combine_card.txt");
  outfile << "imax 1 #number of bins\njmax * #number of processes\nkmax * #number of systematics\n\n";
  outfile << "bin cat\n";
  outfile << "observation " << Form("%i", int(S+B)) << "\n";
  outfile << "bin cat cat\n";
  outfile << "process signal background\n";
  outfile << "process 0 1\n";
  outfile << Form("rate %.6f %.6f\n", S, B);
  outfile << "sigN lnN 1.1 -\n";
  outfile << "bkgN lnN  -  1.2\n";
  outfile.close();

  TString command = "combine -d tmp_combine_card.txt -t -1";
  command += " --rMin -100. --rMax 100. --cl 0.9";
  command += " --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1";
  command += " --cminPreScan --cminPreFit 1 --rAbsAcc 0.0001 --rRelAcc 0.001";
  command += " | grep 'Observed Limit' | awk '{print $NF}'";
  TString output = gSystem->GetFromPipe(command);
  double ul = output.Atof();
  return (ul > 0.) ? 1./ul : 0.;
}

double getAsimovSig(double S, double B) {
  if (S <= 0.0) return 0.0;
  if (B <= 1e-9) return std::sqrt(2.0 * S); // Smooth mathematical limit at B = 0
  return std::sqrt(2.0 * ((S + B) * std::log(1.0 + S / B) - S));
}
double get_sig(double S, double B) {
  if(S <= 0.) return 0.;
  if(use_combine_) return getCombine(S,B);
  if(use_asimov_) return getAsimovSig(S,B);
  if(B <= 0.) return 0.;
  return S/std::sqrt(B);
}

struct Sample_t {
  Sample_t(TFile* f, TTree* t,
           double norm, TString name) : file(f), tree(t),
                                        norm(norm), name(name) {}
  TFile* file;
  TTree* tree;
  double norm;
  TString name;
  bool signal = false;
};

enum {kBoth, kLeft, kRight, kVeto};
struct Var_t {
  Var_t(TString name, double min, double max,
        int type = kBoth) : min(min), max(max),
                            name(name), type(type) {}
  double min;
  double max;
  TString name;
  int type;
  bool discrete = false;
  bool applied = false;
};

struct TestPoint_t {
  TestPoint_t(double min, double max,
              double n_sig_0, double n_bkg_0,
              double n_sig, double n_bkg) : min(min), max(max),
                                            n_sig_0(n_sig_0), n_bkg_0(n_bkg_0),
                                            n_sig(n_sig), n_bkg(n_bkg) {
    sig = get_sig(n_sig, n_bkg);
    grad = get_grad();
  }
  double get_grad(bool asimov_sig = false) { // d(sig) / d(signal)
    if(n_bkg_0 <= 0.) return 0.;
    const double dsig = sig - get_sig(n_sig_0, n_bkg_0);
    const double ds = n_sig_0 - n_sig;
    if(ds == 0.) return 0.;
    const double slope = dsig / ds;
    return slope;
  }
  double max;
  double min;
  double n_sig_0;
  double n_bkg_0;
  double n_sig;
  double n_bkg;
  double sig;
  double grad;
};

TString build_cut_string(std::vector<Var_t>& vars) {
  TString cut_string = "";
  bool first = true;
  for(auto& var : vars) {
    if(!var.applied) continue;
    if(!first) cut_string += " && ";
    else first = false;
    if     (var.type == kBoth) cut_string += Form("%s > %.6g && %s < %.6g",
                                                  var.name.Data(), var.min,
                                                  var.name.Data(), var.max);
    else if(var.type == kLeft) cut_string += Form("%s > %.6g",
                                                  var.name.Data(), var.min);
    else if(var.type == kRight)cut_string += Form("%s < %.6g",
                                                  var.name.Data(), var.max);
    else if(var.type == kVeto) cut_string += Form("(%s < %.6g || %s > %.6g)",
                                                  var.name.Data(), var.min,
                                                  var.name.Data(), var.max);
  }
  if(cut_string == "") return "(1)";
  return cut_string;
}

double evaluate_norm(TString name, TFile* f, const DatasetInfo_t& info) {
  double scale = (name.Contains("cosmic")) ? livetime_ : npot_;
  TTree* t_norm = (TTree*) f->Get(Form("%sdata/Norm", dir_path_.Data()));
  if(t_norm) {
    Long64_t nseen(0), ntotal(0);
    t_norm->SetBranchAddress("nseen", &nseen);
    for(Long64_t entry = 0; entry < t_norm->GetEntries(); ++entry) {
      t_norm->GetEntry(entry);
      ntotal += nseen;
    }
    const Long64_t nexpect = info.ndigi_;
    scale *= nexpect * 1. / nseen; // account for missing events
    if(name == "mumem" || name == "mumep") scale *= 1.e-14; // reasonable Rmue
  } else {
    cout << "Missing norm tree\n";
  }
  scale /= info.ngen_; // divide out the gen size

  return scale;
}
double evaluate_yield(Sample_t& sample, TString cut_string) {
  TCut cut("weight * (" + cut_string + ")");
  TH1F* hist = new TH1F("hist", "hist", 1, -1e20, 1e20);
  sample.tree->Draw("trkp>>hist", cut, "goff");
  double n = sample.norm*hist->Integral();
  delete hist;
  return n;
}

std::pair<double,double> evaluate_yields(std::vector<Sample_t>& samples, TString cut_string) {
  double n_sig = 0.;
  double n_bkg = 0.;
  for(auto& sample : samples) {
    const double yield = evaluate_yield(sample, cut_string);
    if(sample.signal) n_sig += yield;
    else              n_bkg += yield;
  }
  return std::pair<double,double>(n_sig, n_bkg);
}

TestPoint_t find_test_point(Var_t& var, std::vector<Var_t>& vars,
                        std::vector<Sample_t>& samples, double step) {
  const double min_0 = var.min;
  const double max_0 = var.max;
  TString cut_string = build_cut_string(vars);
  const auto [n_sig_0, n_bkg_0] = evaluate_yields(samples, cut_string);
  const double sig_0 = get_sig(n_sig_0, n_bkg_0);
  const bool applied = var.applied;

  double sig_l(0.), n_sig_l, n_bkg_l; // test sensitivities
  double sig_r(0.), n_sig_r, n_bkg_r;
  double min_l = min_0; // cut values for the optimized points
  double max_r = max_0;
  const double tol = 0.1; // tolerance on the efficiency cut
  var.applied = true;

  // Check each cut direction
  const bool is_veto = var.type == kVeto;
  for(int side = 0; side < 2; ++side) {
    if(side == 0 && var.type == kRight)  continue;
    if(side == 1 && var.type == kLeft )  continue;
    double min_t = min_0;
    double max_t = max_0;
    double loss = 0.;
    if(debug_level_ > 1) cout << "    Testing side " << side << endl;
    while(std::fabs(loss - step) > tol*step && (max_t - min_t) > (max_0 - min_0)*1.e-4) {
      double min_point = (min_t + max_t) / 2.;
      if(side == 0) { // left cut
        var.min = min_point;
        var.max = max_0;
      } else {        // right cut
        var.max = min_point;
        var.min = min_0;
      }
      cut_string = build_cut_string(vars);
      const auto [n_sig, n_bkg] = evaluate_yields(samples, cut_string);
      loss = 1. - n_sig / n_sig_0;
      if(debug_level_ > 2) cout << "      Test point [" << var.min << ", " << var.max << "] loss = "
           << loss << endl;
      if(side == is_veto) { // if vetoing, invert the assumptions
        sig_l = (n_bkg > 0.) ? get_sig(n_sig, n_bkg) : 0.;
        min_l = min_point;
        n_sig_l = n_sig; n_bkg_l = n_bkg;
        if(loss < step) min_t = min_point; // increase the cut
        else            max_t = min_point; // decrease the cut
      } else {
        sig_r = (n_bkg > 0.) ? get_sig(n_sig, n_bkg) : 0.;
        max_r = min_point;
        n_sig_r = n_sig; n_bkg_r = n_bkg;
        if(loss < step) max_t = min_point; // increase the cut
        else            min_t = min_point; // decrease the cut
      }
    }
  }

  // Restore the nominal values
  var.min = min_0;
  var.max = max_0;
  var.applied = applied;

  // Return the best point
  if(sig_l > sig_r) return TestPoint_t(min_l, max_0, n_sig_0, n_bkg_0, n_sig_l, n_bkg_l);
  return                   TestPoint_t(min_0, max_r, n_sig_0, n_bkg_0, n_sig_r, n_bkg_r);
}

std::pair<TGraph*,TGraph*> perform_grid_scan(std::vector<Sample_t>& samples,
                                             std::vector<Var_t>& vars,
                                             const double n_bkg_0, const double n_sig_0,
                                             const int ntests = 1e4,
                                             TGraph* scan_points = nullptr
                                             ) {

  TRandom3 rnd(90);
  TGraph* points = new TGraph();
  TGraph* sigs   = new TGraph();
  const double sig_0 = get_sig(n_sig_0, n_bkg_0);
  for(int itest = 0; itest < ntests; ++itest) {
    if(itest % 1000 == 0) cout << "Performing random grid point "
                               << itest << " / " << ntests
                               << endl;
    std::vector<Var_t> vars_test;
    for(const auto& var : vars) {
      auto var_test(var);
      var_test.applied = rnd.Uniform() > 0.5;
      if(var.type != kRight) var_test.min = var.min + rnd.Uniform()*(var.max - var.min);
      if(var.type != kLeft ) var_test.max = var.max - rnd.Uniform()*(var.max - var_test.min);
      vars_test.push_back(var_test);
    }
    TString cut_string = build_cut_string(vars_test);
    const auto [n_sig, n_bkg] = evaluate_yields(samples, cut_string);
    const double sig = (n_bkg > 0.) ? get_sig(n_sig, n_bkg) / sig_0 : 0.;
    const double eff = n_sig / n_sig_0;
    if(n_bkg > 0.) {
      const double rej = 1. - n_bkg / n_bkg_0;
      points->AddPoint(eff, rej);
      sigs  ->AddPoint(eff, sig);
      // Compare to the scan points
      if(scan_points && eff > scan_points->GetPointX(scan_points->GetN()-1)) {
        // int low_point  = 0;
        // int high_point = scan_points->GetN();
        // int point; double scan_eff, scan_rej;
        // while(low_point < high_point) {
        //   point = (high_point + low_point)/2;
        //   scan_eff = scan_points->GetX()[point];
        //   scan_rej = scan_points->GetY()[point];
        //   if(scan_eff < eff) high_point = point - 1;
        //   else               low_point  = point + 1;
        // }
        double scan_eff = eff;
        double scan_rej = scan_points->Eval(eff);
        if((rej <  0.5 && rej > scan_rej*1.01) || // 1% better
           (rej <= 0.5 && (1. - rej) < (1. - scan_rej)*0.99)) {
          cout << __func__ << ": Found a selection above the scan! (eff, rej) = ("
               << eff << ", " << rej << ") vs. scan  (" << scan_eff << ", " << scan_rej
               << ") cuts = " << cut_string << endl;
        }
      }
    }
  }
  return std::pair<TGraph*, TGraph*>(points,sigs);
}

int optimize_selection(const bool mumem = true, const int base_set = 70) {

  // Retrieve the relevant backgrounds
  const double pmin = (mumem) ? 101. : 90.;
  const double pmax = (mumem) ? 107. : 93.;
  init_physics("run1a"); // default to Run 1A
  init_dataset_info();
  std::vector<TString> names = {(mumem) ? "mumem" : "mumep", "cosmic", "rpc_ext", "rpc_int"};
  std::vector<Sample_t> samples;
  for(auto name : names) {
    const auto& info = get_dataset_info(name);
    if(info.name_ == "") {
      cout << "Dataset " << name << " not found!\n";
      return 1;
    }
    TFile* f = TFile::Open(Form("%sConvAna.%s.%s.m%i.hist", hist_path_, hist_func_,
                                info.name_.Data(), hist_mode_), "READ");
    if(!f) {
      cout << "Dataset " << name << " hist file not found!\n";
      return 1;
    }
    TTree* t = (TTree*) f->Get(Form("%sHist/trs_%i/Tree", dir_path_.Data(), base_set));
    if(!t) {
      cout << "Dataset " << name << " tree not found!\n";
      return 1;
    }
    gROOT->cd(); // put the copy in gROOT directory
    // apply momentum window, track charge, and event weight selection
    const int max_events = (name == "mumem" || "mumep") ? 40000 : 1e7; // reduce signal, comes out in the relative S/sqrt(B)
    const double events_scale = (t->GetEntries() > max_events) ? t->GetEntries() * 1./max_events : 1.;
    t = t->CopyTree(Form("trkp < %.2f && trkp > %.2f && (%s * trkq) > 0 && weight > 0.",
                         pmax, pmin, (mumem) ? "-1" : "1"), "", max_events);
    const double scale = events_scale*evaluate_norm(name, f, info);
    samples.push_back(Sample_t(f, t, scale*info.theory_, name));
    samples.back().signal = name == "mumem" || name == "mumep";
    cout << name << ": N(entries) = " << t->GetEntries() << " scale = " << scale << " yield = " << evaluate_yield(samples.back(), "1") << endl;
  }

  // Define the list of variables to test
  std::vector<Var_t> vars;
  vars.push_back(Var_t("trkcostheta",  0.4,   0.9));
  vars.push_back(Var_t("trkt0"      , 501., 1650.));
  vars.push_back(Var_t("trkrmax"    , 400.,  700.));
  vars.push_back(Var_t("trkep"      ,   0.,   1.2));
  vars.push_back(Var_t("trkdt"      ,  -4.,    4.));
  vars.push_back(Var_t("trkfitcon"  ,   0.,    1., kLeft)); // only cut low track quality
  vars.push_back(Var_t("trkqual"    ,  -1.,    1., kLeft)); // only cut on low scores for MVA ID scores
  vars.push_back(Var_t("trkpid"     ,  -1.,    1., kLeft));
  vars.push_back(Var_t("trkonlypid" ,  -1.,    1., kLeft));
  // vars.push_back(Var_t("trkcsmid"   ,  -1.,    1., kLeft));

  // Validate the variables
  for(auto& var : vars) {
    for(auto& sample : samples) {
      if(!sample.tree->GetBranch(var.name)) {
        cout << "Var " << var.name << " is undefined for sample " << sample.name << endl;
        return 1;
      }
    }
  }
  std::vector<Var_t> vars_grid(vars); // save the nominal bounds for a grid test

  // Start running the optimization
  const auto [n_sig_0, n_bkg_0] = evaluate_yields(samples, "1");
  const double sig_0 = get_sig(n_sig_0, n_bkg_0);
  printf("N(signal) = %.1f, N(background) = %.1f, significance = %.2f\n",
         n_sig_0, n_bkg_0, sig_0);

  // File to output results into
  TString f_out = "cut_strings.txt";
  gSystem->Exec(Form("[ -f %s ] && rm %s; touch %s",
                     f_out.Data(), f_out.Data(), f_out.Data()));

  TGraph points, sigs;
  points.AddPoint(1., 0.); // first point is 100% efficiencies
  sigs.AddPoint(1.,1.); // first point is nominal significance

  // Continue running until low efficiency
  double eff = 1.;
  double max_sig = 0.;
  const double eff_step = (use_combine_) ? 0.03 : 0.01; // aim for 1% efficiency steps
  const bool use_grad = true; // use gradient or just significance
  while(eff > 0.3) {
    std::vector<TestPoint_t> candidates;
    size_t index = 0;
    double sig_best(0.), grad_best(-1.e10);
    for(auto& var : vars) {
      if(debug_level_ > 1) cout << "  Testing var " << var.name << endl;
      candidates.push_back(find_test_point(var, vars, samples, eff_step));
      if((use_grad && candidates.back().grad > grad_best) ||
         (!use_grad && candidates.back().sig > sig_best)) {
        sig_best = candidates.back().sig;
        grad_best = candidates.back().grad;
        index = candidates.size() - 1;
      }
    }
    // Update with the best cut
    vars[index].min = candidates[index].min;
    vars[index].max = candidates[index].max;
    vars[index].applied = true;
    TString cut_string = build_cut_string(vars);
    const auto [n_sig, n_bkg] = evaluate_yields(samples, cut_string);
    const double sig = (n_bkg > 0.) ? get_sig(n_sig, n_bkg) / sig_0 : 0.;
    eff = n_sig / n_sig_0;
    printf("Eff = %.3f, significance = %.3f: %s\n",
           eff, sig, cut_string.Data());
    if(n_bkg > 0.) {
      points.AddPoint(eff, 1. - n_bkg / n_bkg_0);
      sigs  .AddPoint(eff, sig);
      if(sig > max_sig) max_sig = sig;
    }

    // Add the result to the file
    gSystem->Exec(Form("echo %.6f %.6g %.4f \"%s\" >> %s",
                       n_sig/n_sig_0, n_bkg/n_bkg_0, sig, cut_string.Data(), f_out.Data()));
    if(n_bkg <= 0.) break; // can't continue
  }

  // Do a random grid test
  const int ntests = 1e3;
  auto [grid_points, grid_sigs] = perform_grid_scan(samples, vars_grid, n_bkg_0, n_sig_0, ntests,
                                                    &points
                                                    );

  // re-scale the significances to be 0 - 1
  if(max_sig > 0.) {
    cout << "Re-scaling significances by 1 / " << max_sig << endl;
    for(int point = 0; point < sigs.GetN(); ++point)
      sigs.SetPointY(point, sigs.GetPointY(point)/max_sig);
    for(int point = 0; point < grid_sigs->GetN(); ++point)
      grid_sigs->SetPointY(point, grid_sigs->GetPointY(point)/max_sig);
  }

  // Plot the ROC curve
  TCanvas c;
  points.SetLineWidth(2);
  points.SetLineColor(kRed);
  points.Draw("AL");
  sigs  .SetLineWidth(2);
  sigs  .SetLineColor(kGreen);
  sigs  .Draw("L");
  points.GetXaxis()->SetRangeUser(0., 1.1);
  points.GetYaxis()->SetRangeUser(0., 1.1);
  c.SetGrid();
  TLegend leg;
  leg.AddEntry(&points, "Background rejection", "L");
  leg.AddEntry(&sigs  , "Significance", "L");
  leg.Draw();
  c.SaveAs("roc.png");

  gStyle->SetOptStat(0);
  TH1* haxis = new TH1F("haxis","",1,0.,1.1);
  haxis->Draw("hist");
  haxis->SetTitle(";Signal efficiency;Background rejection");
  grid_points->Draw("P");
  grid_points->SetMarkerColor(kBlack);
  grid_points->SetMarkerStyle(6);
  points.Draw("L");
  haxis->GetXaxis()->SetRangeUser(0., 1.1);
  haxis->GetYaxis()->SetRangeUser(0., 1.3);
  TLegend leg_p(0.15, 0.80, 0.89, 0.85);
  leg_p.SetNColumns(2);
  leg_p.AddEntry(&points    , "Greedy selection", "L");
  leg_p.AddEntry(grid_points, "Random selection", "P");
  leg_p.Draw();
  c.SaveAs("roc_grid.png");

  return 0;
}
