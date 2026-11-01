// Optimize the event selection

#include "defaults.C"
#include "physics.C"
#include "datasets.C"

int debug_level_ = 0;

const bool use_asimov_  = false; // use Asimov significance instead of S/sqrt(B)
const bool use_combine_ = false;

// Stopping-target designs, as in combine/samples.py TARGETS: unmixed OnSpill CeEndpoint signal and DIO samples
// (dataset, N(generated), events in the full dataset) and the mu- stops per resampled MuBeamCat event
struct TargetInfo_t {
  TString name;
  TString signal; double signal_gen; Long64_t signal_events;
  TString dio;    double dio_gen;    Long64_t dio_events;
  double stops_per_draw;
};
const std::vector<TargetInfo_t> targets_ = {
  {"nominal"       , "cele0b0s5r0101", 1e7, 4247008, "dio00b0s5r0101", 8e6, 2995328, 1435092 * 1000. / 4e9},
  {"st37plus6up6dn", "cele0b0s5r0110", 1e7, 4002040, "dio00b0s5r0110", 8e6, 2764431, 5856806 / 1e7},
  {"st49hole1to5"  , "cele0b0s5r0120", 1e7, 3967395, "dio00b0s5r0120", 8e6, 2717026, 6527950 / 1e7},
  {"st49hole9to41" , "cele0b0s5r0130", 1e7, 3961096, "dio00b0s5r0130", 8e6, 2700249, 5826669 / 1e7},
};
const double mubeam_per_pot_ = 213816 / 1e8; // MuBeamCat events per POT
TString target_hist_path_ = (gSystem->Getenv("CONVANA_TARGET_HIST_PATH")) ? TString(gSystem->Getenv("CONVANA_TARGET_HIST_PATH")) + "/"
                                                                          : TString("/exp/mu2e/data/projects/run1a/mumep_ana/histograms_alt_target_2026-10-31/");


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

// Evaluates the continuous Asimov CLs value for a given mu, s, and b
double GetAsimovCLs(double mu, double s, double b) {
    if (mu <= 0.0) return 1.0;

    double mu_s_plus_b = mu * s + b;

    // 1. Calculate q_mu (Profile Likelihood Ratio for the tested hypothesis)
    // Avoid log(0) issues if background is theoretically tiny
    double q_mu = 0.0;
    if (b > 0) {
        q_mu = 2.0 * (mu_s_plus_b * std::log(mu_s_plus_b / b) - mu * s);
    }
    if (q_mu < 0.0) q_mu = 0.0; // Protection against floating-point underflow
    double sqrt_q_mu = std::sqrt(q_mu);

    // 2. Calculate q_mu_A (Variance proxy under the background-only Asimov assumption)
    double sqrt_q_mu_A = (mu * s) / std::sqrt(b);

    // 3. Convert test statistics to standard normal Gaussian p-values
    double one_minus_cl_sb = 1.0 - ROOT::Math::normal_cdf(sqrt_q_mu);
    double one_minus_cl_b  = ROOT::Math::normal_cdf(sqrt_q_mu_A - sqrt_q_mu);

    // 4. Return CLs ratio
    if (one_minus_cl_b <= 0.0) return 0.0;
    return one_minus_cl_sb / one_minus_cl_b;
}

double getAsimovExpUL(double s, double b) {
  if(b <= 0. || s <= 0.) return 0.;

  // CL target
  const double target_cls = 0.10;

  // Determine a range for the signal strength
  double mu_min = 0.0;
  double mu_max = 1.5*(2.30 + 1.645 * std::sqrt(b)) / s;
  while (true) {
    double cls = GetAsimovCLs(mu_max, s, b);
    if (cls < target_cls) {
      break; // mu_max is high enough
    }
    // mu_max too small --> increase the size
    mu_min = mu_max;
    mu_max *= 2.0;
  }

  // Binary search given a defined mu range
  const double tolerance = 1e-4; // Accuracy for mu
  while ((mu_max - mu_min) > tolerance) {
    double mu_mid = 0.5 * (mu_min + mu_max);
    double cls = GetAsimovCLs(mu_mid, s, b);

    if (cls > target_cls) {
      mu_min = mu_mid; // mu is too small (CLs is too high)
    } else {
      mu_max = mu_mid; // mu is too large (CLs is too low)
    }
  }
  return 0.5 * (mu_min + mu_max);
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
  // For kVeto, [min, max] is the range the veto window may grow into, starting from an empty window at seed
  Var_t(TString name, double min, double max,
        int type = kBoth, double seed = std::nan("")) : min(min), max(max), lo(min), hi(max),
                                                        name(name), type(type) {
    if(type == kVeto) this->min = this->max = (std::isnan(seed)) ? 0.5*(min + max) : seed;
  }
  double min;
  double max;
  double lo; // nominal (search) range
  double hi;
  TString name;
  int type;
  bool discrete = false; // integer-valued: cut edges sit on half-integers and move in whole units
  bool applied = false;
};

// Integer-valued variable taking values in [min, max] (not for kVeto)
Var_t discrete_var(TString name, int min, int max, int type = kBoth) {
  Var_t var(name, min - 0.5, max + 0.5, type);
  var.discrete = true;
  return var;
}

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
    scale *= nexpect * 1. / ntotal; // account for missing events
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
  sample.tree->Draw("trk_p>>hist", cut, "goff");
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

// Widen a veto window [min, max] within [lo, hi] by one edge until the signal loss reaches the step
TestPoint_t find_veto_test_point(Var_t& var, std::vector<Var_t>& vars,
                                 std::vector<Sample_t>& samples, double step) {
  const double min_0 = var.min;
  const double max_0 = var.max;
  TString cut_string = build_cut_string(vars);
  const auto [n_sig_0, n_bkg_0] = evaluate_yields(samples, cut_string);
  const bool applied = var.applied;
  const double tol = 0.1; // tolerance on the efficiency cut
  var.applied = true;

  double sig[2] = {0., 0.}, edge[2] = {min_0, max_0}, n_sig[2] = {n_sig_0, n_sig_0}, n_bkg[2] = {n_bkg_0, n_bkg_0};
  for(int side = 0; side < 2; ++side) { // side 0: lower edge, side 1: upper edge
    double edge_lo = (side == 0) ? var.lo : max_0; // range of the moving edge
    double edge_hi = (side == 0) ? min_0  : var.hi;
    double loss = 0.;
    if(debug_level_ > 1) cout << "    Testing veto side " << side << endl;
    while(std::fabs(loss - step) > tol*step && (edge_hi - edge_lo) > (var.hi - var.lo)*1.e-4) {
      const double point = (edge_lo + edge_hi) / 2.;
      var.min = (side == 0) ? point : min_0;
      var.max = (side == 1) ? point : max_0;
      cut_string = build_cut_string(vars);
      const auto [n_sig_t, n_bkg_t] = evaluate_yields(samples, cut_string);
      loss = 1. - n_sig_t / n_sig_0;
      if(debug_level_ > 2) cout << "      Test veto [" << var.min << ", " << var.max << "] loss = "
                                << loss << endl;
      sig[side] = (n_bkg_t > 0.) ? get_sig(n_sig_t, n_bkg_t) : 0.;
      edge[side] = point;
      n_sig[side] = n_sig_t; n_bkg[side] = n_bkg_t;
      const bool widen = loss < step;
      if(side == 0) { if(widen) edge_hi = point; else edge_lo = point; }
      else          { if(widen) edge_lo = point; else edge_hi = point; }
    }
  }

  // Restore the nominal values
  var.min = min_0;
  var.max = max_0;
  var.applied = applied;

  // Return the best point
  if(sig[0] > sig[1]) return TestPoint_t(edge[0], max_0, n_sig_0, n_bkg_0, n_sig[0], n_bkg[0]);
  return                     TestPoint_t(min_0, edge[1], n_sig_0, n_bkg_0, n_sig[1], n_bkg[1]);
}

// Move an edge of a discrete variable inwards by whole units: the largest move with a signal loss within the
// step, or the smallest move that loses signal if a single unit already exceeds it (loss jumps between values)
TestPoint_t find_discrete_test_point(Var_t& var, std::vector<Var_t>& vars,
                                     std::vector<Sample_t>& samples, double step) {
  const double min_0 = var.min;
  const double max_0 = var.max;
  TString cut_string = build_cut_string(vars);
  const auto [n_sig_0, n_bkg_0] = evaluate_yields(samples, cut_string);
  const bool applied = var.applied;
  var.applied = true;

  double sig[2] = {0., 0.}, edge[2] = {min_0, max_0}, n_sig[2] = {n_sig_0, n_sig_0}, n_bkg[2] = {n_bkg_0, n_bkg_0};
  const int n_max = std::lround(max_0 - min_0) - 1; // largest move that keeps one value
  for(int side = 0; side < 2; ++side) { // side 0: lower edge, side 1: upper edge
    if(side == 0 && var.type == kRight) continue;
    if(side == 1 && var.type == kLeft ) continue;
    if(n_max < 1) continue;
    auto evaluate = [&](int n) { // yields with the edge moved inwards by n units
      var.min = (side == 0) ? min_0 + n : min_0;
      var.max = (side == 1) ? max_0 - n : max_0;
      return evaluate_yields(samples, build_cut_string(vars));
    };
    // Binary search for the largest move with loss <= step, using that the loss grows with the move
    int lo = 0, hi = n_max;
    if(1. - evaluate(n_max).first / n_sig_0 <= step) lo = n_max;
    while(hi - lo > 1) {
      const int mid = (lo + hi) / 2;
      if(1. - evaluate(mid).first / n_sig_0 <= step) lo = mid;
      else                                           hi = mid;
    }
    int n = lo;
    std::pair<double,double> yields = evaluate(n);
    if(yields.first >= n_sig_0 && n < n_max) yields = evaluate(++n); // no signal lost yet: one more unit
    if(debug_level_ > 2) cout << "      Test discrete [" << var.min << ", " << var.max << "] loss = "
                              << 1. - yields.first / n_sig_0 << endl;
    sig[side] = (yields.second > 0.) ? get_sig(yields.first, yields.second) : 0.;
    edge[side] = (side == 0) ? var.min : var.max;
    n_sig[side] = yields.first; n_bkg[side] = yields.second;
  }

  // Restore the nominal values
  var.min = min_0;
  var.max = max_0;
  var.applied = applied;

  // Return the best point
  const bool left = var.type == kLeft || (var.type == kBoth && sig[0] > sig[1]);
  if(left) return TestPoint_t(edge[0], max_0, n_sig_0, n_bkg_0, n_sig[0], n_bkg[0]);
  return          TestPoint_t(min_0, edge[1], n_sig_0, n_bkg_0, n_sig[1], n_bkg[1]);
}

TestPoint_t find_test_point(Var_t& var, std::vector<Var_t>& vars,
                        std::vector<Sample_t>& samples, double step) {
  if(var.type == kVeto) return find_veto_test_point(var, vars, samples, step);
  if(var.discrete) return find_discrete_test_point(var, vars, samples, step);
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
      if(side == 0) {
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
      if(var.type == kVeto) { // random veto window within the search range
        var_test.min = var.lo + rnd.Uniform()*(var.hi - var.lo);
        var_test.max = var_test.min + rnd.Uniform()*(var.hi - var_test.min);
        vars_test.push_back(var_test);
        continue;
      }
      if(var.discrete) { // whole-unit edges, keeping at least one value
        if(var.type != kRight) var_test.min = var.min + rnd.Integer(std::lround(var.max - var.min));
        if(var.type != kLeft ) var_test.max = var.max - rnd.Integer(std::lround(var.max - var_test.min));
        vars_test.push_back(var_test);
        continue;
      }
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

// optimize_crv: drop the ConvAna CRV veto (add the vetoed set, base_set + 1000) and optimize the CRV veto window
// target: stopping-target design of targets_ (mu- --> e- only). The signal and DIO are the target's unmixed OnSpill
//         samples, normalized to its stops per POT (also applied to the RMC), using set 81 (see optimize_target_selection)
int optimize_selection(const bool mumem = false, const int base_set = 40, const bool optimize_crv = false,
                       const TString target = "", const int ngrid_tests = 1e4) {

  // Retrieve the relevant backgrounds
  const double pmin = (mumem) ? 101. : 90.;
  const double pmax = (mumem) ? 107. : 93.;
  init_physics("run1a"); // default to Run 1A
  use_evtana_ = true;
  set_evtana_defaults();
  init_dataset_info();
  std::vector<TString> names = {(mumem) ? "mumem" : "mumep", "cosmic", "rmc_ext_0n", "rmc_int_0n", "rpc_ext", "rpc_int"};

  // Stopping-target study: replace the signal and add the DIO with the target's samples
  const bool do_target = target != "";
  if(do_target) {
    if(!mumem) {
      cout << "The stopping-target study is only defined for mu- --> e-\n";
      return 1;
    }
    const TargetInfo_t* tgt = nullptr;
    for(const auto& t : targets_) if(t.name == target) tgt = &t;
    if(!tgt) {
      cout << "Target \"" << target << "\" not found!\n";
      return 1;
    }
    const double nmu_per_pot = mubeam_per_pot_ * tgt->stops_per_draw;
    const double rate_sig = muon_capture_fraction_*nmu_per_pot;
    const double rate_dio_95 = (1.-muon_capture_fraction_)*nmu_per_pot*dio_frac_95_;
    datasets_["mumem"] = DatasetInfo_t(tgt->signal_gen, tgt->signal_events, 0., 1., rate_sig   , tgt->signal);
    datasets_["dio"]   = DatasetInfo_t(tgt->dio_gen   , tgt->dio_events   , 0., 1., rate_dio_95, tgt->dio   );
    for(TString rmc : {"rmc_ext_0n", "rmc_int_0n"}) datasets_[rmc].theory_ *= nmu_per_pot / nmuons_per_pot_; // RMC follows the stops
    names.insert(names.begin() + 1, "dio");
    printf("Target %s: stops/POT = %.4g (Run 1A default %.4g)\n", target.Data(), nmu_per_pot, nmuons_per_pot_);
  }

  // Selection applied to all samples that isn't optimized
  TString preselection = Form("trk_p < %.2f && trk_p > %.2f && (%s * trk_charge) > 0 && weight > 0.",
                              pmax, pmin, (mumem) ? "-1" : "1");
  if(do_target) { // set 81 has no material, cluster, or CRV requirements; the ST intersection cuts are left out
    preselection += " && trk_cluster > 0. && trk_opa == 0 && trk_ntsda == 0";
    if(!optimize_crv) preselection += " && !(trk_crv_raw_dt > 0. && trk_crv_raw_dt < 150.)"; // closest CRV cluster
  }
  cout << "Preselection: " << preselection << endl;

  std::vector<int> sets = {base_set};
  if(optimize_crv && !do_target) sets.push_back(base_set + 1000); // ConvAna kCRVVetoOffset (set 81 is not CRV vetoed)
  std::vector<Sample_t> samples;
  for(auto name : names) {
    const auto& info = get_dataset_info(name);
    if(info.name_ == "") {
      cout << "Dataset \"" << name << "\" not found!\n";
      return 1;
    }
    const bool target_sample = do_target && (name == "mumem" || name == "dio");
    TFile* f = TFile::Open(Form("%sConvAna.%s.%s.m%i.%s", (target_sample) ? target_hist_path_.Data() : hist_path_,
                                hist_func_, info.name_.Data(), hist_mode_, file_type_.Data()), "READ");
    if(!f) {
      cout << "Dataset " << name << " hist file not found!\n";
      return 1;
    }
    for(int set : sets) { // each set is a separate sample, so each has its own sub-sampling scale
      const char* tpath = Form("%sHist/trs_%i/tree", dir_path_.Data(), set);
      TTree* t = (TTree*) f->Get(tpath);
      if(!t) {
        cout << "Dataset " << name << " tree not found! Path = " << tpath << endl;
        return 1;
      }
      gROOT->cd(); // put the copy in gROOT directory
      // apply momentum window, track charge, and event weight selection
      const int max_events = (name == "mumem" || name == "mumep") ? 100000 : 1e7; // reduce signal, comes out in the relative S/sqrt(B)
      const double events_scale = (t->GetEntries() > max_events) ? t->GetEntries() * 1./max_events : 1.;
      t = t->CopyTree(preselection, "", max_events);
      if(!t) {
        cout << "Dataset " << name << " set " << set << " has no events that pass!\n";
        continue;
      }
      const double scale = events_scale*evaluate_norm(name, f, info);
      samples.push_back(Sample_t(f, t, scale*info.theory_, name));
      samples.back().signal = name == "mumem" || name == "mumep";
      cout << name << " (set " << set << "): N(entries) = " << t->GetEntries() << " scale = " << scale
           << " yield = " << evaluate_yield(samples.back(), "1") << endl;
    }
  }

  // Define the list of variables to test
  std::vector<Var_t> vars;
  vars.push_back(Var_t("trk_cos"       ,  0.4,   0.9));
  vars.push_back(Var_t("trk_t0"        , 500., 1650.));
  vars.push_back(Var_t("trk_rmax"      , 400.,  700.));
  vars.push_back(Var_t("trk_ep"        ,   0.,   1.2));
  vars.push_back(Var_t("trk_dt"        ,  -4.,    4.));
  vars.push_back(Var_t("trk_fitcon"    ,   0.,    1., kLeft)); // only cut low track quality
  vars.push_back(Var_t("trk_qual"      ,  -1.,    1., kLeft)); // only cut on low scores for MVA ID scores
  vars.push_back(Var_t("trk_pid"       ,  -1.,    1., kLeft));
  if(!mumem) vars.push_back(Var_t("trk_trkonlypid",  -1.,    1., kLeft));
  // vars.push_back(Var_t("trk_csmid"     ,  -1.,    1., kLeft));
  if(do_target) { // remaining set 75/79 variables (trk_cos covers tan(dip)); no ST intersection or boundary cuts
    vars.push_back(Var_t("trk_d0"        , -150.,  150.));
    vars.push_back(Var_t("trk_t0err"     ,   0.,    2., kRight));
    vars.push_back(discrete_var("trk_nactive", 10, 60, kLeft));
    // vars.push_back(discrete_var("trk_opa"    ,  0,  1, kRight)); // OPA intersection flag
    // vars.push_back(discrete_var("trk_ntsda"  ,  0,  5, kRight)); // N(TSDA intersections)
  }
  if(optimize_crv) { // track - CRV cluster time vetoes, starting in the cosmic peaks (no cluster: -10000, never vetoed)
    vars.push_back(Var_t("trk_crv_dt"    , -150.,  150., kVeto, -10.)); // via the calo/ST, closest cluster
    vars.push_back(Var_t("trk_crv_raw_dt", -100.,  300., kVeto,  35.)); // raw cluster time, closest cluster
  }

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
  TString tag = (mumem) ? "mumem" : "mumep";
  if(do_target) tag += "_tgt_" + target;
  if(optimize_crv) tag += "_crv";
  TString f_out = "cut_strings_" + tag + ".txt";
  gSystem->Exec(Form("[ -f %s ] && rm %s; touch %s",
                     f_out.Data(), f_out.Data(), f_out.Data()));

  TGraph points, sigs;
  points.AddPoint(1., 0.); // first point is 100% efficiencies
  sigs.AddPoint(1.,1.); // first point is nominal significance

  // Continue running until low efficiency
  double eff = 1.;
  double max_sig = 0.;
  const double eff_step = (use_combine_) ? 0.03 : (mumem) ? 0.01 : 0.005; // aim for 1% efficiency steps
  const bool use_grad = true; // use gradient or just significance
  const double end_eff = (mumem) ? 0.5 : 0.8;
  while(eff > end_eff) {
    std::vector<TestPoint_t> candidates;
    int index = -1;
    double sig_best(0.), grad_best(-1.e10);
    for(auto& var : vars) {
      if(debug_level_ > 1) cout << "  Testing var " << var.name << endl;
      candidates.push_back(find_test_point(var, vars, samples, eff_step));
      if(candidates.back().n_sig >= candidates.back().n_sig_0) continue; // cut can't be tightened further
      if(index < 0 ||
         (use_grad && candidates.back().grad > grad_best) ||
         (!use_grad && candidates.back().sig > sig_best)) {
        sig_best = candidates.back().sig;
        grad_best = candidates.back().grad;
        index = candidates.size() - 1;
      }
    }
    if(index < 0) {
      cout << "No cut can be tightened further, stopping\n";
      break;
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
  auto [grid_points, grid_sigs] = perform_grid_scan(samples, vars_grid, n_bkg_0, n_sig_0, ngrid_tests,
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
  c.SaveAs("roc_" + tag + ".png");

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
  c.SaveAs("roc_grid_ " + tag + ".png");

  return 0;
}

// Stopping-target study: mu- --> e- selection without the ST intersection/boundary cuts, starting from set 81
// (trigger, upstream and multi-track vetoes, p > 100 MeV/c) with the target's unmixed signal and DIO samples
int optimize_target_selection(const TString target = "nominal", const bool optimize_crv = false) {
  return optimize_selection(true, 81, optimize_crv, target);
}
