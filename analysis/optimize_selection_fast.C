// Optimize the event selection
// Faster version of optimize_selection.C: the samples are loaded into memory once, and each cut edge is found from
// cumulative yields over the events sorted by the variable instead of bisecting with TTree::Draw. Additionally:
// - figure of merit S/sqrt(B) (default), with 1 / the expected 90% CL CLs upper limit (robust for small backgrounds)
//   also reported and plotted for comparison, or used as the figure of merit with fom_ = kExpectedUL
// - mu- --> e-: the yields are in the 102.5 - 105 MeV/c signal region, with the cosmic ray and RPC backgrounds taken
//   from 100 - 110 MeV/c as flat in momentum to reduce their fluctuations
// - the selection is optimized on half the events (even event numbers) and evaluated on the other half
// - besides the greedy path, the selection is refined by a local search from several starting points
// Don't load both macros in one session (shared names).

#pragma cling optimize(3)

#include "defaults.C"
#include "physics.C"
#include "datasets.C"
#include <array>
#include <numeric>

int debug_level_ = 0;

// Figure of merit
enum {kSOverSqrtB, kAsimovZ, kExpectedUL};
int fom_ = kSOverSqrtB; // kExpectedUL: 1 / N(signal) at the expected 90% CL CLs upper limit, so S / N_UL(B)
const bool use_combine_ = false;

// Momentum windows for mu- --> e-: the signal region, and a wider window for the backgrounds approximated as flat in
// momentum, scaled to the signal region width
double p_sr_min_   = 102.5;
double p_sr_max_   = 105.;
double p_flat_min_ = 100.;
double p_flat_max_ = 110.;
bool is_flat_bkg(const TString name) { return name.BeginsWith("cosmic") || name.BeginsWith("rpc"); }

// Optimize on the events with even event numbers and evaluate on the odd ones (each half weighted by 2)
bool split_ = true;

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

double getAsimovSig(double S, double B) {
  if (S <= 0.0) return 0.0;
  if (B <= 1e-9) return std::sqrt(2.0 * S); // Smooth mathematical limit at B = 0
  return std::sqrt(2.0 * ((S + B) * std::log(1.0 + S / B) - S));
}

// 90% CL CLs upper limit on the signal events for n observed events with a known background b
double poisson_cls_ul(const int n, const double b, const double cl = 0.9) {
  const double alpha = 1. - cl;
  const double cl_b = ROOT::Math::poisson_cdf(n, b);
  auto cls = [&](double s) { return ROOT::Math::poisson_cdf(n, s + b) / cl_b; };
  double lo = 0., hi = n + 3. + 5.*std::sqrt(n + 1.);
  while(cls(hi) > alpha) hi *= 2.;
  for(int i = 0; i < 100 && hi - lo > 1.e-7*hi; ++i) {
    const double mid = 0.5*(lo + hi);
    if(cls(mid) > alpha) lo = mid;
    else                 hi = mid;
  }
  return 0.5*(lo + hi);
}

// Expected upper limit on the signal events: the average over background-only outcomes, sum_n P(n|b) N_UL(n, b)
double expected_ul_events_exact(const double b) {
  if(b <= 0.) return poisson_cls_ul(0, 0.);
  const int n_lo = std::max(0, (int) std::floor(b - 7.*std::sqrt(b) - 2.));
  const int n_hi = (int) std::ceil(b + 7.*std::sqrt(b) + 5.);
  double ul(0.), norm(0.);
  for(int n = n_lo; n <= n_hi; ++n) {
    const double p = ROOT::Math::poisson_pdf(n, b);
    ul   += p*poisson_cls_ul(n, b);
    norm += p;
  }
  return ul / norm;
}

// Tabulated in log10(b) for speed, computed directly above the table
double expected_ul_events(const double b) {
  static std::vector<double> table;
  const double lmin(-4.), lmax(3.);
  const int npoints = 701;
  if(table.empty()) {
    for(int i = 0; i < npoints; ++i) table.push_back(expected_ul_events_exact(std::pow(10., lmin + i*(lmax - lmin)/(npoints - 1))));
  }
  if(b <= std::pow(10., lmin)) return table.front();
  const double l = std::log10(b);
  if(l >= lmax) return expected_ul_events_exact(b);
  const double x = (l - lmin)/(lmax - lmin)*(npoints - 1);
  const int i = (int) x;
  const double f = x - i;
  return table[i]*(1. - f) + table[i+1]*f;
}

// 1 / expected 90% CL upper limit on the signal strength: S / N_UL(B), independent of the signal normalization
double get_inv_limit(double S, double B) {
  if(S <= 0.) return 0.;
  return S / expected_ul_events(std::max(B, 0.));
}

double get_sig(double S, double B) {
  if(S <= 0.) return 0.;
  if(use_combine_) return getCombine(S,B);
  if(fom_ == kExpectedUL) return get_inv_limit(S, B);
  if(fom_ == kAsimovZ) return getAsimovSig(S,B);
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
  double get_grad() { // d(sig) / d(signal)
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

// %.8g so the written cut values resolve the event values the cuts were placed between
TString build_cut_string(const std::vector<Var_t>& vars) {
  TString cut_string = "";
  bool first = true;
  for(auto& var : vars) {
    if(!var.applied) continue;
    if(!first) cut_string += " && ";
    else first = false;
    if     (var.type == kBoth) cut_string += Form("%s > %.8g && %s < %.8g",
                                                  var.name.Data(), var.min,
                                                  var.name.Data(), var.max);
    else if(var.type == kLeft) cut_string += Form("%s > %.8g",
                                                  var.name.Data(), var.min);
    else if(var.type == kRight)cut_string += Form("%s < %.8g",
                                                  var.name.Data(), var.max);
    else if(var.type == kVeto) cut_string += Form("(%s < %.8g || %s > %.8g)",
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

// Yield of a cut string from the sample's tree (slow, used for printouts and the cross-check)
double draw_yield(Sample_t& sample, TString cut_string) {
  TCut cut("weight * (" + cut_string + ")");
  TH1F* hist = new TH1F("hist", "hist", 1, -1e20, 1e20);
  sample.tree->Draw("trk_p>>hist", cut, "goff");
  double n = sample.norm*hist->Integral();
  delete hist;
  return n;
}

std::pair<double,double> draw_yields(std::vector<Sample_t>& samples, TString cut_string) {
  double n_sig = 0.;
  double n_bkg = 0.;
  for(auto& sample : samples) {
    const double yield = draw_yield(sample, cut_string);
    if(sample.signal) n_sig += yield;
    else              n_bkg += yield;
  }
  return std::pair<double,double>(n_sig, n_bkg);
}

//------------------------------------------------------------------------------------
// In-memory events of all samples: the variables, the normalized signal and background weights, and for each
// variable the events ordered by value
struct Events_t {
  std::vector<std::vector<double>> vals; // [variable][event]
  std::vector<double> w_sig;
  std::vector<double> w_bkg;
  std::vector<std::vector<int>> order; // [variable] event indices sorted by value
  int size() const { return w_sig.size(); }
};

void sort_events(Events_t& events) {
  events.order.assign(events.vals.size(), std::vector<int>(events.size()));
  for(size_t ivar = 0; ivar < events.vals.size(); ++ivar) {
    auto& order = events.order[ivar];
    const auto& vals = events.vals[ivar];
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return vals[a] < vals[b]; });
  }
}

bool load_events(std::vector<Sample_t>& samples, const std::vector<Var_t>& vars, Events_t& events,
                 std::vector<double>& event_numbers) {
  events.vals.assign(vars.size(), std::vector<double>());
  for(auto& sample : samples) {
    TTree* t = sample.tree;
    const Long64_t n = t->GetEntries();
    t->SetEstimate(n + 1);
    auto column = [&](TString expr, std::vector<double>& out) {
      if(t->Draw(expr, "", "goff") != n) {
        cout << "Failed to read " << expr << " for sample " << sample.name << endl;
        return false;
      }
      out.insert(out.end(), t->GetV1(), t->GetV1() + n);
      return true;
    };
    std::vector<double> weight;
    if(!column("weight", weight)) return false;
    if(!column("event", event_numbers)) return false;
    for(double w : weight) {
      events.w_sig.push_back((sample.signal) ? w*sample.norm : 0.);
      events.w_bkg.push_back((sample.signal) ? 0. : w*sample.norm);
    }
    for(size_t ivar = 0; ivar < vars.size(); ++ivar) {
      if(!column(vars[ivar].name, events.vals[ivar])) return false;
    }
  }
  sort_events(events);
  return true;
}

// Split the events into halves by event number parity, each weighted by 2 to estimate the full sample yields
void split_events(const Events_t& all, const std::vector<double>& event_numbers, Events_t& train, Events_t& test) {
  train.vals.assign(all.vals.size(), std::vector<double>());
  test .vals.assign(all.vals.size(), std::vector<double>());
  for(int i = 0; i < all.size(); ++i) {
    Events_t& half = (std::fmod(event_numbers[i], 2.) == 0.) ? train : test;
    half.w_sig.push_back(2.*all.w_sig[i]);
    half.w_bkg.push_back(2.*all.w_bkg[i]);
    for(size_t ivar = 0; ivar < all.vals.size(); ++ivar) half.vals[ivar].push_back(all.vals[ivar][i]);
  }
  sort_events(train);
  sort_events(test);
}

bool passes(const Var_t& var, const double v) {
  if(!var.applied) return true;
  if(var.type == kBoth)  return v > var.min && v < var.max;
  if(var.type == kLeft)  return v > var.min;
  if(var.type == kRight) return v < var.max;
  return v < var.min || v > var.max; // kVeto
}

// Signal and background yields passing the cuts, optionally storing which events pass
std::pair<double,double> evaluate_yields(const std::vector<Var_t>& vars, const Events_t& events,
                                         std::vector<char>* pass = nullptr) {
  double n_sig = 0.;
  double n_bkg = 0.;
  const int n = events.size();
  if(pass) pass->assign(n, 0);
  for(int i = 0; i < n; ++i) {
    bool ok = true;
    for(size_t ivar = 0; ok && ivar < vars.size(); ++ivar) ok = passes(vars[ivar], events.vals[ivar][i]);
    if(!ok) continue;
    if(pass) (*pass)[i] = 1;
    n_sig += events.w_sig[i];
    n_bkg += events.w_bkg[i];
  }
  return std::pair<double,double>(n_sig, n_bkg);
}

// Candidate positions y of one cut edge that keeps events with u > y, moving up from y0 towards y_lim. u is the
// variable oriented so the edge tightens upwards (u = -value for an upper edge). entries: (u, signal weight,
// background weight) of the events the edge acts on, sorted by u; fixed: yields kept wherever the edge is
struct EdgePoint_t {
  double y;
  double n_sig;
  double n_bkg;
};

std::vector<EdgePoint_t> scan_edge(const std::vector<std::array<double,3>>& entries, const double y0, const double y_lim,
                                   const bool discrete, const double s_fixed, const double b_fixed) {
  const int n = entries.size();
  std::vector<double> s_above(n + 1, 0.), b_above(n + 1, 0.); // yields of entries [i, n)
  for(int i = n - 1; i >= 0; --i) {
    s_above[i] = s_above[i+1] + entries[i][1];
    b_above[i] = b_above[i+1] + entries[i][2];
  }
  std::vector<EdgePoint_t> scan;
  auto add = [&](double y, int first) { scan.push_back({y, s_fixed + s_above[first], b_fixed + b_above[first]}); };

  if(discrete) { // whole-unit moves, keeping at least one value
    const int n_max = std::lround(y_lim - y0) - 1;
    for(int move = 1; move <= n_max; ++move) {
      const double y = y0 + move;
      const int first = std::upper_bound(entries.begin(), entries.end(), y,
                                         [](double val, const std::array<double,3>& e) { return val < e[0]; })
                        - entries.begin();
      add(y, first);
    }
    return scan;
  }

  // Continuous: just below each signal value, so all the background beneath it is removed for the same signal loss
  for(int i = 0; i < n; ) {
    int j = i;
    bool has_signal = false;
    while(j < n && entries[j][0] == entries[i][0]) has_signal |= entries[j++][1] > 0.;
    const double u = entries[i][0];
    if(has_signal && u > y0) {
      const double y = std::max(y0, (i > 0) ? 0.5*(entries[i-1][0] + u) : y0);
      if(y >= y_lim) break;
      add(y, i);
    }
    i = j;
  }
  return scan;
}

// The largest move with a signal loss within the step, or the first move if all exceed it. A move that loses no
// signal is taken one further, since the loss jumps between points
int choose_edge(const std::vector<EdgePoint_t>& scan, const double n_sig_0, const double step) {
  if(scan.empty()) return -1;
  int index = 0;
  for(int i = 0; i < (int) scan.size() && 1. - scan[i].n_sig / n_sig_0 <= step; ++i) index = i;
  if(scan[index].n_sig >= n_sig_0 && index + 1 < (int) scan.size()) ++index;
  return index;
}

// Tighten variable ivar by about one efficiency step, from the events passing the current cuts (pass). Each side's
// edge moves inwards (kVeto: the window widens outwards) and the side with the larger significance is returned
TestPoint_t find_test_point(const int ivar, const std::vector<Var_t>& vars, const Events_t& events,
                            const std::vector<char>& pass, const double n_sig_0, const double n_bkg_0,
                            const double step) {
  const Var_t& var = vars[ivar];
  const auto& vals = events.vals[ivar];
  const auto& order = events.order[ivar];
  const int n = order.size();
  TestPoint_t best(var.min, var.max, n_sig_0, n_bkg_0, n_sig_0, n_bkg_0); // no change
  bool found = false;
  for(int side = 0; side < 2; ++side) { // side 0: lower edge, side 1: upper edge
    if(side == 0 && var.type == kRight) continue;
    if(side == 1 && var.type == kLeft ) continue;
    const bool veto = var.type == kVeto;
    // Orientation: a cut edge tightens upwards for the lower edge, a veto edge for the upper edge
    const double sign = (veto) ? ((side == 0) ? -1. : 1.) : ((side == 0) ? 1. : -1.);
    double y0, y_lim;
    if(veto) { y0 = (side == 0) ? -var.min : var.max; y_lim = (side == 0) ? -var.lo  : var.hi; }
    else     { y0 = (side == 0) ?  var.min : -var.max; y_lim = (side == 0) ?  var.max : -var.min; }

    std::vector<std::array<double,3>> entries;
    entries.reserve(n);
    double s_fixed(0.), b_fixed(0.);
    for(int k = 0; k < n; ++k) {
      const int i = order[(sign > 0.) ? k : n - 1 - k]; // ascending u
      if(!pass[i]) continue;
      const double v = vals[i];
      bool movable(false), fixed(false);
      if(veto) { // events beside the window on this side move, those on the other side stay
        movable = (side == 0) ? v < var.min : v > var.max;
        fixed   = (side == 0) ? v > var.max : v < var.min;
      } else {   // the other edge of a two-sided cut stays where it is
        movable = var.type != kBoth || ((side == 0) ? v < var.max : v > var.min);
      }
      if(movable)    entries.push_back({sign*v, events.w_sig[i], events.w_bkg[i]});
      else if(fixed) { s_fixed += events.w_sig[i]; b_fixed += events.w_bkg[i]; }
    }

    const auto scan = scan_edge(entries, y0, y_lim, var.discrete && !veto, s_fixed, b_fixed);
    const int index = choose_edge(scan, n_sig_0, step);
    if(index < 0) continue;
    const auto& point = scan[index];
    const double edge = sign*point.y;
    TestPoint_t test((side == 0) ? edge : var.min, (side == 1) ? edge : var.max,
                     n_sig_0, n_bkg_0, point.n_sig, point.n_bkg);
    if(debug_level_ > 2) cout << "      Test side " << side << " [" << test.min << ", " << test.max << "] loss = "
                              << 1. - test.n_sig / n_sig_0 << " sig = " << test.sig << endl;
    if(!found || test.sig > best.sig) best = test;
    found = true;
  }
  return best;
}

// Loosen an applied cut (not a veto) until about a step of signal is gained, from the events passing the other cuts.
// Returns false if it can't be loosened
bool find_loose_point(const int ivar, const std::vector<Var_t>& vars, const Events_t& events, const double step,
                      TestPoint_t& result) {
  const Var_t& var = vars[ivar];
  if(!var.applied || var.type == kVeto) return false;
  const auto [n_sig_c, n_bkg_c] = evaluate_yields(vars, events);
  std::vector<Var_t> others(vars);
  others[ivar].applied = false;
  std::vector<char> pass;
  evaluate_yields(others, events, &pass);
  const auto& vals = events.vals[ivar];
  const auto& order = events.order[ivar];
  const int n = order.size();
  const double target = n_sig_c*(1. + step);
  bool found = false;
  for(int side = 0; side < 2; ++side) { // side 0: lower edge, side 1: upper edge
    if(side == 0 && var.type == kRight) continue;
    if(side == 1 && var.type == kLeft ) continue;
    const double sign  = (side == 0) ? 1. : -1.;
    const double y0    = (side == 0) ? var.lo  : -var.hi;  // loosest edge
    const double y_cur = (side == 0) ? var.min : -var.max; // current edge
    if(y_cur <= y0) continue;
    std::vector<std::array<double,3>> entries;
    entries.reserve(n);
    double s_loose(0.), b_loose(0.);
    for(int k = 0; k < n; ++k) {
      const int i = order[(sign > 0.) ? k : n - 1 - k]; // ascending u
      if(!pass[i]) continue;
      const double v = vals[i];
      if(var.type == kBoth && !((side == 0) ? v < var.max : v > var.min)) continue; // other edge stays
      entries.push_back({sign*v, events.w_sig[i], events.w_bkg[i]});
      if(sign*v > y0) { s_loose += events.w_sig[i]; b_loose += events.w_bkg[i]; }
    }
    // Edge positions from the loosest to the current one; take the tightest that gains the step
    auto scan = scan_edge(entries, y0, y_cur, var.discrete, 0., 0.);
    scan.insert(scan.begin(), EdgePoint_t{y0, s_loose, b_loose});
    int index = 0;
    for(int i = 0; i < (int) scan.size() && scan[i].n_sig >= target; ++i) index = i;
    if(scan[index].n_sig <= n_sig_c) continue; // no gain
    const double edge = sign*scan[index].y;
    TestPoint_t test((side == 0) ? edge : var.min, (side == 1) ? edge : var.max,
                     n_sig_c, n_bkg_c, scan[index].n_sig, scan[index].n_bkg);
    if(!found || test.sig > result.sig) result = test;
    found = true;
  }
  return found;
}

//------------------------------------------------------------------------------------
// Selections and the search strategies

struct Selection_t {
  std::vector<Var_t> vars;
  double n_sig;
  double n_bkg;
};

void apply_point(std::vector<Var_t>& vars, const int ivar, const TestPoint_t& point) {
  vars[ivar].min = point.min;
  vars[ivar].max = point.max;
  vars[ivar].applied = true;
}

// Tightening test points of all variables from the current selection, skipping those that can't be tightened
std::vector<std::pair<int,TestPoint_t>> tighten_candidates(const std::vector<Var_t>& vars, const Events_t& events,
                                                           const double step) {
  std::vector<char> pass;
  const auto [n_sig, n_bkg] = evaluate_yields(vars, events, &pass);
  std::vector<std::pair<int,TestPoint_t>> candidates;
  for(size_t ivar = 0; ivar < vars.size(); ++ivar) {
    if(debug_level_ > 1) cout << "  Testing var " << vars[ivar].name << endl;
    const auto point = find_test_point(ivar, vars, events, pass, n_sig, n_bkg, step);
    if(point.n_sig >= point.n_sig_0) continue; // cut can't be tightened further
    candidates.push_back({(int) ivar, point});
  }
  return candidates;
}

// The candidate with the largest gain per signal lost, or a random one of the top_k with rnd
int choose_greedy(std::vector<std::pair<int,TestPoint_t>>& candidates, TRandom3* rnd = nullptr, const int top_k = 1) {
  if(candidates.empty()) return -1;
  std::stable_sort(candidates.begin(), candidates.end(),
                   [](const std::pair<int,TestPoint_t>& a, const std::pair<int,TestPoint_t>& b) { return a.second.grad > b.second.grad; });
  if(!rnd) return 0;
  return rnd->Integer(std::min<int>(top_k, candidates.size()));
}

// Greedy path from a starting selection until the signal efficiency (relative to n_sig_ref) falls below end_eff
std::vector<Selection_t> run_greedy(std::vector<Var_t> vars, const Events_t& events, const double n_sig_ref,
                                    const double step, const double end_eff,
                                    TRandom3* rnd = nullptr, const int top_k = 1) {
  std::vector<Selection_t> path;
  double eff = 1.;
  while(eff > end_eff) {
    auto candidates = tighten_candidates(vars, events, step);
    const int index = choose_greedy(candidates, rnd, top_k);
    if(index < 0) {
      if(!rnd) cout << "No cut can be tightened further, stopping\n";
      break;
    }
    apply_point(vars, candidates[index].first, candidates[index].second);
    path.push_back({vars, candidates[index].second.n_sig, candidates[index].second.n_bkg});
    eff = path.back().n_sig / n_sig_ref;
    if(path.back().n_bkg <= 0.) break; // can't continue
  }
  return path;
}

Selection_t best_on_path(const std::vector<Selection_t>& path) {
  int best = 0;
  for(int i = 1; i < (int) path.size(); ++i) {
    if(get_sig(path[i].n_sig, path[i].n_bkg) > get_sig(path[best].n_sig, path[best].n_bkg)) best = i;
  }
  return path[best];
}

// Local search on the significance: each iteration tries tightening any cut by a step, loosening any cut by a step,
// and loosening one cut then tightening the greedy choice, and takes the largest improvement
Selection_t refine(std::vector<Var_t> vars, const Events_t& events, const double step, const int max_iter = 100) {
  auto [n_sig, n_bkg] = evaluate_yields(vars, events);
  for(int iter = 0; iter < max_iter; ++iter) {
    double best_sig = get_sig(n_sig, n_bkg)*(1. + 1.e-6);
    bool improved = false;
    std::vector<Var_t> best_vars;
    double best_s(0.), best_b(0.);
    auto consider = [&](const std::vector<Var_t>& test, double s, double b) {
      const double sig = get_sig(s, b);
      if(sig <= best_sig) return;
      best_sig = sig; best_vars = test; best_s = s; best_b = b; improved = true;
    };
    for(auto& [ivar, point] : tighten_candidates(vars, events, step)) {
      auto test(vars);
      apply_point(test, ivar, point);
      consider(test, point.n_sig, point.n_bkg);
    }
    for(size_t ivar = 0; ivar < vars.size(); ++ivar) {
      TestPoint_t loose(0., 0., 0., 0., 0., 0.);
      if(!find_loose_point(ivar, vars, events, step, loose)) continue;
      auto test(vars);
      apply_point(test, ivar, loose);
      consider(test, loose.n_sig, loose.n_bkg);
      auto candidates = tighten_candidates(test, events, step);
      const int index = choose_greedy(candidates);
      if(index < 0) continue;
      apply_point(test, candidates[index].first, candidates[index].second);
      consider(test, candidates[index].second.n_sig, candidates[index].second.n_bkg);
    }
    if(!improved) break;
    vars = best_vars; n_sig = best_s; n_bkg = best_b;
    if(debug_level_ > 0) printf("  Refine iteration %i: significance = %.4g\n", iter, get_sig(n_sig, n_bkg));
  }
  return {vars, n_sig, n_bkg};
}

// Starting selections from the existing mu- --> e- cut-sets without their ST intersection and boundary cuts
// (tan(dip) as cos(theta) = tan(dip) / sqrt(1 + tan(dip)^2), N(active) >= N as > N - 0.5)
struct SeedCut_t {
  TString name;
  double min;
  double max;
};
double tandip_to_cos(const double t) { return t / std::sqrt(1. + t*t); }
std::vector<std::pair<TString,std::vector<SeedCut_t>>> seed_selections() {
  return {
    {"set79", {{"trk_qual", 0.145, 1.}, {"trk_pid", 0.574, 1.}, {"trk_t0err", 0., 0.81},
               {"trk_cos", tandip_to_cos(0.56), tandip_to_cos(0.8787)}, {"trk_d0", -60., 94.},
               {"trk_nactive", 16.5, 60.5}, {"trk_t0", 540., 1650.}}},
    {"set75", {{"trk_qual", 0.155, 1.}, {"trk_pid", 0.54, 1.}, {"trk_t0err", 0., 0.85},
               {"trk_cos", tandip_to_cos(0.575), tandip_to_cos(0.85)},
               {"trk_nactive", 19.5, 60.5}, {"trk_t0", 540., 1650.}}},
  };
}

std::vector<Var_t> apply_seed(std::vector<Var_t> vars, const std::vector<SeedCut_t>& cuts) {
  for(const auto& cut : cuts) {
    bool found = false;
    for(auto& var : vars) {
      if(var.name != cut.name) continue;
      if(var.type != kRight) var.min = cut.min;
      if(var.type != kLeft ) var.max = cut.max;
      var.applied = found = true;
    }
    if(!found) cout << "Seed variable " << cut.name << " is not in the variable list, skipping it\n";
  }
  return vars;
}

std::pair<TGraph*,TGraph*> perform_grid_scan(const Events_t& events,
                                             const std::vector<Var_t>& vars,
                                             const double n_bkg_0, const double n_sig_0,
                                             const int ntests = 1e4,
                                             TGraph* scan_points = nullptr,
                                             const double step = 0. // only report points losing more signal than this
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
    const auto [n_sig, n_bkg] = evaluate_yields(vars_test, events);
    const double sig = get_sig(n_sig, n_bkg) / sig_0;
    const double eff = n_sig / n_sig_0;
    if(n_bkg > 0.) {
      const double rej = 1. - n_bkg / n_bkg_0;
      points->AddPoint(eff, rej);
      sigs  ->AddPoint(eff, sig);
      // Compare to the scan points
      if(scan_points && eff > scan_points->GetPointX(scan_points->GetN()-1) && 1. - eff > step) {
        double scan_eff = eff;
        double scan_rej = scan_points->Eval(eff);
        if((rej <  0.5 && rej > scan_rej*1.01) || // 1% better
           (rej <= 0.5 && (1. - rej) < (1. - scan_rej)*0.99)) {
          cout << __func__ << ": Found a selection above the scan! (eff, rej) = ("
               << eff << ", " << rej << ") vs. scan  (" << scan_eff << ", " << scan_rej
               << ") cuts = " << build_cut_string(vars_test) << endl;
        }
      }
    }
  }
  return std::pair<TGraph*, TGraph*>(points,sigs);
}

// optimize_crv: drop the ConvAna CRV veto (add the vetoed set, base_set + 1000) and optimize the CRV veto window
// target: stopping-target design of targets_ (mu- --> e- only). The signal and DIO are the target's unmixed OnSpill
//         samples, normalized to its stops per POT (also applied to the RMC), using set 81 (see optimize_target_selection_fast)
// n_random: randomized greedy starts (each step one of the top 3 candidates), refined like the other starts
int optimize_selection_fast(const bool mumem = false, const int base_set = 40, const bool optimize_crv = false,
                            const TString target = "", const int ngrid_tests = 1e4, const int n_random = 5) {

  // Retrieve the relevant backgrounds
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
  printf("Figure of merit: %s\n", (fom_ == kExpectedUL) ? "1 / expected 90% CL CLs upper limit"
                                 : (fom_ == kAsimovZ)  ? "Asimov significance" : "S / sqrt(B)");

  // Selection applied to all samples that isn't optimized, besides the momentum window
  TString preselection = Form("(%s * trk_charge) > 0 && weight > 0.", (mumem) ? "-1" : "1");
  if(do_target) { // set 81 has no material, cluster, or CRV requirements; the ST intersection cuts are left out
    preselection += " && trk_cluster > 0. && trk_opa == 0 && trk_ntsda == 0";
    if(!optimize_crv) preselection += " && !(trk_crv_raw_dt > 0. && trk_crv_raw_dt < 150.)"; // closest CRV cluster
  }
  if(mumem) printf("Momentum window: %.1f - %.1f MeV/c, flat backgrounds from %.1f - %.1f MeV/c\n",
                   p_sr_min_, p_sr_max_, p_flat_min_, p_flat_max_);
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
    // Momentum window: the flat backgrounds are taken from the wide window, scaled to the signal region width
    const bool flat = mumem && is_flat_bkg(name);
    const double p_lo = (!mumem) ? 90. : (flat) ? p_flat_min_ : p_sr_min_;
    const double p_hi = (!mumem) ? 93. : (flat) ? p_flat_max_ : p_sr_max_;
    const double flat_scale = (flat) ? (p_sr_max_ - p_sr_min_) / (p_flat_max_ - p_flat_min_) : 1.;
    const TString sample_preselection = Form("trk_p > %.2f && trk_p < %.2f && ", p_lo, p_hi) + preselection;
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
      t = t->CopyTree(sample_preselection, "", max_events);
      if(!t) {
        cout << "Dataset " << name << " set " << set << " has no events that pass!\n";
        continue;
      }
      const double scale = events_scale*evaluate_norm(name, f, info);
      samples.push_back(Sample_t(f, t, scale*info.theory_*flat_scale, name));
      samples.back().signal = name == "mumem" || name == "mumep";
      cout << name << " (set " << set << "): N(entries) = " << t->GetEntries() << " scale = " << scale*flat_scale
           << " yield = " << draw_yield(samples.back(), "1") << endl;
      if(flat) { // check the flat approximation at the preselection
        const TString in_sr = Form("trk_p > %.2f && trk_p < %.2f", p_sr_min_, p_sr_max_);
        const double total = draw_yield(samples.back(), "1");
        const double frac = (total > 0.) ? draw_yield(samples.back(), in_sr) / total : 0.;
        printf("  Flat check: %.3f of the %s yield is in the signal region (flat: %.3f, %lld entries in it)\n",
               frac, name.Data(), flat_scale, t->GetEntries(in_sr));
      }
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

  // Load the events into memory, split into the optimization (train) and evaluation (test) halves
  TStopwatch timer;
  Events_t all, train_split, test_split;
  std::vector<double> event_numbers;
  if(!load_events(samples, vars, all, event_numbers)) return 1;
  if(split_) split_events(all, event_numbers, train_split, test_split);
  const Events_t& train = (split_) ? train_split : all;
  const Events_t& test  = (split_) ? test_split  : all;
  printf("Loaded %i events (%i train, %i test) in %.1f s\n", all.size(), train.size(), test.size(), timer.RealTime());
  expected_ul_events(1.); // build the upper limit table

  // Baselines of each half
  const auto [n_sig_0, n_bkg_0] = evaluate_yields(vars, train);
  const auto [n_sig_0_test, n_bkg_0_test] = evaluate_yields(vars, test);
  const double sig_0 = get_sig(n_sig_0, n_bkg_0);
  const double sig_0_test = get_sig(n_sig_0_test, n_bkg_0_test);
  const double lim_0 = get_inv_limit(n_sig_0, n_bkg_0);
  const double lim_0_test = get_inv_limit(n_sig_0_test, n_bkg_0_test);
  printf("Train: N(signal) = %.3g, N(background) = %.3g, significance = %.3g\n", n_sig_0, n_bkg_0, sig_0);
  printf("Test : N(signal) = %.3g, N(background) = %.3g, significance = %.3g\n", n_sig_0_test, n_bkg_0_test, sig_0_test);
  // efficiency, relative significance, and relative 1 / limit of a selection on the test half
  auto evaluate_test = [&](const std::vector<Var_t>& cuts) {
    const auto [s, b] = evaluate_yields(cuts, test);
    return std::array<double,3>{s / n_sig_0_test, get_sig(s, b) / sig_0_test, get_inv_limit(s, b) / lim_0_test};
  };

  // File to output results into
  TString tag = (mumem) ? "mumem" : "mumep";
  if(do_target) tag += "_tgt_" + target;
  if(optimize_crv) tag += "_crv";
  TString f_out = "cut_strings_" + tag + ".txt";
  TString f_val = "cut_validation_" + tag + ".txt";
  TString f_starts = "cut_strings_" + tag + "_starts.txt";
  for(TString f : {f_out, f_val, f_starts}) gSystem->Exec(Form("[ -f %s ] && rm %s; touch %s", f.Data(), f.Data(), f.Data()));
  gSystem->Exec(Form("echo \"# eff_train sig_train inv_limit_train eff_test sig_test inv_limit_test"
                     " (significance and 1 / expected limit relative to the preselection)\" >> %s", f_val.Data()));

  TGraph points, sigs, sigs_test, lims, lims_test;
  points.AddPoint(1., 0.); // first point is 100% efficiencies
  for(TGraph* g : {&sigs, &sigs_test, &lims, &lims_test}) g->AddPoint(1., 1.); // first point is the preselection

  // Greedy path, optimized on the train half
  timer.Start();
  const double eff_step = (use_combine_) ? 0.03 : (mumem) ? 0.01 : 0.005; // aim for 1% efficiency steps
  const double end_eff = (mumem) ? 0.3 : 0.8;
  const auto path = run_greedy(vars, train, n_sig_0, eff_step, end_eff);
  double max_sig(0.), max_lim(0.);
  for(const auto& point : path) {
    const TString cut_string = build_cut_string(point.vars);
    const double sig = get_sig(point.n_sig, point.n_bkg) / sig_0;
    const double lim = get_inv_limit(point.n_sig, point.n_bkg) / lim_0;
    const double eff = point.n_sig / n_sig_0;
    const auto [eff_test, sig_test, lim_test] = evaluate_test(point.vars);
    printf("Eff = %.3f, significance = %.3f, 1/limit = %.3f (test: %.3f, %.3f, %.3f): %s\n",
           eff, sig, lim, eff_test, sig_test, lim_test, cut_string.Data());
    if(point.n_bkg > 0.) {
      points.AddPoint(eff, 1. - point.n_bkg / n_bkg_0);
      sigs  .AddPoint(eff, sig);
      sigs_test.AddPoint(eff_test, sig_test);
      if(sig > max_sig) max_sig = sig;
    }
    lims.AddPoint(eff, lim);
    lims_test.AddPoint(eff_test, lim_test);
    if(lim > max_lim) max_lim = lim;

    // Add the result to the files
    gSystem->Exec(Form("echo %.6f %.6g %.4f \"%s\" >> %s",
                       eff, point.n_bkg/n_bkg_0, sig, cut_string.Data(), f_out.Data()));
    gSystem->Exec(Form("echo %.6f %.4f %.4f %.6f %.4f %.4f >> %s", eff, sig, lim, eff_test, sig_test, lim_test, f_val.Data()));
  }
  printf("Greedy optimization took %.1f s\n", timer.RealTime());
  if(path.empty()) {
    cout << "No selection found\n";
    return 1;
  }

  // Cross-check the in-memory yields of the final selection against the trees
  {
    const auto& final_vars = path.back().vars;
    const auto [n_sig_tr, n_bkg_tr] = evaluate_yields(final_vars, train);
    const auto [n_sig_te, n_bkg_te] = evaluate_yields(final_vars, test);
    const double norm = (split_) ? 0.5 : 1.; // each half estimates the full yields
    const double n_sig_f = (split_) ? norm*(n_sig_tr + n_sig_te) : n_sig_tr;
    const double n_bkg_f = (split_) ? norm*(n_bkg_tr + n_bkg_te) : n_bkg_tr;
    const auto [n_sig_d, n_bkg_d] = draw_yields(samples, build_cut_string(final_vars));
    printf("Cross-check of the final selection: N(signal) = %.6g (trees: %.6g), N(background) = %.6g (trees: %.6g)\n",
           n_sig_f, n_sig_d, n_bkg_f, n_bkg_d);
  }

  // Local search from several starting points
  timer.Start();
  std::vector<std::pair<TString,Selection_t>> starts;
  starts.push_back({"greedy", best_on_path(path)});
  if(do_target) {
    for(const auto& [name, cuts] : seed_selections()) {
      auto seed_vars = apply_seed(vars, cuts);
      const auto [s, b] = evaluate_yields(seed_vars, train);
      starts.push_back({name, {seed_vars, s, b}});
    }
  }
  for(int irun = 0; irun < n_random; ++irun) {
    TRandom3 rnd(1000 + irun);
    const auto random_path = run_greedy(vars, train, n_sig_0, eff_step, end_eff, &rnd, 3);
    if(!random_path.empty()) starts.push_back({Form("random%i", irun), best_on_path(random_path)});
  }
  printf("\n%-10s %8s %8s %8s %8s %8s %8s %8s %8s\n", "Start", "eff", "sig", "refined:", "eff", "sig", "1/limit",
         "test sig", "test 1/l");
  int best_start = -1;
  double best_start_sig = 0.;
  std::vector<Selection_t> refined;
  for(int i = 0; i < (int) starts.size(); ++i) {
    const auto& [name, start] = starts[i];
    refined.push_back(refine(start.vars, train, eff_step));
    const auto& result = refined.back();
    const double sig = get_sig(result.n_sig, result.n_bkg) / sig_0;
    const double lim = get_inv_limit(result.n_sig, result.n_bkg) / lim_0;
    const auto [eff_test, sig_test, lim_test] = evaluate_test(result.vars);
    printf("%-10s %8.3f %8.3f %8s %8.3f %8.3f %8.3f %8.3f %8.3f\n", name.Data(), start.n_sig / n_sig_0,
           get_sig(start.n_sig, start.n_bkg) / sig_0, "", result.n_sig / n_sig_0, sig, lim, sig_test, lim_test);
    gSystem->Exec(Form("echo %s %.6f %.4f %.4f %.6f %.4f %.4f \"%s\" >> %s", name.Data(), result.n_sig / n_sig_0, sig, lim,
                       eff_test, sig_test, lim_test, build_cut_string(result.vars).Data(), f_starts.Data()));
    if(sig > best_start_sig) { best_start_sig = sig; best_start = i; }
  }
  printf("Local search took %.1f s\n", timer.RealTime());
  if(best_start >= 0) {
    const auto [eff_test, sig_test, lim_test] = evaluate_test(refined[best_start].vars);
    printf("Best selection (start %s): eff = %.3f, significance = %.3f (test: %.3f, %.3f, 1/limit %.3f)\n  %s\n",
           starts[best_start].first.Data(), refined[best_start].n_sig / n_sig_0, best_start_sig, eff_test, sig_test,
           lim_test, build_cut_string(refined[best_start].vars).Data());
  }

  // Do a random grid test
  timer.Start();
  auto [grid_points, grid_sigs] = perform_grid_scan(train, vars_grid, n_bkg_0, n_sig_0, ngrid_tests,
                                                    &points, eff_step
                                                    );
  printf("Grid scan took %.1f s\n", timer.RealTime());

  // re-scale the significances and 1 / limits to be 0 - 1, each by its train maximum
  if(max_sig > 0.) {
    cout << "Re-scaling significances by 1 / " << max_sig << endl;
    for(TGraph* g : {&sigs, &sigs_test, grid_sigs})
      for(int point = 0; point < g->GetN(); ++point) g->SetPointY(point, g->GetPointY(point)/max_sig);
  }
  if(max_lim > 0.) {
    cout << "Re-scaling 1 / limits by 1 / " << max_lim << endl;
    for(TGraph* g : {&lims, &lims_test})
      for(int point = 0; point < g->GetN(); ++point) g->SetPointY(point, g->GetPointY(point)/max_lim);
  }

  // Plot the ROC curve
  TCanvas c;
  points.SetLineWidth(2);
  points.SetLineColor(kRed);
  points.Draw("AL");
  sigs  .SetLineWidth(2);
  sigs  .SetLineColor(kGreen);
  sigs  .Draw("L");
  sigs_test.SetLineWidth(2);
  sigs_test.SetLineStyle(kDashed);
  sigs_test.SetLineColor(kGreen+2);
  sigs_test.Draw("L");
  lims.SetLineWidth(2);
  lims.SetLineColor(kBlue);
  lims.Draw("L");
  lims_test.SetLineWidth(2);
  lims_test.SetLineStyle(kDashed);
  lims_test.SetLineColor(kBlue+2);
  lims_test.Draw("L");
  points.GetXaxis()->SetRangeUser(0., 1.1);
  points.GetYaxis()->SetRangeUser(0., 1.1);
  points.SetTitle(";Signal efficiency;Rejection, relative S/#sqrt{B} and 1/limit");
  c.SetGrid();
  TLegend leg(0.12, 0.12, 0.55, 0.40);
  leg.AddEntry(&points   , "Background rejection", "L");
  leg.AddEntry(&sigs     , "S/#sqrt{B} (train)", "L");
  leg.AddEntry(&sigs_test, "S/#sqrt{B} (test)", "L");
  leg.AddEntry(&lims     , "1 / expected 90% CL limit (train)", "L");
  leg.AddEntry(&lims_test, "1 / expected 90% CL limit (test)", "L");
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
int optimize_target_selection_fast(const TString target = "nominal", const bool optimize_crv = false,
                                   const int ngrid_tests = 1e4, const int n_random = 5) {
  return optimize_selection_fast(true, 81, optimize_crv, target, ngrid_tests, n_random);
}
