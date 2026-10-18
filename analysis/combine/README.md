# Combine Fit Workflow (mumep_ana)

The Combine models are built in Python with mu2ecombine (`Mu2eEvtAna/combine`, see its README): the component shapes
are fit, the model of all selection sets and control regions is assembled, and Combine cards are written. Combine
does the statistics.

## Organization

- `samples.py`: the Mu2eEvtAna (EventNtuple) datasets, their rates (as in `../datasets.C` and `../physics.C`), the
  exposures (Run 1A at 1BB or 2BB, Run 2, MDS3c/d), and the rate uncertainties.
- `build_model.py`: builds and writes one model, the shape model and the cut-and-count model.
- `full_loop.py`: builds the model and optionally tests it with Combine.
- `selection_test.py`: sensitivity test of a tightened selection (with `make_selection_hists.C`).
- `sensitivity_vs_fit_range.py`: expected sensitivity vs. the momentum fit range.
- `envelope_scans.py`: likelihood scans of r for each envelope function and with the index profiled.
- `sensitivity_vs_livetime.py`: expected sensitivity vs. a livetime scale factor.
- `selection_comparisons.py`: Asimov expected limits of several models (shape and cut-and-count, with and without
  the constrained nuisances, toys); `--single-test` for only the shape limit with all uncertainties.

## Setup

Set up Combine as in the top-level README (`pyenv rootana`, `env_standalone_mu2e.sh`). mu2ecombine is taken from
`PYTHONPATH`, `MU2ECOMBINE_DIR`, or a `Mu2eEvtAna` clone next to this repository:

```bash
export PYTHONPATH=/path/to/Mu2eEvtAna/combine:$PYTHONPATH
```

## Building a model

```bash
./full_loop.py -p mumep -s "40" -t v1                        # template model of set 40
./full_loop.py -p mumep -s "40" -t v1_env --envelope --plots # data-driven envelope background
./full_loop.py -p mumem -s "75" -c "1000 2000 4000"           # with control regions
./full_loop.py -p mumep -s "40" --2bb --test-combine          # 2BB exposure, then limits and a fit
```

To develop or debug one component, build and fit only it: `--only <process>` (e.g. `signal`, `cosmic`,
`rmc_ext_0n`, `env`) prints its fit result, draws it to `<output>/debug/`, and exits without writing a model.
`--channel` picks one channel; the `FIT_*` overrides apply as usual, e.g.

```bash
FIT_PDF_TYPE_COSMIC=poly2 ./build_model.py -p mumep -s 40 --only cosmic
```

`full_loop.py` takes the `build_model.py` options (`./build_model.py -h`) plus `--test-combine` (text2workspace, the
expected 90% CL limit, FitDiagnostics, and pre/post-fit figures), `--fit-args`, `--skip-model`, and `--dry-run`
(print the Combine commands).

The output, `models/<process>_<sets>[_cr<offset>...][_<tag>]/` (or `-o`):

| File | Content |
|---|---|
| `combined.txt`, `<channel>.txt` | Datacards of all channels / one channel; paths are relative to the card |
| `shapes/<channel>.root` | The workspace `w`: observable `obs_<set>`, templates as RooDataHist, PDFs, data |
| `model.json` | Yields, MC uncertainties, nuisances, fit results, and envelope choices |
| `commands.sh` | text2workspace and common Combine commands |
| `figures/` | Expected stacks and the MC shape fits (`--plots`) |
| `cc/` | The cut-and-count model, in the same layout |

Every channel is a RooWorkspace on its observable `obs_<set>` (mu2ecombine `Model(..., workspace=True)`), also when
all its processes are templates. Any model of the same channels can therefore generate data that another one fits,
e.g. `combine -M GenerateOnly --saveToys` with a template model and `--toysFile` with an envelope or function model.

The card headers record the signal reference branching fraction (the unit of `r`), N(POT), livetime, N(muons), and
the signal efficiency of each channel.

### Channels and control regions

Each selection set `S` is the channel `set_S`. `-c "1000 2000 4000"` adds, for each set and offset `O`, the channel
`crO_S` from set `S + O`. The offsets are the region flags of `Mu2eEvtAna/src/ConvAna.cc`, and they combine (e.g. 3000
is CRV-tagged and time-inverted):

- 1000: inverted CRV veto (CRV-tagged tracks), mainly constraining cosmics (`csmN`)
- 2000: inverted time cut; for set 75 the early time region (t0 <= 540 ns), for sets 20/40 t0 outside 475-1650 ns,
  mainly constraining the prompt backgrounds (`rpcN`, `pbrN`)
- 4000: charge-inverted tracks (e+ for `mumem`, e- for `mumep`), constraining RMC

Process and nuisance names are shared by all channels, so the control regions constrain them in the joint fit. The
CRV-tagged and time-inverted regions keep the signal region's model, including the signal that leaks into them. A
charge-inverted region has no signal, gets DIO only for e- tracks, and models both RMC sources. The momentum ranges
are 100-110 (mumem), 87-97 (mumep), and 95-105 MeV/c (mumem charge-inverted), in 0.2 MeV/c bins.

### Process shapes

| Process | Shape |
|---|---|
| signal | The MC histogram (`--functions`: a double-sided Crystal Ball fit) |
| dio | The MC histogram; above 102 MeV/c the generated spectrum (set 10) convolved with the momentum response |
| cosmic | Chebychev (2 coefficients) fit to the set plus its CRV-tagged version (set + 1000); rate from the set |
| pbar, rpc_ext / rpc_int | Chebychev fits, 1 / 2 coefficients |
| rmc (mumep) | 0n: Gaussian + linear polynomial; 1n: exp(quadratic), i.e. a broad Gaussian |
| rmc (mumem) | Turn-on x power law x exponential: (p - x0)^a (120 - p)^b exp(-l p), x0 = 96 (94 charge-inverted) |

- A non-RMC background filling less than 75% of the bins of its histogram is modeled as uniform.
- The fitted shapes are binned templates by default (`Fit(..., binned=True)`). With `--functions`
  (`CONVANA_HIST_PDFS=0`), the fitted PDFs enter a RooWorkspace instead.
- Processes with no simulated events in the range are left out.
- mumep splits each RMC source into its 0 and 1 neutron knockout samples; mumem uses the 0n sample.
- Overrides from the environment: `FIT_PDF_TYPE[_<COMP>]` (`auto`, `analytic` (no uniform fallback), `hist`,
  `uniform`, `polyN`, `cb`), `FIT_TAIL_MODEL[_DIO]` (`convolution` or `none`), `FIT_SHAPE_SETS[_<COMP>]` (sets whose
  histograms give the shape; the rate stays the set's), and `FIT_CONTROL_SETS[_COSMIC]`. `<COMP>` is `SIGNAL`, `DIO`,
  `COSMIC`, `RPC_EXT`, `RPC_INT`, `PBAR`, `RMC_EXT`, `RMC_INT`, or `RMC_EXT_0N` ... `RMC_INT_1N` for mumep.

### Background envelope

`--envelope` (`CONVANA_USE_ENV=1`) replaces the MC backgrounds of each signal-region channel by a data-driven
envelope: function families fit to the data outside the blind window (90-93 MeV/c for mumep, 102-105.5 for mumem) and
profiled discretely by Combine (`RooMultiPdf`, index `<channel>_env_cat`). The families and the chi^2/F-test
selection are mu2ecombine's defaults, which follow the configuration of `legacy/create_envelope.C` described below. Control
regions keep the MC backgrounds. Building the envelope takes a few minutes.

### Data

`--data` selects `data_obs`:

- `asimov` (default): the expected backgrounds.
- `toy[:seed]`: the expected backgrounds (MC, without an envelope), Poisson-fluctuated per bin.
- `mds3c`, `mds3d`: the MDS ensembles, with their exposure and a 1e-13 signal reference.
- `<file>:<histogram>`: any histogram; `{set}` in the name is replaced by each channel's set.
- `--data-model <model dir>`: the data of another model, e.g. an envelope model fit to the Asimov data of a template
  model.

### Normalization and uncertainties

- `--exposure` (`run1a`, `run2`, `mds3c`, `mds3d`) and `--2bb` (`CONVANA_USE_2BB=1`: the inputs re-weighted to 2BB,
  `ConvAna.*.m2.root`, at the 2BB duty cycle and N(POT) per microbunch). The signal is normalized to `--reference`
  (default 1e-15).
- Rate uncertainties (`lumi`, `csmN`, `dioN`, `rpcN`, `pbrN`, `rmcN`, `intN`): lnN, correlated across all channels,
  in the nuisance group `rate`.
- `--mc-stat`: the MC statistical uncertainty of each process yield, `mcstat_<channel>_<process>` (group `mcstat`).
- Yield scale factors, frozen at 1: `yieldScale` for all processes and `scale_<process>`. Set them when running
  Combine, e.g. `--setParameters scale_dio=2,scale_cosmic=0.5` or `yieldScale=2`.

### Cut-and-count model

`cc/` holds a one-bin model per channel: the momentum window maximizing S/sqrt(S + B) in a signal region, the full
range in a control region, with the shape model's expected yields (an envelope enters at its expected yield, fixed).
The window is recorded in the card header. `--no-cc` skips it.

## Selection tests

`selection_test.py` runs a full sensitivity test of a selection:

```bash
./selection_test.py "trk_t0 > 540" t0_540 -s "40"
```

1. `make_selection_hists.C` writes the histograms to `test_hists/<tag>/` (see below).
2. `build_model.py` builds a template model (`<tag>_hists`) and an envelope model (`<tag>_env`, `--envelope
   --functions`) whose data is the template model's Asimov data set (`--data-model`).
3. The limit of the envelope model on that data set (its observed limit) is the sensitivity estimate, and does not
   depend on which envelope function is the default. The template model's own Asimov limit is the reference.

The limits are written to `test_hists/<tag>/summary.txt`, with the logs in `logs/` and the Combine outputs in
`combine/`. Use `-p`/`-s`/`-c` for other processes, sets, and control regions, `--model-args` and `--fit-args` to pass
options to `build_model.py` and Combine, and `--skip-hists`/`--skip-models` to reuse earlier steps.

## Sensitivity vs. fit range

`sensitivity_vs_fit_range.py` scans the Asimov expected 90% CL limit over momentum fit ranges:

```bash
./sensitivity_vs_fit_range.py --full-model models/mumep_40_r0106_env --low-edges "87 88 89" --high-edges "96 97"
./sensitivity_vs_fit_range.py -p mumem -s 75 --low-edges "100 101 102 103" --high-edges "105 110" --jobs 4
```

- The full-range model is built once from the `build_model.py` options (in `<out dir>/full/`), or taken from
  `--full-model`, and narrowed to each range with mu2ecombine's `narrow_model`, without refitting: the MC templates
  keep their full-range shapes bin by bin, fitted PDFs and the envelope are renormalized over the range, and the yields
  are scaled to the fraction inside it. Only the range changes between points.
- The envelope keeps all of its full-range functions in every range (their parameters still float in Combine), so a
  range cannot look better because it was refit with fewer functions. Its starting yield is the fraction of its
  default function inside the range, and it still floats freely.
- Only the signal-region channels (`set_<S>`) are narrowed; control regions keep their range. Range edges are moved
  to the nearest bin edges. Rate and MC statistical uncertainties stay as evaluated over the full range.
- `--refit` instead rebuilds each range with `build_model.py --range`: the MC processes keep their full-range shapes,
  but an envelope is refit to the data in each range (a few minutes per range).
- Defaults: lower edges from the range minimum to its midpoint in `--step` (0.5 MeV/c) steps, the upper edge at the
  range maximum. Outputs go to `range_scans/<model name>/`: a model directory per range, logs,
  `sensitivity_vs_fit_range.csv`, and the figures:
  - `sensitivity_vs_fit_range.png`: median expected limit vs. the lower edge, one curve per upper edge (up to 6), with
    the +-1/2 sigma bands for a single curve
  - `sensitivity_vs_fit_range_2d.png`: median expected limit vs. (lower edge, upper edge), with the best range
    outlined, when at least two lower and two upper edges were scanned
- The median expected limit is the observed limit on the Asimov data set (`-t -1`); the +-1/2 sigma bands are the
  expected quantiles. A range whose Asimov limit differs from the expected median by more than 5% is flagged
  (`asimov mismatch`) and left out of the figures: a profiled fit missed the global minimum.
- The fits use minimizer strategy 1, `--noFitAsimov`, and 1% limit accuracy (see the comments in the script).
  `--fit-args` adds Combine arguments, `--jobs` runs the limits in parallel (one core each, at most 20),
  `--no-combine` makes the models and prints the Combine commands, and `--plot-only` remakes the figures from the
  table of an earlier scan.

## Envelope likelihood scans

`envelope_scans.py` scans r for each envelope function with the index frozen, then with the index profiled, and
overlays the curves as -2 Delta ln L:

```bash
./envelope_scans.py models/mumep_40_env/combined.txt -r 200 [--fit-args "..."]
```

- The index is the card's discrete parameter (`<channel>_env_cat`; `--index` for a card with several), and the
  functions and their titles are read from the channel workspace.
- `--fit-args` adds Combine arguments to every fit; `--plot-only` remakes the figure from existing scans; `--dry-run`
  prints the Combine commands. The Combine outputs and the figure `envelope_<model name>.png` are written to the
  working directory.
- Combine stores plain NLL values, so the curves are 2(nll0 + nll + deltaNLL - minimum), with the 68% and 95% CL levels
  at 1 and 4. The profiled scan saves the chosen index at each point, and the script checks that it follows the lowest
  single-function curve.
- The fits do not use `MINIMIZER_multiMin_maskConstraints`: masking the constraints in the loop over the functions
  leaves an NLL offset that accumulates over the points of the profiled scan.

## Sensitivity vs. livetime

`sensitivity_vs_livetime.py` scales all yields with the frozen `yieldScale` factor, generates a background-only
Asimov data set with the generation card at each scale, and fits it with the fit card:

```bash
python3 sensitivity_vs_livetime.py models/mumep_40_env/combined.txt --gen-card models/mumep_40_hists/combined.txt
```

## Selection Variations

`make_selection_hists.C` writes dataset histogram files for a tightened selection, for fast
selection sensitivity tests without reprocessing the EventNtuples (EventNtuple inputs only):

```bash
root.exe -q -b 'make_selection_hists.C("40", "trk_qual > 0.9 && trk_pid > 0.95", "selection_hists/qual_pid")'
./full_loop.py -p mumep -s "40" -t qual_pid --hist-path selection_hists/qual_pid
```

- The model inputs of each requested set (`obs`, `t0`, and a few track variables) are re-filled from
  the set's selection tree, `Ana/Hist/trs_<set>/tree`, for tracks passing the extra cut, with the
  binning of the original histograms. The cut uses the tree branch names (`trk_p`, `trk_t0`,
  `trk_qual`, `trk_pid`, `trk_d0`, `trk_rmax`, ...) and can only tighten the set's selection. Only
  sets with a selection tree can be re-filled (e.g. 40, 75, 81, and their control-region sets 1040,
  2040, ...), so pass the control-region sets too when the model uses them (`"40 1040 2040"`).
- The output files keep the input names and layout, with the `Norm` tree for the normalization, and
  `obs`/`t0` of all other sets copied unchanged, so shape and control-region sets still resolve.
  `--hist-path` (or `CONVANA_HIST_PATH`) points the model building at them in place of the default histogram directory.
- Several variations are made in one pass over the trees with a cut file, one `<tag> : <cut>` per
  line, each written to `<out_dir>/<tag>/`:
  `root.exe -q -b 'make_selection_hists.C("40", "cuts.txt", "selection_hists")'`
- With no cut, the outputs reproduce the original `obs` and `t0` histograms exactly.

## Macro-based model (legacy/)

The previous workflow, ROOT macros (`*_fit.C`, `build_model.C`, `create_envelope.C`, `narrow_model.C`, ...) driven by
`full_loop.sh`, `perform_fits.sh`, and `selection_test.sh`, is kept in `legacy/` for now, together with the versions
of the scan tools that read its outputs. Run them from `legacy/`: they write their cards to `datacards/` and
workspaces to `workspaces/` in the working directory. The sections below describe them.

### Background Envelope (mumep)

With `use_env_ = true` in `defaults.C`, `build_model.C` replaces the mumep background model with a
data-driven envelope:

- The processes named in `explicit_processes` (in `build_model.C`) are summed into a `RooAddPdf`
  with every parameter frozen, including their expected rates. That sum is held fixed while the
  envelope functions are fit, and each of those processes is still written to the workspace and
  listed in the data card with its own rate and uncertainties.
- Every other background is absorbed into the envelope, which is fit to the data sidebands with
  `p_blind_min` to `p_blind_max` excluded.
- The envelope enters the workspace as three objects: the `RooMultiPdf`
  (`<process>_<selection>_env_pdf`), the discrete index it is profiled over
  (`<process>_<selection>_env_cat`), and its freely floating yield
  (`<process>_<selection>_env_pdf_norm`). Its shape parameters are left floating.
- In the data card the envelope appears as the process `env` with rate `1` (Combine multiplies it
  by the `_norm` variable), carries no rate uncertainties, and the discrete index is declared with
  a `<index> discrete` line.

#### Function families

Each family is scanned over a range of orders by `add_family`, and the orders that pass are added to
the envelope. Which families are considered is set by the flags at the top of `create_envelope.C`:

| Flag | Family | Orders scanned | Order used when `force_fit_order_` |
| --- | --- | --- | --- |
| `use_poly_family_` | Chebychev | 1-4 | 3 |
| `use_bernstein_family_` | Bernstein (`use_generic/fast_bernstein_` pick the implementation) | 1-4 | 4 |
| `use_exp_family_` | sum of exponentials | 1-3 | 2 |
| `use_power_family_` | sum of power laws | 1-3 | 2 |
| `use_laurent_family_` | Laurent | 1-6 | 1 |
| `use_inv_poly_family_` | inverse polynomial | 1 | 1 |
| `use_gaussian_family_` | broad Gaussian (`RooGaussian`), skipped when exp(polynomial) order 2 is in the envelope | 1 function | - |
| `use_gaus_poly_family_` | Gaussian + polynomial | 0-3 | 1 |
| `use_gaus_expo_family_` | Gaussian + exponential | 1-2 | 1 |
| `use_gaus_power_family_` | Gaussian + power law | 1-2 | 1 |
| `use_endpoint_family_` | endpoint power law (E0 - p)^a x (1 + polynomial) | 0-2 | 0 |
| `use_exp_poly_family_` | exp(polynomial) | 1-4 | 3 |

How well each family describes the expected mumep background (the sum of the fitted component
shapes, tag `evt_r0105_env_v02`, set 40, 87-96 MeV/c; ~95% radiative muon capture positrons), fit
to the sidebands as the envelope does, with the miss in the blinded 90-93 MeV/c region in units of
its statistical precision (sqrt(2015) events):

| Function | Shape parameters | Blinded-region bias |
| --- | --- | --- |
| endpoint power law, order 0 (E0 = 100.5 MeV/c, a = 3.1) | 2 | -0.05 sigma |
| Bernstein, order 4 | 4 | -0.08 sigma |
| exp(polynomial), order 3 | 3 | +0.18 sigma |
| Gaussian + polynomial, order 0 | 3 | -0.35 sigma |
| Gaussian + power law, order 1 | 4 | -0.69 sigma |
| broad Gaussian (= exp(polynomial), order 2) | 2 | -0.91 sigma |
| sums of exponentials or power laws, Chebychev, inverse polynomial, Laurent | | poor sideband fits |

The endpoint power law follows the physics of this background, a spectrum falling to the radiative
capture endpoint just above the fit range, and needs only two parameters; exp(polynomial) order 3
adds a genuinely different shape with the opposite-sign bias. Gaussian + exponential/power law
order 1 and exp(polynomial) order 2 reduce to the broad Gaussian here, so the broad Gaussian is
only added when exp(polynomial) order 2 was not accepted (or that family is off).

`force_fit_order_` keeps only the fixed order listed above; otherwise the orders are accepted on a
chi^2 test (`chisq_p_min_`) and `enforce_ftest_` stops the scan once a higher order stops improving
the fit (`ftest_p_max_`). The two are mutually exclusive by construction, so the F-test is skipped
when the order is forced.

With `verbose_ > 0`, `add_family` prints the acceptance policy for each family and one line per
order tested -- parameter count, bins used, chi^2/dof, p(chi^2), the F-test p-value where it
applies, and whether the order was accepted or rejected and why -- followed by a line when a scan
stops early. This is what to read to see why a family stopped at a given order or why a function was
not accepted:

```
add_family: ===== Exponential family, orders 1 - 3: keeping orders with p(chi^2) > 0.001, stopping once p(F) > 0.05 =====
add_family:   Exponential order 1:  1 par,  35 bins, chi^2/dof =    39.18 /  33 =  1.187, p = 2.12e-01, ACCEPTED (p(chi^2) = 0.212 > 0.001)
add_family:   Exponential order 2:  3 par,  35 bins, chi^2/dof =    39.19 /  31 =  1.264, p = 1.48e-01, F-test vs order 1: p = 1.00e+00, rejected (no significant improvement over order 1)
add_family:   Exponential: stopping the scan at order 2 (higher orders not tested)
add_family:   Exponential: 1 / 2 tested order(s) added to the envelope
```

The Gaussian + polynomial/exponential/power law families are `RooGenericPdf`s that need a numeric
normalization on every likelihood call, which makes them far slower to fit than the analytic
families; the broad Gaussian is a `RooGaussian` with an analytic normalization.

#### Output

`build_background_envelope` writes `<figdir>/<process>_<selection>_env_functions.png` (and a `_log`
version) showing every accepted function against the data it was fit to, with the blinded region
shaded, a chi^2/dof per function in the legend, and per-bin pulls underneath.

The comparison file `workspaces/comp_<process>_<selection>_<tag>.root` keeps the fit model in
`category/` (signal, envelope `env`, and any explicitly modeled backgrounds), so sums over that
directory stay meaningful. The MC templates of the backgrounds the envelope absorbed are written to
`mc_templates/` at their expected yields, with their raw/normalized/smoothed histograms in
`workflow/`.

### Envelope Likelihood Scans

`envelope_scans.py` (or the `envelope_scans.sh` wrapper) scans r for each envelope function with
the index frozen, then with the index profiled, and overlays the curves as -2 Delta ln L:

```bash
python3 envelope_scans.py datacards/combine_mumep_40_evt_r0104_env.txt -r 200 [--fit-args "..."]
```

- Defaults from the card: the index `<process>_<selection>_env_cat` (checked against the card's
  `discrete` line), the number of functions and their names (from the workspace), the Combine output
  tag `envelope_<card name>`, and the figure `envelope_<card name>.png`.
- `--fit-args` adds Combine arguments to every fit; `--plot-only` remakes the figure from existing
  scans; `--dry-run` prints the Combine commands.
- Combine stores plain NLL values, so the curves are 2(nll0 + nll + deltaNLL - minimum), with the
  68% and 95% CL levels at 1 and 4. The profiled scan saves the chosen index at each point, and the
  script checks that it follows the lowest single-function curve.
- The fits do not use `MINIMIZER_multiMin_maskConstraints`: masking the constraints in the loop over
  the functions leaves an NLL offset that accumulates over the points of the profiled scan.

### Control-Region Set Lists

Component fits can merge shape inputs from multiple control-region sets.

- `cosmic_fit`: `control_region_sets` argument
- `dio_fit`: `shape_sets` argument

If not provided, `cosmic_fit` defaults to `{selection + 1000}`.

### Sensitivity vs. Fit Range

`sensitivity_vs_fit_range.py` scans the Asimov expected 90% CL limit of one set over momentum fit
ranges, starting from the per-set shape card `build_model` writes for the full range:

```bash
python3 sensitivity_vs_fit_range.py datacards/combine_mumem_75_<tag>.txt \
    --low-edges "100 101 102 103" --high-edges "105 110" --jobs 4
```

- Each range gets a model from `narrow_model.C` instead of a refit: the PDFs keep their full-range
  shapes, each card rate is scaled by the fraction of its PDF inside the range, and Combine
  renormalizes the PDFs over the narrowed observable. Inside the range this reproduces the
  full-range model bin by bin, so only the range changes between points, and a point costs one
  short ROOT job plus one Combine call instead of the component fits and `build_model`.
- Range edges are moved to the nearest bin edges. The lnN priors, shape-systematic normalization
  effects, and MC statistical uncertainties stay as evaluated over the full range. Processes with
  no yield in a range are dropped by Combine, and their nuisances are removed from the groups.
- The floating background envelope (mumep) is supported: its card rate is scaled by the fraction
  of the current envelope function inside the range, so the model starts at the envelope's yield
  in the range, while its `<pdf>_norm` still floats and Combine still profiles the function index
  (each function renormalized over the range). The functions keep their full-range sideband fits.
- Only single-channel 1D shape cards are supported (no control regions or 2D fits).
- Defaults: lower edges from the range minimum to its midpoint in `--step` (0.5 MeV/c) steps, the
  upper edge at the range maximum. Outputs go to `range_scans/<card name>/`: the narrowed cards and
  workspaces, logs, `sensitivity_vs_fit_range.csv`, and the figures:
  - `sensitivity_vs_fit_range.png`: median expected limit vs. the lower edge, one curve per upper
    edge (up to 6), with the +-1/2 sigma bands for a single curve
  - `sensitivity_vs_fit_range_2d.png`: median expected limit vs. (lower edge, upper edge), with the
    best range outlined, when at least two lower and two upper edges were scanned
- The median expected limit is the observed limit on the Asimov data set (`-t -1`), which Combine
  computes directly instead of on the coarser grid of its `Expected 50.0%` quantile; the +-1/2 sigma
  bands are the expected quantiles. The table records both (`median` and `50.0`). A range whose
  Asimov limit differs from the expected median by more than 5% is flagged (`asimov mismatch`) and
  left out of the figures: it means a profiled fit missed the global minimum.
- The fits use minimizer strategy 1. With strategy 0, the profiled fits of the limit search can stay
  in a worse minimum of the envelope shape parameters after the large-r points, biasing the Asimov
  limit (down to ~0). They also use `--noFitAsimov`, so the expected bands come from the same
  pre-fit Asimov data set as the median: the regenerated Asimov of a background-only refit can follow
  a different envelope function than the pre-fit one in narrow ranges. `--fit-args` adds Combine
  arguments to every fit.
- `--no-combine` builds the models and prints the Combine commands only. `--plot-only` remakes the
  figures from the `sensitivity_vs_fit_range.csv` of an earlier scan, without ROOT or Combine; for a
  table from before the `median` column, it reads the Asimov limits from the scan's Combine logs.
