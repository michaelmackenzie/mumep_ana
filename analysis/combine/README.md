# Combine Fit Workflow (mumep_ana)

This directory contains fit and model-building macros used to produce workspaces and datacards for statistical analysis.

## Organization

- `model_io_utils.C`
  - Shared histogram I/O and normalization-tree scaling helpers.
  - Supports loading shapes from one set or merged control-region set lists.
- `fit_workspace_utils.C`
  - Shared fit orchestration helpers.
  - PDF selection (`default/auto/hist/uniform/poly/polyN/cb/analytic`).
  - Shared right-tail smoothing (`exp/power/convolution/none`).
- `*_fit.C`
  - Component-specific fits that now call shared fit/model helpers.
- `background_model.C`, `signal_model.C`
  - Shared model assembly and workspace-reading helpers.
- `create_envelope.C`
  - Data-driven background envelope: fits function families to the data sidebands and packs the
    surviving functions into a `RooMultiPdf` that Combine profiles over.
- `perform_fits.sh`
  - Executes per-component fits for a selection.
- `full_loop.sh`
  - End-to-end: perform fits, build model, merge cards, run combine.

## Background Envelope (mumep)

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

### Function families

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

### Output

`build_background_envelope` writes `<figdir>/<process>_<selection>_env_functions.png` (and a `_log`
version) showing every accepted function against the data it was fit to, with the blinded region
shaded, a chi^2/dof per function in the legend, and per-bin pulls underneath.

The comparison file `workspaces/comp_<process>_<selection>_<tag>.root` keeps the fit model in
`category/` (signal, envelope `env`, and any explicitly modeled backgrounds), so sums over that
directory stay meaningful. The MC templates of the backgrounds the envelope absorbed are written to
`mc_templates/` at their expected yields, with their raw/normalized/smoothed histograms in
`workflow/`.

## Envelope Likelihood Scans

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

## Control-Region Set Lists

Component fits can merge shape inputs from multiple control-region sets.

- `cosmic_fit`: `control_region_sets` argument
- `dio_fit`: `shape_sets` argument

If not provided, `cosmic_fit` defaults to `{selection + 1000}`.

## Combine Control Regions

Unlike the set lists above, which only feed the component shape fits, a Combine control region is
an extra channel in the merged card. Its processes and nuisances (`rmcN`, `lumi`, ...) share
names with the signal region, so the control region constrains them in the joint fit.

- `./full_loop.sh ... --control-regions "1000 2000 4000"` adds, for each selection `S` and offset
  `O`, the channel `crO_S` built from set `S + O`. The offsets are the region flags of
  `Mu2eEvtAna/src/ConvAna.cc`, and they combine (e.g. 3000 is CRV-tagged and time-inverted):
  - 1000: inverted CRV veto (CRV-tagged tracks), mainly constraining cosmics (`csmN`)
  - 2000: inverted time cut; for set 75 the early time region (t0 <= 540 ns), for sets 20/40
    t0 outside 475-1650 ns, mainly constraining the prompt backgrounds (`rpcN`, `pbrN`)
  - 4000: charge-inverted tracks (e+ for `mumem`, e- for `mumep`), constraining RMC
- `regions.C` defines the regions. The CRV-tagged and time-inverted regions keep the signal
  region's model, including the signal that leaks into them. A charge-inverted region has no
  signal, gets DIO only for e- tracks, and always models both RMC sources; `perform_fits.sh`
  skips the signal and systematics fits there.
- The cosmic shape fit adds the CRV-tagged version of a set (`S + 1000`) to its input unless the
  set already is CRV-tagged.
- 2D (p, t0) fits are not available for time-inverted regions.
- The momentum windows are set separately per process in `defaults.C`: `xmin_ci_em_`/`xmax_ci_em_`
  (mumem control region) and `xmin_ci_ep_`/`xmax_ci_ep_` (mumep control region).

## Data Card Systematics

`tools/write_datacard.C` writes the systematics of the shape and counting cards:

- Rate uncertainties (`lumi`, `csmN`, `dioN`, `rpcN`, `pbrN`, `rmcN`, `intN`): lnN, named without
  the selection, so they are correlated across all channels of a merged card.
- Shape uncertainties: `shape` lines in the shape card. The counting card carries their
  normalization effect as asymmetric lnN.
- MC statistical uncertainties, off by default (`./full_loop.sh ... --mc-stat`, or the
  `include_mc_stat` argument of `build_model`): one lnN per process and channel,
  `<process>_<selection>_<component>_mcstat`, from the bin errors of the process's input histogram.
  Only processes whose PDF is the MC histogram itself get one. A fit or smoothed shape has none;
  the component fit records which kind it saved as `<process>_<selection>_<component>_mc_template`
  in its workspace, and workspaces from before this flag count as fit/smoothed.
- Nuisance groups `rate`, `mcstat`, and `shape`, e.g. for `--freezeNuisanceGroups` or impact
  breakdowns.

## Sensitivity vs. Fit Range

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

## Selection Variations

`make_selection_hists.C` writes dataset histogram files for a tightened selection, for fast
selection sensitivity tests without reprocessing the EventNtuples (EventNtuple inputs only):

```bash
root.exe -q -b 'make_selection_hists.C("40", "trk_qual > 0.9 && trk_pid > 0.95", "selection_hists/qual_pid")'
CONVANA_HIST_PATH=selection_hists/qual_pid ./full_loop.sh --evt-ana -p mumep -s "40" -t qual_pid
```

- The model inputs of each requested set (`obs`, `t0`, and a few track variables) are re-filled from
  the set's selection tree, `Ana/Hist/trs_<set>/tree`, for tracks passing the extra cut, with the
  binning of the original histograms. The cut uses the tree branch names (`trk_p`, `trk_t0`,
  `trk_qual`, `trk_pid`, `trk_d0`, `trk_rmax`, ...) and can only tighten the set's selection. Only
  sets with a selection tree can be re-filled (e.g. 40, 75, 81, and their control-region sets 1040,
  2040, ...), so pass the control-region sets too when the model uses them (`"40 1040 2040"`).
- The output files keep the input names and layout, with the `Norm` tree for the normalization, and
  `obs`/`t0` of all other sets copied unchanged, so shape and control-region sets still resolve.
  `CONVANA_HIST_PATH` points the model building at them in place of the default histogram directory.
- Several variations are made in one pass over the trees with a cut file, one `<tag> : <cut>` per
  line, each written to `<out_dir>/<tag>/`:
  `root.exe -q -b 'make_selection_hists.C("40", "cuts.txt", "selection_hists")'`
- With no cut, the outputs reproduce the original `obs` and `t0` histograms exactly.

`selection_test.sh` runs a full sensitivity test of a selection:

```bash
./selection_test.sh "trk_t0 > 540" t0_540
```

1. `make_selection_hists.C` writes the histograms to `test_hists/<tag>/`.
2. `full_loop.sh` builds a histogram-based model (`<tag>_hists`: `CONVANA_HIST_PDFS=1`,
   `CONVANA_USE_ENV=0`) and an envelope model (`<tag>_env`: `CONVANA_HIST_PDFS=0`,
   `CONVANA_USE_ENV=1`) from them, in parallel, without the systematic fits (`--no-sys`), since
   the test histograms only have the nominal set histograms.
3. Combine generates the Asimov data set of the histogram-based model and fits the envelope model to
   it (`AsymptoticLimits` on that data set). The limit on it is the sensitivity estimate, and does
   not depend on which envelope function is the default. The histogram model's own Asimov limit is
   computed as a reference.

The limits are written to `test_hists/<tag>/summary.txt`, with the logs in `logs/` and the Combine
outputs in `combine/`. Use `-s`/`-c` for other sets and control regions, `--loop-args` and
`--fit-args` to pass options to `full_loop.sh` and Combine, and `--skip-hists`/`--skip-models` to
reuse earlier steps.

### Model configuration from the shell

`defaults.C` reads these environment variables, so the model setup can be chosen without editing it:

- `CONVANA_HIST_PATH`: histogram directory, in place of the default one
- `CONVANA_HIST_PDFS`: `1` for histogram PDFs, `0` for fit functions (`hist_pdfs_`)
- `CONVANA_USE_ENV`: `1` to model the mumep background with the envelope (`use_env_`)

Unset variables keep the values in `defaults.C`.

## Driver Script Configuration

`perform_fits.sh` and `full_loop.sh` support environment overrides:

- Global:
  - `FIT_PDF_TYPE`
  - `FIT_TAIL_MODEL`
  - `FIT_SHAPE_SETS`
  - `FIT_CONTROL_SETS`
- Component-specific:
  - `FIT_PDF_TYPE_<COMP>`
  - `FIT_TAIL_MODEL_<COMP>`
  - `FIT_SHAPE_SETS_<COMP>`
  - `FIT_CONTROL_SETS_<COMP>`
- Set-list controls:
  - any of the shape/control variables above can be comma-separated integer lists

Supported `<COMP>` values:
`SIGNAL DIO COSMIC RPC_EXT RPC_INT PBAR RMC_EXT RMC_INT`

For `mumep` the RMC components are split by neutron knockout, so the RMC overrides use
`RMC_EXT_0N RMC_EXT_1N RMC_INT_0N RMC_INT_1N` instead of `RMC_EXT`/`RMC_INT`.

Example:

```bash
export FIT_PDF_TYPE=default
export FIT_PDF_TYPE_SIGNAL=cb
export FIT_TAIL_MODEL_DIO=convolution
export FIT_SHAPE_SETS_DIO=20,30
export FIT_CONTROL_SETS_COSMIC=1020,1030
./full_loop.sh --process mumem --selections "20 30" --tag test
```
