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
| `use_poly_family_` | Bernstein | 1-4 | 3 |
| `use_poly_family_` | Chebychev | 3 | 3 |
| `use_exp_family_` | exponential | 1-3 | 2 |
| `use_power_family_` | power law | 1-3 | 2 |
| `use_laurent_family_` | Laurent | 1-6 | 1 |
| `use_inv_poly_family_` | inverse polynomial | 1 | 1 |
| `use_gaus_poly_family_` | Gaussian + polynomial | 0-3 | 1 |
| `use_gaus_expo_family_` | Gaussian + exponential | 1-2 | 1 |
| `use_gaus_power_family_` | Gaussian + power law | 1-2 | 1 |

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

The Gaussian families are `RooGenericPdf`s that need a numeric normalization on every likelihood
call, which makes them far slower to fit than the analytic families.

### Output

`build_background_envelope` writes `<figdir>/<process>_<selection>_env_functions.png` (and a `_log`
version) showing every accepted function against the data it was fit to, with the blinded region
shaded, a chi^2/dof per function in the legend, and per-bin pulls underneath.

## Control-Region Set Lists

Component fits can merge shape inputs from multiple control-region sets.

- `cosmic_fit`: `control_region_sets` argument
- `dio_fit`: `shape_sets` argument

If not provided, `cosmic_fit` defaults to `{selection + 1000}`.

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
