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
- `perform_fits.sh`
  - Executes per-component fits for a selection.
- `full_loop.sh`
  - End-to-end: perform fits, build model, merge cards, run combine.

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

Example:

```bash
export FIT_PDF_TYPE=default
export FIT_PDF_TYPE_SIGNAL=cb
export FIT_TAIL_MODEL_DIO=convolution
export FIT_SHAPE_SETS_DIO=20,30
export FIT_CONTROL_SETS_COSMIC=1020,1030
./full_loop.sh --process mumem --selections "20 30" --tag test
```
