# PID plots

Plots for the particle identification (PID) section of the mu- -> e- and mu- -> e+ conversion papers.

The PID is the v1 model, a small dense network (4 -> 5 -> 10 -> 5 -> 1, ReLU, sigmoid output). It was trained in
[MLTrain/TrkPID](https://github.com/Mu2e/MLTrain/tree/main/TrkPID) (`TrackPIDTrain.py`, version 1) on flat-spectrum
e- (signal) and mu- (background) tracks from muon stops. The samples were
`nts.mu2e.Flat{e,Mu}MinusOnSpill-reco-ntuple.MDC2025-004.root`, skimmed with `make_inputs.py`. The model has four
inputs:
- E(cluster) / p(track), with p taken at the tracker entrance
- t(track) - t(cluster), unbiased
- the track fit consistency
- the slope of straw-hit time vs. track time (`trkdtdt_slope`, added by `make_inputs.py`)

Offline runs it in `ArtAnalysis/TrkDiag` (`TrackPID_module.cc`, `data/TrackPID_v1.dat`). The Run 1A selection
requires a calorimeter cluster and a score above 0.67 (`Mu2eEvtAna` `ConvAna::Run1ATrackID`).

## Environment

```bash
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh
pyenv rootana 2.5.0
cd mumep_ana/studies/pid
```

## Training-sample plots: `make_training_plots.py`

This script reads the skimmed training ntuples (default `/exp/mu2e/data/users/mmackenz/trkpid/data/`). It
reproduces the training track selection and train/test split, which match the training csv files row by row.
It then evaluates the deployed weights (`ArtAnalysis/TrkDiag/data/TrackPID_v1.dat`, which agree with the trained
Keras model to < 1e-6). The training ntuples predate v1, so their `trkpid` branch holds the v0 score and is not used.

```bash
python make_training_plots.py          # ~1 min; the selection is cached in outputs/
```

Outputs:
- `figures/training/`:
  - `pid_features*`: input feature distributions
  - `pid_score`: score for the test and training samples
  - `pid_roc`: mu- vs. e- efficiency
  - `pid_eff_vs_cut`: efficiencies vs. the score cut
  - `pid_eff_vs_p`: efficiencies at the working point vs. momentum
- `outputs/pid_summary.json`: efficiencies at the working point, overall and in the 100-110 and 87-97 MeV/c
  windows, and the test-sample AUC

## ConvAna plots: `make_convana_plots.py`

These plots use the Offline v1 score (the `trkpid` branch of the MDC2025au-001 ntuples), as processed by ConvAna. They
add the effect of pileup, e+ tracks (not in the training), and the conversion signal and cosmic samples. They use
the track sets filled before the PID requirement: set 4 (e-, p > 95 MeV/c) and set 5 (e+, p > 75 MeV/c).

The flat samples are in `Mu2eEvtAna/scripts/datasets.C` but are not processed by default:
- `fele0b0s5r0101`, `fele0b1s5r0101`: flat e-, without and with pileup
- `fpos0b0s5r0101`, `fpos0b1s5r0101`: flat e+, without and with pileup
- `fmum0b0s5r0101`: flat mu-, without pileup

Process them with, e.g.:
```bash
# in a muse setup of the work area, from the directory that should receive the histograms
root -q -b '${MUSE_WORK_DIR}/Mu2eEvtAna/scripts/make_histograms.C(1, "fele0b0s5r0101", 1, "cnv_ana", 10, -1)'
```
Then make the plots in `figures/convana/`:
```bash
python make_convana_plots.py [--hist-dir <dir>]   # default: $CONVANA_HIST_PATH or the run1a histogram area
```
