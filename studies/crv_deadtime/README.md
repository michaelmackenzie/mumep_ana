# CRV veto timing window

A CRV veto rejects cosmic-ray tracks with a CRV cluster in time with the track. It also rejects signal tracks
that have an accidental in-time CRV cluster, which is dead time. This study compares the track - CRV time
differences in the ConvAna track trees (raw, propagated via the calorimeter, propagated via the stopping target, and
the closer of the two propagated times) and scans the veto window on each to maximize S / sqrt(B). It uses the
mu- -> e+ selection (ConvAna set 40, 90 < p < 93 MeV/c).

## Summary and takeaways

Last update: 2026-10-27, adding the upstream track window (below). CRV window results: 2026-10-26, adding the
separate propagation hypotheses (MDC2025au_best_v1_1 ntuples, `r0101`
ConvAna histograms, mode 1, re-made in `histograms_crv_dt_2026-10-26`)

| Veto | Signal loss | Cosmic eff. | All bkg. eff. | S/sqrt(B) | S/sqrt(B_cosmic) |
|---|---|---|---|---|---|
| ConvAna: propagated (-50, 60) OR raw (-25, 0) | 0.156% | 0.414 | 0.9920 | 1.0024 | 1.553 |
| **Proposed: raw (0, 172.5) ns** | 0.152% | 0.377 | 0.9919 | 1.0026 | 1.627 |
| Raw optimum (10, 162.5) ns, no margin | 0.135% | 0.377 | 0.9919 | 1.0027 | 1.628 |
| Propagated (-62.5, 85) ns | 0.173% | 0.377 | 0.9916 | 1.0025 | 1.627 |
| Via calorimeter (-95, 85) ns | 0.158% | 0.377 | 0.9919 | 1.0025 | 1.627 |
| Via stopping target (-62.5, 112.5) ns | 0.156% | 0.377 | 0.9916 | 1.0026 | 1.627 |
| Split: via calo (-95, 37.5) OR via target (65, 115) ns | 0.143% | 0.377 | 0.9919 | 1.0026 | 1.627 |

All the proposals except the ConvAna veto remove every CRV-tagged cosmic ray in the MC. Each includes the 10 ns
margin on each edge; the optima without the margin are 0.135% (raw), 0.142% (via calo), 0.137% (via target),
0.155% (closer of the two), and 0.120% (split).

S / sqrt(B) values are relative to no CRV veto.

- **Proposed window: veto if 0 < t(track) - t(CRV, raw) < 172.5 ns.** It removes every CRV-tagged cosmic ray in
  the MC, 101 of 162 tracks. The ConvAna veto misses 6 of them, which are just outside its propagated-time window
  (4 at -52 to -50 ns, 2 at +72 ns). The proposed window has about the same signal loss as the ConvAna veto
  (0.15%).
- **The raw time difference is the better veto variable.** The propagated difference
  (`trk_crv_dt`) takes whichever of two path hypotheses gives the smaller |dt|: through the calorimeter or through
  the stopping target. That gives each accidental CRV cluster two chances to land near zero, so the signal accidentals
  have a peak at 0 (`crv_dt_trk_crv_dt`). The raw difference has a flat accidental distribution. Rejecting all
  tagged cosmic rays then costs 0.135% of signal on the raw time, against 0.155% on the propagated time.
- **The separate propagation hypotheses don't need separate handling.** On its own, each hypothesis has a flat
  accidental distribution like the raw time, and needs a window as wide (~150 ns) to contain all the CRV-tagged
  cosmic rays: -83 to 73 ns via the calorimeter, -52 to 101 ns via the target, against 11 to 160 ns raw. For the
  cosmic rays the two hypotheses are tightly correlated (`crv_dt_2d_trk_crv_cal_dt_trk_crv_st_dt`): each track lies
  on one of two bands of fixed offset, so the ~150 ns spread is a real timing spread that neither propagation
  removes, not an ambiguity in the path. Taking the closer of the two (as ConvAna does) is the worst option, at
  0.155% signal loss.
- **Two separate windows save little.** The exact search for a calorimeter window OR a target window that vetoes
  every tagged cosmic ray at the least signal loss gives 0.120% (0.143% with margins), against 0.135% (0.152%) for
  the single raw window. The target window is there only for the 6 cosmic-ray tracks in the late tail (via target
  79-101 ns), so the gain rests on 6 MC events and adds a window. The single raw window stays the proposal.
- **The all-background S/sqrt(B) can't choose a window.** In 90-93 MeV/c, RMC is ~98% of the background, so every
  veto choice moves S/sqrt(B) by < 0.3%, and the differences between windows (~0.03%) are within the MC
  statistics. So the proposal maximizes S/sqrt(B) with only the cosmic-ray background (`--fom cosmic`). The momentum
  fit separates RMC by its spectrum shape, and cosmic rays are the background the veto targets. The
  all-background S/sqrt(B) is shown in every plot and table for comparison (`--fom all` uses it for the proposal).
- **Signal dead time is ~8.5e-6 per ns of raw window** (Mix1BB pileup), i.e. ~0.009% per 10 ns. The window
  margin is cheap: the 10 ns added to each side costs 0.017% of signal.
- **The windows are limited by MC statistics.** Only 101 cosmic-ray tracks have a CRV cluster, so the optimum edges
  sit on the outermost of them. That is why the proposal adds a 10 ns margin (`--margin`). 38% of the
  cosmic-ray tracks in this selection have no CRV cluster at all, so no window can remove them.
- The stored variables reproduce the ConvAna veto exactly. On the signal and cosmic-ray samples, the
  `trk_crv_dt` (-50, 60) OR `trk_crv_raw_dt` (-25, 0) windows match the set 40 / 1040 split.

### Upstream (reflected) track window

Last update: 2026-10-27 (ConvAna set 43, `histograms_crv_dt_2026-10-27`, `python make_crv_plots.py --set 43`)

Set 43 has the multi-track veto applied but not the CRV or upstream vetoes. Efficiencies are relative to set 43,
so the cosmic efficiencies differ from the table above (179 cosmic-ray MC tracks, 111 with a CRV cluster).

| Veto | Signal loss | Cosmic eff. | All bkg. eff. | S/sqrt(B) | S/sqrt(B_cosmic) |
|---|---|---|---|---|---|
| ConvAna: CRV OR upstream (40, 110) | 0.238% | 0.374 | 0.9900 | 1.0026 | 1.631 |
| ConvAna upstream (40, 110) only | 0.082% | 0.905 | 0.9980 | 1.0002 | 1.050 |
| Raw CRV (0, 172.5) only | 0.152% | 0.380 | 0.9912 | 1.0029 | 1.620 |
| **Raw CRV (0, 172.5) OR upstream (40, 110)** | 0.234% | 0.341 | 0.9899 | 1.0028 | 1.709 |
| Raw CRV (10, 162.5) OR upstream (47.5, 105), no margin | 0.201% | 0.341 | 0.9900 | 1.0030 | 1.710 |
| Raw CRV (0, 172.5) OR upstream (37.5, 115), 10 ns margin | 0.245% | 0.341 | 0.9898 | 1.0027 | 1.709 |

- **Keep the ConvAna upstream window (40, 110) ns, and add it to the raw CRV window.** Together they veto every
  cosmic-ray track tagged by either, at 0.234% signal loss. That is less than the current ConvAna veto (0.238%),
  which misses 6 of these tracks. Cosmic-ray efficiency falls from 0.380 (raw CRV alone) to 0.341, and
  S/sqrt(B_cosmic) rises from 1.62 to 1.71.
- **The upstream veto adds 7 cosmic-ray tracks that have no CRV cluster.** 17 of the 179 tracks have an upstream
  track 43-103 ns earlier. 10 of them are also CRV-tagged, at raw CRV dt 90-155 ns, with t(upstream) - t(CRV)
  ≈ 45-55 ns (`crv_dt_2d_trk_crv_raw_dt_trk_upstream_dt`): the reflected cosmic mu+ of the late raw peak. The other
  7 are the only gain. The remaining 61 cosmic-ray tracks (~8 weighted events) have neither tag.
- **Don't widen the upstream window past 110 ns.** A signal track has an upstream track at ~1.2e-5/ns below
  110 ns and ~4e-5/ns above it (`crv_dt_trk_upstream_dt`). Without vetoes (set 45) the rate is flat at ~4e-5/ns.
  The step comes from the multi-track veto, which already removes most signal tracks with an upstream track less
  than ~110 ns earlier, so the upstream window is cheap only below that. A 10 ns margin costs 0.011% of signal,
  and the cosmic-ray tracks end at 103 ns. A single upstream window that also covers the CRV-tagged track at
  193 ns ((42.5, 195), the `trk_upstream_dt` optimum) costs 0.38% and isn't useful.
- The gain from the upstream veto rests on 7 MC tracks, and its cost depends on the multi-track veto. If that veto
  changes, rescan the upstream window.
- On set 43, the rebuilt ConvAna CRV and upstream vetoes reproduce the set 40 and 1040 yields exactly for signal
  and cosmic rays. So the closest-track `trk_upstream_dt` is enough to rebuild the any-track ConvAna veto.

Caveats:
- Signal dead time scales with the CRV accidental rate, which here is the Mix1BB rate. Expect roughly twice the
  dead time at 2BB.
- Each variable uses the single cluster closest in time. A window that excludes 0 could miss a cluster inside the
  window when another cluster is closer to the track. The scan is exact for the stored closest-cluster variables,
  but a veto in ConvAna should loop over all clusters.

## Inputs

The study reads the ConvAna track trees of set 40 (pass the CRV veto) and set 1040 (fail it, `kCRVVetoOffset`),
so the ConvAna veto is removed. Tracks without a CRV cluster have -10000 and are never vetoed. The variables
are filled in `Mu2eEvtAna/src/Mu2eEvtAna.cc` (`FillTree`):
- `trk_crv_dt`: track t0 - CRV cluster time propagated to the tracker, via the calorimeter front or the stopping
  target, whichever is closer in time. It uses the cluster closest in time, which is the same cluster that the
  ConvAna veto matches.
- `trk_crv_raw_dt`: track t0 - raw CRV cluster time, for the cluster closest in time.
- `trk_crv_cal_dt`, `trk_crv_st_dt`: track t0 - CRV cluster time propagated via the calorimeter front or via the
  stopping target, each for the cluster closest in that time. `trk_crv_dt` is the one with the smaller |dt|.
  These were added on 2026-10-26; older ConvAna files don't have them.
- `trk_upstream_dt`: track t0 - t0 of the upstream-going track closest in time, among those more than 25 ns earlier.
  ConvAna vetoes the track if any upstream-going track is 40-110 ns earlier (a cosmic-ray reflection in the
  solenoid field).
- `trk_crv_x`, `trk_crv_y`, `trk_crv_z`, `trk_crv_sector`: position and sector type of the CRV cluster used for
  `trk_crv_raw_dt`, to test a time-of-flight correction. These were added on 2026-10-26 and are in the
  `histograms_crv_dt_2026-10-27` files; the study doesn't use them yet.

Sets 40 and 1040 already have the upstream and multi-track vetoes applied (both set the `kUpstream` bit of the e+
ID), so they can't be used to study the upstream window. ConvAna set 43 is set 40 without the CRV, upstream, and
multi-track vetoes. With `--set 43`, the study applies the multi-track veto (`|trk_multitrk_dt| > 150 ns`),
emulates the ConvAna CRV and upstream vetoes, and checks that the emulation reproduces sets 40 and 1040. It then also
scans `trk_upstream_dt`, and searches for a raw CRV window OR an upstream window that vetoes every tagged cosmic ray.

The samples and their normalization come from `analysis/datasets.C` and `evaluate_norm` in
`analysis/optimize_selection.C`. The signal is `mumep`. The backgrounds are cosmic rays, RMC external and internal
(0n and 1n), and RPC external and internal.

## Running

```bash
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh
pyenv rootana 2.5.0
cd mumep_ana/studies/crv_deadtime
python make_crv_plots.py [--set 40|43] [--step 2.5] [--margin 10] [--fom cosmic|all] [--pmin 90 --pmax 93]
```

This takes about 1 min. The histogram directory follows `CONVANA_HIST_PATH`, as in `analysis/defaults.C`. Until
the files with set 43 and the per-hypothesis branches replace those in the default directory, use
`export CONVANA_HIST_PATH=/exp/mu2e/data/projects/run1a/mumep_ana/histograms_crv_dt_2026-10-27`
(`histograms_crv_dt_2026-10-26` has the same sets 40 and 1040, but no set 43).

Outputs in `figures/`, one set per variable `<var>` = `trk_crv_dt`, `trk_crv_raw_dt`, `trk_crv_cal_dt`, or
`trk_crv_st_dt`:
- `crv_dt_<var>`: time-difference distributions for signal, cosmic rays, and RMC + RPC, with the ConvAna and
  proposed windows. The left panel is wide (20 ns bins) and the right panel is zoomed (5 ns bins).
- `crv_dt_2d_<x>_<y>`: one time difference vs. another (closer propagated vs. raw, and via calorimeter vs. via
  target) for signal and cosmic rays. Cosmic rays vetoed by ConvAna are outlined.
- `crv_scan1d_<var>`: cosmic and all-background efficiency, signal loss, S/sqrt(B), and S/sqrt(B_cosmic), as each
  window edge is scanned with the other fixed at the proposal. The bands show the MC statistical uncertainty of
  the change from no veto.
- `crv_scan2d_<var>`: signal loss, background efficiency, S/sqrt(B), and S/sqrt(B_cosmic) as functions of
  both window edges.

With `--set 43` the figures go to `figures/set43/` and the summary to `outputs/crv_windows_set43.json`, and there
are also `trk_upstream_dt` figures and a raw CRV vs. upstream track time 2D plot.

`outputs/crv_windows.json` holds the windows, the split (calorimeter OR target) windows, and the efficiencies of
each veto option.

## Optimizing with the full selection

`analysis/optimize_selection.C` can include the CRV windows in the greedy cut optimization:
```bash
cd <scratch dir>  # writes cut_strings_mumep_crv.txt and roc_*mumep_crv.png into the current directory
root -q -b '<path>/mumep_ana/analysis/optimize_selection.C(false, 40, true)'
```
With `optimize_crv = true`, it adds set 1040 to remove the ConvAna CRV veto. It then adds `trk_crv_dt` and
`trk_crv_raw_dt` as veto (`kVeto`) variables, whose windows grow outward from the cosmic-ray peaks. The greedy
optimization steps in 0.5% signal efficiency, while the whole CRV veto costs ~0.15% of signal. So use this study's
scans to set the window, and the optimizer to check whether the other cuts change.
