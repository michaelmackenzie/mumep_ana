# RMC plots for the mu- -> e- sensitivity paper

These are the plots of the radiative muon capture (RMC) background for the mu- -> e- search. They cover the measured
RMC photon spectrum and the e+ control region that constrains the RMC rate and spectrum shape in the e- signal region.

## Environment

The Combine setup from the top-level README is needed. The constraint plots evaluate the post-fit model with the
Combine library:
```bash
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh
pyenv rootana 2.5.0
source /path/to/combine/HiggsAnalysis/CombinedLimit/env_standalone_mu2e.sh
export PYTHONPATH=/path/to/Mu2eEvtAna/combine:$PYTHONPATH   # mu2ecombine, as for analysis/combine
cd mumep_ana/studies/eminus_rmc
```

## RMC photon spectrum: `make_rmc_spectrum_plot.py`

```bash
python make_rmc_spectrum_plot.py [--result <RMCFits>/data/triumf_plestid_fits/data_1992_Al.root]
```

This plot shows the RMC photon spectrum on aluminum, as measured from the TRIUMF 1992 Al data in the RMCFits analysis
([TRIUMF_RMC_Fits](https://github.com/michaelmackenzie/TRIUMF_RMC_Fits), `scripts/fit_plestid.C`). The spectrum is
the sum of two phase-space shapes, x(1 - x)^n with x = k / kmax:

| Final state | kmax [MeV] | n | Fraction (k > 57 MeV) |
|---|---|---|---|
| No nucleon emission (27Mg) | 101.87 | 2 | 0.099 +0.019 -0.020 |
| One neutron emitted (26Mg) | 95.45 | 3.5 | 1 - the above |

The fractions are defined for photons above the 57 MeV fit threshold, so the spectrum is drawn above it, normalized
to unity there. The result is stored in the script. `--result` reads it from an RMCFits output file instead.

Output: `figures/rmc_photon_spectrum.{pdf,png}`, with linear and logarithmic (80-104 MeV) panels.

## e+ control region constraint: `make_rmc_constraint_plots.py`

```bash
python make_rmc_constraint_plots.py [--study-dir ../../analysis/combine/rmc_studies/mumem_75]
```

This script reads the outputs of `analysis/combine/rmc_constraint_study.py` and `rmc_sr_prescription.py`. It runs no
fits. Make those outputs first (see `analysis/combine/README.md`):
```bash
cd ../../analysis/combine
./rmc_constraint_study.py -s 75 --cr-range 94 105   # Combine fits on the background-only Asimov data set
./rmc_sr_prescription.py -s 75                      # signal-region prescription from the control region
```

The e- signal region (set 75, 100-110 MeV/c) sees RMC only from photons within ~2 MeV of the endpoint. The e+
control region (set 4075, 94-105 MeV/c) is dominated by RMC and measures three quantities:
- the RMC rate (`rmcNorm`, free)
- the spectrum power n (`RMCPower`: n = 2 +- 0.5 at theta = +-1, shape only)
- a lowered endpoint (`RMCEndpoint`: kmax - 0.5 MeV at theta = -1, shape only, one-sided)

All results are for the Run 1A exposure. The plots use control region fits only. The signal region is fit
separately: the control region result is extrapolated to it as an RMC rate and shape uncertainty
(`rmc_sr_prediction`).

Outputs in `figures/` (pdf and png):
- `rmc_spectrum_variations`: the extrapolation lever arm, i.e. the shapes of the RMC momentum distribution in the control region (e+) and the signal
  region (e-) for n = 1.5 and 2.5 and for kmax - 0.5 MeV, relative to the nominal.
- `rmc_cr_fit`: the control region Asimov data with the post-fit model (RMC rate and spectrum free) and its 68% band,
  from sampling the fit covariance. The RMC yield there is known to about +-6%.
- `rmc_cr_likelihood`: -2 Delta ln L profiles in the control region.
  - Spectrum power n, with the RMC rate free and fixed.
  - RMC rate, with the spectrum free and with n fixed.
  - Their 68%/95% CL contours.

  The curves with the other parameter free or fixed are nearly identical, so the rate and the shape are nearly
  uncorrelated.
- `rmc_cr_endpoint`: the endpoint shift in the control region. It is constrained when n is fixed (> -0.49 MeV at 68%
  CL, > -0.80 MeV at 95% CL one-sided). It is degenerate with n when both are free (the 2D contours).
- `rmc_sr_prediction`: the signal-region RMC predicted from the control region, without a simultaneous fit.
  - The yield / nominal, sampled from the control region likelihood with the endpoint fixed or lowered by up to
    0.5 MeV (the `rmc_sr_prescription.py` procedure, repeated with its seed). The script prints the quantiles next
    to `prescription.json` as a check. The dotted lines mark the recommended log-normal range (0.44-1.61).
  - The signal-region spectrum for the lowered endpoint, from the same samples. It uses the model's RMCPower and
    RMCEndpoint templates, morphed as Combine does. Its 68% band is split into the rate part (the yield factor) and
    the shape part.
  - The shape part is compared with the power-law shape of the recommendation (RMCPower alone, width 1.05). The two
    widths are similar (event-weighted rms 2.6% vs 2.4%). However, the lowered endpoint makes the e- shape harder
    in the signal region, so the CR shape band is offset from nominal: +7% to +54% at 102.0-102.2 MeV/c, vs -11% to
    +22% for the power law. In power-law units the offset is theta = +1.6 on top of a +-1.15 spread. The power law
    centered on nominal therefore does not cover it. The rate uncertainty dominates the yield in every bin.
- `rmc_sr_recommendation`: the signal-region external RMC with the final recommendation, as the signal-region fit
  will use it.
  - Rate: the asymmetric lnN (0.44/1.61).
  - Shape: the RMCPower variation, shape only, with a Gaussian prior of width 1.05 on theta (n = 2 +- 0.52).
  - Both are sampled with Combine's interpolations: vertical template morphing for the shape, and the asymmetric lnN
    interpolation for the rate.
  - The plot shows the nominal with the combined 68% band and the power-law shape up/down (theta = +-1.05). The ratio panel separates
    the rate-only and shape-only bands. The rate dominates; the shape adds about +-10-20% at the highest momenta.

The control region is Asimov data, so the plots show the expected constraint, not a measurement.
