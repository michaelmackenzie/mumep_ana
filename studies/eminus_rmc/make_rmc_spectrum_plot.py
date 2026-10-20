#! /usr/bin/env python3
"""The radiative muon capture (RMC) photon spectrum on aluminum, from the phase-space forward-convolution fit to the
TRIUMF 1992 Al data in the RMC review paper (RMC-Review, analysis_triumf.tex, Table tab:triumf_plestid_results),
normalized with the TRIUMF RMC branching fraction.

    python make_rmc_spectrum_plot.py

The spectrum is the sum of two phase-space shapes, x (1 - x)^n with x = k / kmax: capture without nucleon emission
(27Mg final state, kmax = 101.87 MeV, n = 2) and with one emitted neutron (26Mg, kmax = 95.45 MeV, n = 3.5). The
paper gives their fractions of the spectrum above k = 57 MeV:
    B(no emission) = (9.9 +- 1.9 (stat) +4.4 -3.4 (syst)) %,   B(one neutron) = 1 - B(no emission)
and TRIUMF measured the RMC branching fraction above 57 MeV, relative to muon capture, as (1.40 +- 0.11) x 1e-5. The
spectrum is drawn above 57 MeV as dN_gamma/dk per muon capture. The band is the 68% range of the spectrum, sampling the
fraction (stat and syst in quadrature, each side as a half-Gaussian) and the branching fraction (Gaussian).

Outputs in figures/: rmc_photon_spectrum.{pdf,png} (linear, 57-104 MeV) and rmc_photon_spectrum_log.{pdf,png}
(logarithmic, 80-104 MeV).
"""

import numpy as np

from rmc_tools import BLACK, BLUE, GRAY, ORANGE, save, setup_style

THRESHOLD = 57.   # MeV: the fractions and the branching fraction are defined for photons above it
# RMC review paper, phase-space fit to the TRIUMF 1992 Al data: the no-emission fraction above THRESHOLD
FRACTION_0N = {"value": 0.099, "stat": 0.019, "syst_up": 0.044, "syst_down": 0.034}
# TRIUMF RMC branching fraction above THRESHOLD, relative to muon capture
BR_RMC = {"value": 1.40e-5, "error": 0.11e-5}
COMPONENTS = (("0n", 101.87, 2.0, "No nucleon emission", ORANGE, "--"),
              ("1n", 95.45, 3.5, "One neutron emitted", BLUE, "-."))
N_SAMPLES = 100000


def shape(k, kmax, n):
    """x (1 - x)^n, normalized to unity above THRESHOLD."""
    x = np.clip(k / kmax, 0., 1.)
    f = np.where((k >= THRESHOLD) & (k < kmax), x * (1. - x) ** n, 0.)
    grid = np.linspace(THRESHOLD, kmax, 20001)
    xg = grid / kmax
    return f / np.trapz(xg * (1. - xg) ** n, grid)


def sample_fraction(rng, n):
    """The no-emission fraction with stat and syst in quadrature, a half-Gaussian of each width on each side."""
    f = FRACTION_0N
    up, down = np.hypot(f["stat"], f["syst_up"]), np.hypot(f["stat"], f["syst_down"])
    z = rng.normal(0., 1., n)
    return np.clip(f["value"] + np.where(z >= 0., up, down) * z, 0., 1.), up, down


def main():
    rng = np.random.default_rng(1)
    b = FRACTION_0N["value"]
    fractions, up, down = sample_fraction(rng, N_SAMPLES)
    rates = np.clip(rng.normal(BR_RMC["value"], BR_RMC["error"], N_SAMPLES), 0., None)
    print(f"No-emission fraction (k > {THRESHOLD:g} MeV): {b:.3f} +{up:.3f} -{down:.3f} (stat + syst); "
          f"B(RMC, k > {THRESHOLD:g} MeV) = ({BR_RMC['value'] * 1e5:.2f} +- {BR_RMC['error'] * 1e5:.2f}) x 1e-5")

    k = np.linspace(THRESHOLD, 104., 2351)
    g = {name: shape(k, kmax, n) for name, kmax, n, *_ in COMPONENTS}
    total = lambda f, r: r * (f * g["0n"] + (1. - f) * g["1n"])
    samples = rates[:, None] * (fractions[:, None] * g["0n"][None, :] + (1. - fractions[:, None]) * g["1n"][None, :])
    band_lo, band_hi = np.percentile(samples, [16., 84.], axis=0)
    nominal = total(b, BR_RMC["value"])

    plt = setup_style()
    for log in (False, True):
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.fill_between(k, band_lo, band_hi, color=GRAY, alpha=0.5, linewidth=0,
                        label="Uncertainty (68%, stat. + syst.)")
        ax.plot(k, nominal, color=BLACK, linewidth=2, label="Total")
        for (name, kmax, n, label, color, style), f in zip(COMPONENTS, (b, 1. - b)):
            ax.plot(k, BR_RMC["value"] * f * g[name], color=color, linestyle=style, linewidth=2, label=label)
        ax.set_xlabel(r"Photon energy $k$ [MeV]")
        ax.set_ylabel(r"$dN_\gamma/dk$ per muon capture [MeV$^{-1}$]")
        if log:
            ax.set_yscale("log")
            ax.set_ylim(1e-11, 1e-6)
            ax.set_xlim(80., 104.)
            ax.legend(loc="lower left", fontsize=11)
        else:
            ax.set_xlim(THRESHOLD, 104.)
            ax.set_ylim(0., 1.45 * nominal.max())
            ax.legend(loc="upper right", fontsize=11)
        ax.text(0.97, 0.95 if log else 0.6, r"RMC on $^{27}$Al, fit to TRIUMF data", transform=ax.transAxes, va="top",
                ha="right", fontsize=12)
        fig.tight_layout()
        save(fig, "rmc_photon_spectrum" + ("_log" if log else ""))
        plt.close(fig)

if __name__ == "__main__":
    main()
