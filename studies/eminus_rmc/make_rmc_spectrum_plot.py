#! /usr/bin/env python3
"""The radiative muon capture (RMC) photon spectrum on aluminum, as measured from the TRIUMF RMC data (1992 Al
dataset) in the RMCFits analysis (https://github.com/michaelmackenzie/TRIUMF_RMC_Fits, scripts/fit_plestid.C).

    python make_rmc_spectrum_plot.py [--result <RMCFits>/data/triumf_plestid_fits/data_1992_Al.root]

The measured spectrum is the sum of two phase-space shapes, x (1 - x)^n with x = k / kmax, for capture without
nucleon emission (27Mg final state, kmax = 101.87 MeV, n = 2) and with one emitted neutron (26Mg, kmax = 95.45 MeV,
n = 3.5). Their fractions are measured for photons above the k = 57 MeV threshold of the fit, so the spectrum is
drawn above it, normalized to unity there. The band is the 68% uncertainty of the fractions (they sum to one).

The result is stored below (RESULT). --result reads it from an RMCFits output file instead (the RooWorkspace
values br_0, br_unc_m_0, br_unc_p_0 of the no-emission component).

Outputs in figures/: rmc_photon_spectrum.{pdf,png} (linear and logarithmic panels).
"""

import argparse

import numpy as np

from rmc_tools import BLACK, BLUE, GRAY, ORANGE, save, setup_style

THRESHOLD = 57.   # MeV: the fractions are defined for photons above it
# RMCFits 1992 Al result (fit range 57-110 MeV), 2026-08-05
RESULT = {"fraction_0n": 0.098918, "minus": 0.019703, "plus": 0.018531}
COMPONENTS = (("0n", 101.87, 2.0, r"No nucleon emission ($k_\mathrm{max}$ = 101.87 MeV)", ORANGE, "--"),
              ("1n", 95.45, 3.5, r"One neutron emitted ($k_\mathrm{max}$ = 95.45 MeV)", BLUE, "-."))


def shape(k, kmax, n):
    """x (1 - x)^n, normalized to unity above THRESHOLD."""
    x = np.clip(k / kmax, 0., 1.)
    f = np.where((k >= THRESHOLD) & (k < kmax), x * (1. - x) ** n, 0.)
    grid = np.linspace(THRESHOLD, kmax, 20001)
    xg = grid / kmax
    return f / np.trapz(xg * (1. - xg) ** n, grid)


def read_result(path):
    import ROOT
    f = ROOT.TFile.Open(path)
    w = f.Get("workspace")
    out = {"fraction_0n": w.var("br_0").getVal(), "minus": w.var("br_unc_m_0").getVal(),
           "plus": w.var("br_unc_p_0").getVal()}
    f.Close()
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--result", default=None, help="RMCFits output file (data_1992_Al.root) to read the result from")
    args = parser.parse_args()
    result = read_result(args.result) if args.result else RESULT
    b, lo, hi = result["fraction_0n"], result["fraction_0n"] - result["minus"], result["fraction_0n"] + result["plus"]
    print(f"No-emission fraction (k > {THRESHOLD:g} MeV): {b:.4f} -{result['minus']:.4f} +{result['plus']:.4f}")

    k = np.linspace(THRESHOLD, 104., 4701)
    g = {name: shape(k, kmax, n) for name, kmax, n, *_ in COMPONENTS}
    total = lambda f: f * g["0n"] + (1. - f) * g["1n"]
    curves = np.array([total(lo), total(hi)])
    band_lo, band_hi = curves.min(axis=0), curves.max(axis=0)

    plt = setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for ax, log in zip(axes, (False, True)):
        ax.fill_between(k, band_lo, band_hi, color=GRAY, alpha=0.5, linewidth=0, label="Uncertainty (68%)")
        ax.plot(k, total(b), color=BLACK, linewidth=2, label="Total")
        for (name, kmax, n, label, color, style), f in zip(COMPONENTS, (b, 1. - b)):
            ax.plot(k, f * g[name], color=color, linestyle=style, linewidth=2, label=label)
        ax.set_xlabel(r"Photon energy $k$ [MeV]")
        ax.set_ylabel(r"$(1/N)\,dN/dk$ [MeV$^{-1}$], $k > 57$ MeV")
        ax.set_xlim(THRESHOLD, 104.)
        if log:
            ax.set_yscale("log")
            ax.set_ylim(1e-6, 0.1)
            ax.set_xlim(80., 104.)
        else:
            ax.set_ylim(0., 1.45 * total(b).max())
    axes[0].legend(loc="upper right", fontsize=11)
    axes[1].text(0.97, 0.95, r"RMC on $^{27}$Al, fit to TRIUMF data", transform=axes[1].transAxes, va="top",
                 ha="right", fontsize=12)
    axes[1].text(0.03, 0.05, f"No-emission fraction ($k > {THRESHOLD:g}$ MeV):\n"
                             f"${b:.3f}^{{+{result['plus']:.3f}}}_{{-{result['minus']:.3f}}}$",
                 transform=axes[1].transAxes, va="bottom", fontsize=12)
    fig.tight_layout()
    save(fig, "rmc_photon_spectrum")
    plt.close(fig)


if __name__ == "__main__":
    main()
