#! /usr/bin/env python3
"""Paper plots of the RMC constraint from the e+ control region for the mu- --> e- search, from the outputs of
analysis/combine/rmc_constraint_study.py and rmc_sr_prescription.py (no Combine fits are run here).

    python make_rmc_constraint_plots.py [--study-dir ../../analysis/combine/rmc_studies/mumem_75]

The radiative muon capture (RMC) background of the e- signal region (set S, 100-110 MeV/c) comes from photons within
~2 MeV of the endpoint. The e+ control region (set S + 4000, 94-105 MeV/c) is dominated by RMC and measures its rate
(rmcNorm) and the shape of the photon spectrum: the power n of x (1 - x)^n (RMCPower, theta = +-1 is n = 2 +- 0.5)
and a lowered endpoint (RMCEndpoint, theta = -1 is kmax - 0.5 MeV). All results are on the background-only Asimov
data set of the Run 1A exposure. Only control region fits are used: the signal region is fit separately, with the
RMC rate and shape uncertainties extrapolated from the control region (rmc_sr_prediction).

Outputs in figures/:
  rmc_spectrum_variations   how the photon spectrum variations change the e+ control region and e- signal region
                            RMC momentum shapes
  rmc_cr_fit                the control region (Asimov) data with the post-fit model and its 68% band
  rmc_cr_likelihood         -2 Delta ln L of the spectrum power n and of the RMC rate in the control region, and their
                            2D contours
  rmc_cr_endpoint           the endpoint shift in the control region, with the power fixed or free, and its 2D contours
                            with the power
  rmc_sr_prediction         the signal-region RMC yield and shape predicted from the control region (the prescription
                            for analyses that do not fit the control region simultaneously)
  rmc_sr_recommendation     the signal-region RMC with the recommended uncertainties (lnN rate, RMCPower shape with
                            its prior width), sampled with Combine's interpolations: the input to the signal-region fit
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np

from rmc_tools import (BLACK, BLUE, DEFAULT_STUDY, GRAY, GREEN, ORANGE, PINK, chi2_levels, combine_module, save,
                       setup_style)

study = combine_module("rmc_constraint_study")
prescription = combine_module("rmc_sr_prescription")
POWER, POWER_SHIFT, ENDPOINT_SHIFT = 2., study.POWER_SHIFT, study.ENDPOINT_SHIFT
to_n = lambda theta: POWER + POWER_SHIFT * np.asarray(theta)
to_shift = lambda theta: ENDPOINT_SHIFT * np.asarray(theta)
P_LABEL = r"$p$ [MeV/$c$]"
N_LABEL = r"Photon spectrum power $n$"
RATE_LABEL = r"RMC rate / nominal"
SHIFT_LABEL = r"Photon spectrum endpoint shift [MeV]"
DNLL = r"$-2\Delta\ln L$"
CR_PLOT_RANGE = (95., 102.)   # MeV/c: the control region momentum range shown (the fit uses all of it)


class Study:
    def __init__(self, path):
        self.path = Path(path)
        self.summary = json.load(open(self.path / "summary.json"))
        self.labels = self.summary["labels"]
        self.model = self.path / "free"
        channels = list(json.load(open(self.model / "model.json"))["channels"])
        self.sr = next(c for c in channels if c.startswith("set_"))
        self.cr = next(c for c in channels if c.startswith("cr"))
        self.sr_set, self.cr_set = self.sr.split("_")[1], str(4000 + int(self.sr.split("_")[1]))

    def scan(self, kind, param="RMCPower"):
        return study.read_scan(str(self.model / f"higgsCombine.{kind}.MultiDimFit.mH120.root"), param)

    def scan_2d(self, kind, px, py):
        return study.read_scan_2d(str(self.model / f"higgsCombine.{kind}.MultiDimFit.mH120.root"), px, py)

    def cr_range(self):
        edges = json.load(open(self.model / "model.json"))["channels"][self.cr]["observable"]["edges"]
        return edges[0], edges[-1]


def region_labels(s):
    lo, hi = s.cr_range()
    return {"sr": rf"Signal region ($e^-$, 100-110 MeV/$c$)",
            "cr": rf"$e^+$ control region ({lo:g}-{hi:g} MeV/$c$)"}


def stairs_band(ax, edges, lo, hi, **kw):
    ax.fill_between(edges, np.append(lo, lo[-1]), np.append(hi, hi[-1]), step="post", **kw)


#------------------------------------------------------------------------------------------------------
def plot_variations(s, plt):
    """Normalized RMC shape / nominal for the photon spectrum variations, in the control and signal regions."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    labels = region_labels(s)
    for ax, ch, region in ((axes[0], s.cr, "cr"), (axes[1], s.sr, "sr")):
        edges, power, _ = study.rmc_variation_decomposition(str(s.model), ch, "RMCPower")
        _, endpoint, _ = study.rmc_variation_decomposition(str(s.model), ch, "RMCEndpoint")
        nom = power["total"]["nominal"]
        ok = nom > 1e-3 * nom.max()
        ratio = lambda c: np.where(ok, (c / c.sum()) / np.where(ok, nom / nom.sum(), 1.), np.nan)
        ax.axhline(1., color=BLACK, linewidth=1.5, label=rf"Nominal ($n$ = 2, $k_\mathrm{{max}}$ = 101.87 MeV)")
        ax.stairs(ratio(power["total"]["up"]), edges, baseline=None, color=BLUE, linewidth=2,
                  label=rf"$n$ = {POWER + POWER_SHIFT:g}")
        ax.stairs(ratio(power["total"]["down"]), edges, baseline=None, color=ORANGE, linewidth=2, linestyle="--",
                  label=rf"$n$ = {POWER - POWER_SHIFT:g}")
        ax.stairs(ratio(endpoint["total"]["down"]), edges, baseline=None, color=GREEN, linewidth=2, linestyle="-.",
                  label=rf"$k_\mathrm{{max}}$ $-$ {ENDPOINT_SHIFT:g} MeV")
        last = edges[1:][ok].max()
        ax.set_xlim(edges[0], min(edges[-1], last + 0.5))
        ax.set_ylim(0., 2.6)
        ax.set_xlabel(P_LABEL)
        ax.set_ylabel("RMC shape / nominal")
        ax.text(0.03, 0.95, labels[region], transform=ax.transAxes, va="top", fontsize=12)
    axes[0].legend(loc="upper left", bbox_to_anchor=(0., 0.88), fontsize=11)
    fig.tight_layout()
    save(fig, "rmc_spectrum_variations", label_axes=axes)
    plt.close(fig)


#------------------------------------------------------------------------------------------------------
def plot_cr_fit(s, plt):
    """The control region (Asimov) data with the post-fit model stack, its 68% band, and the data / model ratio."""
    import ROOT
    info = json.load(open(s.model / "model.json"))["channels"][s.cr]
    edges = np.array(info["observable"]["edges"])
    processes = list(info["processes"])
    nominal, samples = study.postfit_yields(str(s.model), str(s.model / "multidimfit.cr_singles.root"), s.cr, edges,
                                            processes)
    f = ROOT.TFile.Open(str(s.model / "shapes" / f"{s.cr}.root"))
    dh = f.Get("w").data(f"{s.cr}_data_obs")
    data = np.array([(dh.get(i), dh.weight())[1] for i in range(dh.numEntries())])
    f.Close()
    centers = 0.5 * (edges[1:] + edges[:-1])
    zero = np.zeros(len(centers))
    rmc = [p for p in processes if p.startswith("rmc")]
    others = [p for p in processes if p not in rmc and p != "cosmic"]
    total = sum(nominal[p] for p in processes)
    t_lo, t_hi = np.percentile(sum(samples[p] for p in processes), [16., 84.], axis=0)
    stack = [("Other (RPC, $\\bar{p}$)", sum((nominal[p] for p in others), zero), GREEN),
             ("Cosmic rays", nominal.get("cosmic", zero), BLUE),
             ("RMC (external conversion)", nominal.get("rmc_ext", zero), ORANGE),
             ("RMC (internal conversion)", nominal.get("rmc_int", zero), PINK)]

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7, 6.4), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    cumulative = zero.copy()
    for label, y, color in stack:
        ax.stairs(cumulative + y, edges, baseline=cumulative, fill=True, color=color, alpha=0.85, label=label)
        cumulative = cumulative + y
    ax.stairs(total, edges, color=BLACK, linewidth=1.2)
    stairs_band(ax, edges, t_lo, t_hi, facecolor="none", edgecolor=BLACK, hatch="////", linewidth=0,
                label="Fit uncertainty (68%)")
    ax.errorbar(centers, data, yerr=np.sqrt(np.clip(data, 0., None)), fmt="o", color="black", markersize=3.5,
                label="Asimov data")
    width = edges[1] - edges[0]
    ax.set_ylabel(f"Events / {width:.1f} MeV/$c$")
    shown = (centers > CR_PLOT_RANGE[0]) & (centers < CR_PLOT_RANGE[1])
    ax.set_ylim(0., 1.55 * max(data[shown].max(), t_hi[shown].max()))
    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index("Asimov data")] + [i for i, l in enumerate(labels) if l != "Asimov data"][::-1]
    ax.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper right", fontsize=10.5)
    ok = shown & (total > 0)
    safe = np.where(ok, total, 1.)
    # the band runs one bin past the shown range (clipped by the axis), so the step fill reaches its edge
    band_ok = (edges[:-1] < CR_PLOT_RANGE[1] + 1e-6) & (edges[1:] > CR_PLOT_RANGE[0]) & (total > 0)
    safe = np.where(band_ok, total, 1.)
    stairs_band(rax, edges, np.where(band_ok, t_lo / safe, np.nan), np.where(band_ok, t_hi / safe, np.nan), color=GRAY,
                alpha=0.5, linewidth=0)
    rax.errorbar(centers[ok], data[ok] / total[ok], yerr=np.sqrt(np.clip(data[ok], 0., None)) / total[ok], fmt="o",
                 color="black", markersize=3)
    rax.axhline(1., color=BLACK, linewidth=1)
    rax.set_ylim(0.5, 1.5)
    rax.set_ylabel("Data / fit")
    rax.set_xlabel(r"$e^+$ " + P_LABEL)
    rax.set_xlim(*CR_PLOT_RANGE)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.06)
    save(fig, "rmc_cr_fit", label_axes=[ax])
    plt.close(fig)
    rmc_s = sum(samples[p] for p in rmc).sum(axis=1)
    rmc_n = sum(nominal[p] for p in rmc).sum()
    print(f"Control region: RMC {rmc_n:.1f} events, 68% band {np.percentile(rmc_s, 16) / rmc_n - 1:+.3f} / "
          f"{np.percentile(rmc_s, 84) / rmc_n - 1:+.3f}")


#------------------------------------------------------------------------------------------------------
def contour(ax, pts, tx, ty):
    """68% and 95% CL contours (2 parameters) of a 2D scan, with transformed axes."""
    x, y = tx(pts[:, 0]), ty(pts[:, 1])
    ax.tricontour(x, y, pts[:, 2], levels=[2.30], colors=[BLACK], linestyles=["-"], linewidths=2)
    ax.tricontour(x, y, pts[:, 2], levels=[5.99], colors=[BLUE], linestyles=["--"], linewidths=2)
    best = np.argmin(pts[:, 2])
    ax.plot(x[best], y[best], marker="+", color=ORANGE, markersize=14, markeredgewidth=2, linestyle="none",
            label="Best fit")
    ax.plot([], [], color=BLACK, linewidth=2, label="68% CL")
    ax.plot([], [], color=BLUE, linewidth=2, linestyle="--", label="95% CL")


def plot_likelihood(s, plt):
    """In the control region alone: -2 Delta ln L of n (RMC rate free / fixed), of the RMC rate (n free / fixed),
    and their 2D contours."""
    labels = region_labels(s)
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.2))
    ax = axes[0]
    for kind, color, style, label in (("cr_scan", BLACK, "-", r"RMC rate and $k_\mathrm{max}$ free"),
                                      ("cr_scan_power_rate_fixed", BLUE, "--", "RMC rate fixed")):
        sc = s.scan(kind)
        if sc is not None:
            ax.plot(to_n(sc[0]), sc[1], color=color, linestyle=style, linewidth=2, label=label)
    ax.set_xlabel(N_LABEL)
    ax.set_xlim(0., 4.)
    ax.legend(loc="upper left", bbox_to_anchor=(0.05, 1.), fontsize=10.5)
    ax.text(0.97, 0.03, labels["cr"], transform=ax.transAxes, ha="right", fontsize=10.5)
    ax = axes[1]
    for kind, color, style, label in (("cr_scan_norm", BLACK, "-", r"Spectrum ($n$, $k_\mathrm{max}$) free"),
                                      ("cr_scan_norm_power_fixed", BLUE, "--", r"$n$ fixed to 2")):
        sc = s.scan(kind, "rmcNorm")
        if sc is not None:
            ax.plot(sc[0], sc[1], color=color, linestyle=style, linewidth=2, label=label)
    ax.set_xlabel(RATE_LABEL)
    ax.set_xlim(0.75, 1.25)
    ax.legend(loc="upper center", fontsize=10.5)
    ax.text(0.03, 0.03, labels["cr"], transform=ax.transAxes, fontsize=10.5)
    for ax in axes[:2]:
        ax.set_ylabel(DNLL)
        ax.set_ylim(0., 8.)
        chi2_levels(ax)
    ax = axes[2]
    pts = s.scan_2d("cr_scan_2d", "RMCPower", "rmcNorm")
    if pts is not None and len(pts) > 10:
        contour(ax, pts, to_n, lambda v: v)
    ax.set_xlabel(N_LABEL)
    ax.set_ylabel(RATE_LABEL)
    ax.set_ylim(0.8, 1.25)
    ax.legend(loc="upper right", fontsize=10.5)
    ax.text(0.03, 0.03, labels["cr"], transform=ax.transAxes, fontsize=10.5)
    fig.tight_layout()
    save(fig, "rmc_cr_likelihood", label_axes=axes)
    plt.close(fig)


def plot_endpoint(s, plt):
    """The endpoint shift in the CR: -2 Delta ln L with the power fixed or free, and the 2D contours with the power."""
    labels = region_labels(s)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    ax = axes[0]
    for kind, color, style, label in (("cr_scan_endpoint", BLACK, "-", r"$n$ free"),
                                      ("cr_scan_endpoint_power_fixed", BLUE, "--", r"$n$ fixed to 2")):
        sc = s.scan(kind, "RMCEndpoint")
        if sc is not None:
            ax.plot(to_shift(sc[0]), sc[1], color=color, linestyle=style, linewidth=2, label=label)
    ax.set_xlim(-2.5, 0.)
    ax.set_ylim(0., 8.)
    ax.set_xlabel(SHIFT_LABEL)
    ax.set_ylabel(DNLL)
    chi2_levels(ax, ((1., "68% CL"), (2.71, "95% CL (one-sided)")))
    ax.legend(loc="upper right", fontsize=11)
    ax.text(0.97, 0.75, labels["cr"] + "\nRMC rate free", transform=ax.transAxes, ha="right", fontsize=10.5)
    ax = axes[1]
    pts = s.scan_2d("cr_scan_2d_spectrum", "RMCPower", "RMCEndpoint")
    if pts is not None and len(pts) > 10:
        contour(ax, pts, to_n, to_shift)
    ax.set_xlabel(N_LABEL)
    ax.set_ylabel(SHIFT_LABEL)
    ax.set_ylim(-2.5, 0.25)
    ax.legend(loc="lower right", fontsize=11)
    fig.tight_layout()
    save(fig, "rmc_cr_endpoint", label_axes=axes)
    plt.close(fig)


#------------------------------------------------------------------------------------------------------
def sr_templates(s):
    """The signal-region external RMC templates of the model: edges, nominal, RMCPower up/down, RMCEndpoint down (the
    endpoint cannot move up: its up template is the nominal)."""
    edges, power, _ = study.rmc_variation_decomposition(str(s.model), s.sr, "RMCPower")
    _, endpoint, _ = study.rmc_variation_decomposition(str(s.model), s.sr, "RMCEndpoint")
    return (edges, power["rmc_ext"]["nominal"], power["rmc_ext"]["up"], power["rmc_ext"]["down"],
            endpoint["rmc_ext"]["down"])


def shape_rms(shapes, nominal):
    """Event-weighted rms of the normalized shape / nominal around its mean (the rmc_sr_prescription.py measure)."""
    f = nominal / nominal.sum()
    ok = f > 1e-4
    ratio = (shapes / shapes.sum(axis=1, keepdims=True))[:, ok] / f[ok]
    return float(np.sqrt(np.mean(np.sum(f[ok] * (ratio - ratio.mean(axis=0)) ** 2, axis=1))))


def plot_sr_prediction(s, plt, n_samples):
    """The signal-region RMC predicted from the control region likelihood (the rmc_sr_prescription.py procedure,
    repeated with its seed): the yield / nominal with the endpoint fixed or lowered, and the spectrum for the lowered
    endpoint split into its rate (yield factor) and shape parts, compared with the shape part of the recommendation
    (the RMCPower variation alone with its prior width). The shapes use the model templates and Combine's vertical
    morphing (each nuisance's shift added to the nominal)."""
    path = s.path / "sr_prescription" / "prescription.json"
    if not path.exists():
        print(f"WARNING: no {path}: run analysis/combine/rmc_sr_prescription.py; skipping rmc_sr_prediction")
        return
    pre = json.load(open(path))
    inputs = pre["inputs"]
    ratio = inputs["sr_cr_ratio"]
    prior_max = pre["endpoint_prior_MeV"]
    kappa = pre["recommendation"]["lnN"]
    width = pre["recommendation"]["rmc_power_theta_sigma"]
    scan = prescription.read_scan_2d(str(s.model / "higgsCombine.cr_scan_2d_spectrum.MultiDimFit.mH120.root"))
    rng = np.random.default_rng(1)   # as rmc_sr_prescription.py: the same samples in the same order
    posterior = {}
    for prior in dict.fromkeys((0., prior_max)):
        tp, te, _ = prescription.sample_spectrum(scan, prior, n_samples, rng)
        spectrum = prescription.interp(tp, ratio[0], ratio[1]) * np.power(ratio[2], -te)
        rate = np.exp(rng.normal(0., inputs["rate_sigma"], n_samples) + rng.normal(0., inputs["external_fraction_sigma"],
                                                                                     n_samples))
        posterior[prior] = (tp, te, spectrum * rate)
        q = np.percentile(spectrum * rate, [16., 50., 84.])
        ref = pre["scenarios"][str(prior)]["total"]
        print(f"Signal-region RMC yield factor, endpoint prior {prior:g} MeV: {q[0]:.3f} / {q[1]:.3f} / {q[2]:.3f}"
              f" (prescription.json: {ref[0]:.3f} / {ref[1]:.3f} / {ref[2]:.3f})")

    # signal-region spectra of the lowered-endpoint posterior: yield factor x shape (power and endpoint morphed)
    edges, nominal, p_up, p_down, e_down = sr_templates(s)
    yield0 = nominal.sum()
    tp, te, factor = posterior[prior_max]
    sub = slice(0, min(50000, n_samples))
    tp, te, factor = tp[sub], te[sub], factor[sub]
    t = lambda x: np.asarray(x, dtype=float)[:, None]
    shapes = nominal + 0.5 * t(tp) * ((p_up - p_down) + smooth_step(t(tp)) * (p_up + p_down - 2. * nominal)) \
                     + 0.5 * t(te) * ((nominal - e_down) + smooth_step(t(te)) * (e_down - nominal))
    shapes = np.clip(shapes, 0., None)
    shapes *= yield0 / shapes.sum(axis=1, keepdims=True)
    rec = morph(np.random.default_rng(2).normal(0., width, len(factor)), nominal, p_up, p_down)
    rec *= yield0 / rec.sum(axis=1, keepdims=True)
    q = lambda x: np.percentile(x, [16., 84.], axis=0)
    bands = {"rate": q(factor[:, None] * nominal[None, :]), "shape": q(shapes), "total": q(factor[:, None] * shapes),
             "recommendation": q(rec)}
    rms_post, rms_rec = shape_rms(shapes, nominal), shape_rms(rec, nominal)
    print(f"Signal-region shape (event-weighted rms of shape / nominal): CR posterior with the endpoint lowered by "
          f"<= {prior_max:g} MeV {100 * rms_post:.1f}%, recommendation (RMCPower, width {width:.2f}) {100 * rms_rec:.1f}%")
    safe = np.where(nominal > 0, nominal, 1.)
    for i in (0, len(nominal) // 10, len(nominal) // 5):
        print(f"  {edges[i]:.1f}-{edges[i + 1]:.1f} MeV/c: shape band {bands['shape'][0][i] / safe[i]:.3f}-"
              f"{bands['shape'][1][i] / safe[i]:.3f}, recommendation {bands['recommendation'][0][i] / safe[i]:.3f}-"
              f"{bands['recommendation'][1][i] / safe[i]:.3f}")

    fig = plt.figure(figsize=(14, 7.6))
    grid = fig.add_gridspec(3, 2, height_ratios=[2.6, 1.1, 1.1], hspace=0.07, wspace=0.25)
    ax = fig.add_subplot(grid[:, 0])
    bins = np.linspace(0., 2.5, 101)
    for (prior, (_, _, f)), color, style in zip(posterior.items(), (BLUE, BLACK), ("--", "-")):
        label = r"Endpoint fixed" if prior == 0 else rf"Endpoint lowered by $\leq$ {prior:g} MeV"
        ax.hist(f, bins=bins, histtype="step", linewidth=2, density=True, color=color, linestyle=style, label=label)
    for k in kappa:
        ax.axvline(k, ymax=0.7, color=ORANGE, linestyle=":", linewidth=2)
    ax.plot([], [], color=ORANGE, linestyle=":", linewidth=2, label=f"Uncertainty: {kappa[0]:.2f} - {kappa[1]:.2f}")
    ax.set_xlabel("Signal-region RMC yield / nominal")
    ax.set_ylabel("Probability density")
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0., 1.5 * ax.get_ylim()[1])
    ax.legend(loc="upper right", fontsize=11)

    ax = fig.add_subplot(grid[0, 1])
    rax = fig.add_subplot(grid[1, 1], sharex=ax)
    sax = fig.add_subplot(grid[2, 1], sharex=ax)
    stairs_band(ax, edges, *bands["total"], color=GRAY, alpha=0.55, linewidth=0, label="Rate and shape (68%)")
    ax.stairs(nominal, edges, color=BLACK, linewidth=2, label=rf"Nominal, {yield0:.2f} events")
    ax.set_yscale("log")
    ok = nominal > 1e-3 * nominal.max()
    ax.set_xlim(edges[0], min(edges[-1], edges[1:][ok].max()))
    ax.set_ylim(0.3 * nominal[ok].min(), 15. * nominal.max())
    ax.set_ylabel(f"Events / {edges[1] - edges[0]:.1f} MeV/$c$")
    ax.legend(loc="upper right", fontsize=10.5)
    ax.text(0.03, 0.05, "External RMC, $e^-$ signal region\n" + rf"CR likelihood, endpoint lowered by $\leq$ {prior_max:g} MeV",
            transform=ax.transAxes, fontsize=10.5)
    ax.tick_params(labelbottom=False)
    stairs_band(rax, edges, bands["total"][0] / safe, bands["total"][1] / safe, color=GRAY, alpha=0.55, linewidth=0,
                label="Rate and shape")
    curves = (("rate", rax, BLUE, "--", "Rate only"), ("shape", sax, ORANGE, "-", "Shape only (CR likelihood)"),
              ("recommendation", sax, GREEN, ":", rf"Power law, $n$ = 2 $\pm$ {POWER_SHIFT * width:.2f}"))
    for key, axis, color, style, label in curves:
        for i, side in enumerate(bands[key]):
            axis.stairs(side / safe, edges, baseline=None, color=color, linestyle=style, linewidth=2,
                        label=label if i == 0 else None)
    for axis, (lo, hi) in ((rax, (0., 2.6)), (sax, (0.7, 1.8))):
        axis.axhline(1., color=BLACK, linewidth=1)
        axis.set_ylim(lo, hi)
        axis.set_ylabel("Ratio")
        axis.legend(loc="upper left", ncol=2, fontsize=9.5, columnspacing=1.)
    rax.tick_params(labelbottom=False)
    sax.set_xlabel(r"$e^-$ " + P_LABEL)
    save(fig, "rmc_sr_prediction", label_axes=[fig.axes[0], ax])
    plt.close(fig)


#------------------------------------------------------------------------------------------------------
def smooth_step(x):
    """Combine's smooth step (CMSHistFunc, CombineMathFuncs): +-1 beyond |x| = 1, a quintic in between."""
    x = np.asarray(x, dtype=float)
    xc = np.clip(x, -1., 1.)
    return np.where(np.abs(x) >= 1., np.sign(x), 0.125 * xc * (xc**2 * (3. * xc**2 - 10.) + 15.))


def morph(theta, nominal, up, down):
    """Combine's vertical template morphing of a shape nuisance (CMSHistFunc): quadratic for |theta| < 1, linear
    beyond, negative bins cropped. theta: (n,) -> (n, bins)."""
    t = np.asarray(theta, dtype=float)[:, None]
    out = nominal + 0.5 * t * ((up - down) + smooth_step(t) * (up + down - 2. * nominal))
    return np.clip(out, 0., None)


def asym_lnN(nu, kappa_lo, kappa_hi):
    """Combine's asymmetric lnN factor (ProcessNormalization): kappa_hi^nu / kappa_lo^-nu beyond |nu| = 0.5, with a
    smooth interpolation of log(kappa) in between."""
    nu = np.asarray(nu, dtype=float)
    lhi, llo = np.log(kappa_hi), -np.log(kappa_lo)
    alpha = smooth_step(2. * nu)
    log_kappa = np.where(np.abs(nu) >= 0.5, np.where(nu >= 0, lhi, llo), 0.5 * (lhi + llo) + 0.5 * alpha * (lhi - llo))
    return np.exp(log_kappa * nu)


def plot_sr_recommendation(s, plt, n_samples):
    """The signal-region external RMC with the recommended uncertainties: the asymmetric lnN on the rate and the
    RMCPower variation (shape only) with a Gaussian prior of the recommended width, sampled as Combine would apply
    them, with the 68% bands of the rate alone, the shape alone, and both."""
    path = s.path / "sr_prescription" / "prescription.json"
    if not path.exists():
        print(f"WARNING: no {path}: run analysis/combine/rmc_sr_prescription.py; skipping rmc_sr_recommendation")
        return
    rec = json.load(open(path))["recommendation"]
    kappa_lo, kappa_hi = rec["lnN"]
    width = rec["rmc_power_theta_sigma"]
    edges, curves, _ = study.rmc_variation_decomposition(str(s.model), s.sr, "RMCPower")
    nominal, up, down = (curves["rmc_ext"][k] for k in ("nominal", "up", "down"))
    yield0 = nominal.sum()
    rng = np.random.default_rng(1)
    theta = rng.normal(0., width, n_samples)
    rate = asym_lnN(rng.normal(0., 1., n_samples), kappa_lo, kappa_hi)
    shapes = morph(theta, nominal, up, down)
    shapes *= yield0 / shapes.sum(axis=1, keepdims=True)   # shape only: the morphed template keeps the nominal yield
    q = lambda x: np.percentile(x, [16., 84.], axis=0)
    bands = {"rate": q(rate[:, None] * nominal[None, :]), "shape": q(shapes), "total": q(rate[:, None] * shapes)}
    total_yield = np.percentile(rate * shapes.sum(axis=1), [16., 50., 84.])
    print(f"Signal-region external RMC: nominal {yield0:.3f} events, 16/50/84% {total_yield[0]:.3f} / "
          f"{total_yield[1]:.3f} / {total_yield[2]:.3f}")
    edge_n = lambda t: morph([t], nominal, up, down)[0] * yield0 / morph([t], nominal, up, down)[0].sum()

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7, 6.4), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]})
    stairs_band(ax, edges, *bands["total"], color=GRAY, alpha=0.55, linewidth=0,
                label="Rate and shape uncertainty (68%)")
    ax.stairs(nominal, edges, color=BLACK, linewidth=2, label="Nominal")
    ax.stairs(edge_n(width), edges, baseline=None, color=BLUE, linewidth=1.8, linestyle="--",
              label="Power-law shape up")
    ax.stairs(edge_n(-width), edges, baseline=None, color=ORANGE, linewidth=1.8, linestyle="-.",
              label="Power-law shape down")
    ax.set_yscale("log")
    ok = nominal > 1e-3 * nominal.max()
    hi = min(edges[-1], edges[1:][ok].max())
    ax.set_xlim(edges[0], hi)
    ax.set_ylim(0.3 * nominal[ok].min(), 15. * nominal.max())
    ax.set_ylabel(f"Events / {edges[1] - edges[0]:.1f} MeV/$c$")
    ax.legend(loc="upper right", fontsize=10.5)
    safe = np.where(nominal > 0, nominal, 1.)
    stairs_band(rax, edges, bands["total"][0] / safe, bands["total"][1] / safe, color=GRAY, alpha=0.55, linewidth=0,
                label="Rate and shape")
    for key, color, style, label in (("rate", BLUE, "--", "Rate only"), ("shape", ORANGE, "-", "Shape only")):
        for i, side in enumerate(bands[key]):
            rax.stairs(side / safe, edges, baseline=None, color=color, linestyle=style, linewidth=1.8,
                       label=label if i == 0 else None)
    rax.axhline(1., color=BLACK, linewidth=1)
    rax.set_ylim(0., 2.6)
    rax.set_ylabel("Ratio to nominal")
    rax.set_xlabel(r"$e^-$ " + P_LABEL)
    rax.legend(loc="upper left", ncol=3, fontsize=10)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.06)
    save(fig, "rmc_sr_recommendation", label_axes=[ax])
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--study-dir", default=str(DEFAULT_STUDY),
                        help="rmc_constraint_study.py output directory (default: %(default)s)")
    parser.add_argument("--samples", type=int, default=200000,
                        help="posterior samples, as rmc_sr_prescription.py (default: %(default)s)")
    args = parser.parse_args()
    os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    import ROOT
    ROOT.gROOT.SetBatch(True)
    s = Study(args.study_dir)
    plt = setup_style()
    plot_variations(s, plt)
    plot_cr_fit(s, plt)
    plot_likelihood(s, plt)
    plot_endpoint(s, plt)
    plot_sr_prediction(s, plt, args.samples)
    plot_sr_recommendation(s, plt, args.samples)


if __name__ == "__main__":
    main()
