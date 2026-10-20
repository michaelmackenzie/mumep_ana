#! /usr/bin/env python3
"""How well the e+ (charge-inverted) control region measures the RMC of the mu- --> e- search: its rate in the high
momentum region and the shape of its photon spectrum.

    ./rmc_constraint_study.py [-s 75] [--cr-range 94 105] [-t <tag>] [--jobs 4]

The model (build_model.py with the control region) measures both from the control region data:
  rmcNorm   the RMC rate: a free scale of every RMC process in every channel (no rmcN or lumi on the RMC)
  RMCPower     the photon spectrum power, shape only (Mu2eEvtAna variations of x(1-x)^n, n = 2 +- 0.5, on the 0n
               samples, normalized to the nominal yield of each process in each channel): theta = +-1 is the example shift
  RMCEndpoint  the photon spectrum endpoint, shape only: theta = -1 lowers kmax = 101.866 MeV by 0.5 MeV; the endpoint
               cannot move higher, so theta is limited to [-5, 0] and the data give a one-sided limit
  intN         the internal conversion fraction keeps its 4.5% prior (the internal/external ratio)

Models (in <out dir>/<model>/), with the signal region (set S, 100-110 MeV/c) and the e+ control region (set S + 4000,
--cr-range), and channel masks to fit either one alone:
  free      RMCPower and RMCEndpoint with flat priors
  prior     RMCPower and RMCEndpoint with unit-Gaussian priors (the example shifts taken as 1 sigma)

Combine on the background-only Asimov data set (r = 0), for the signal region alone (sr), the control region alone
(cr, r frozen), and both (sr_cr):
  - MultiDimFit --algo singles: the profile-likelihood interval of RMCPower, and the fit covariance (rmcNorm, intN)
  - likelihood scans of RMCPower and RMCEndpoint (sr, cr, sr_cr), and in the control region of RMCPower with the rate
    fixed, of RMCEndpoint with the power fixed, of rmcNorm with the power free or fixed, of rmcNorm vs RMCPower and of
    RMCEndpoint vs RMCPower (2D)
  - AsymptoticLimits (90% CL): sr and sr_cr with the spectrum free, with the unit priors, and frozen
Outputs: summary.txt and summary.json; rmc_power_scan.png (RMCPower scans, n on the top axis); rmc_cr_fit.png (the
control region alone: RMCPower and rmcNorm likelihood curves, their 2D contours, and the RMC shape for n = 1.5, 2,
2.5); rmc_cr_data_model.png (the control region data with the fitted model and its 68% band, from sampling the fit
covariance); rmc_cr_rmc_model.png (the fitted RMC model alone with its band); rmc_power_templates.png (the nominal and
+-1 templates in each channel); rmc_cr_endpoint.png (the RMCEndpoint likelihood curves, its 2D contours with the
power, and the RMC shape with the lowered endpoint); rmc_endpoint_templates.png.
"""

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
POWER_SHIFT = 0.5      # the RMCPower variation: n = 2 +- 0.5 at theta = +-1
ENDPOINT_SHIFT = 0.5   # the RMCEndpoint variation: kmax lowered by 0.5 MeV at theta = -1
NUISANCES = ["RMCPower", "RMCEndpoint", "rmcNorm", "intN"]


def run(command, cwd, log, dry_run=False):
    print(f"(cd {cwd} && {command})")
    if dry_run:
        return 0, ""
    with open(log, "w") as f:
        result = subprocess.run(command, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True)
        f.write(command + "\n\n" + result.stdout)
    return result.returncode, result.stdout


def parse_singles(text):
    """{parameter: (best fit, -err, +err)} from MultiDimFit --algo singles."""
    out = {}
    for m in re.finditer(r"(?m)^\s*(\w+)\s*:\s*([-+]?[\d.eE+-]+)\s+([-+][\d.eE+-]+)/([-+][\d.eE+-]+)\s*\(68%\)", text):
        out[m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    return out


def parse_limits(text):
    out = {}
    m = re.findall(r"Observed Limit:\s*r\s*<\s*([-\d.eE+]+)", text)
    out["median"] = float(m[-1]) if m else None
    for q in ("16.0", "50.0", "84.0"):
        m = re.findall(rf"Expected\s*{re.escape(q)}%:\s*r\s*<\s*([-\d.eE+]+)", text)
        out[q] = float(m[-1]) if m else None
    return out


def read_scan(path, param="RMCPower"):
    """(x, -2 Delta ln L) of a 1D MultiDimFit grid scan."""
    import ROOT
    if not os.path.exists(path):
        return None
    f = ROOT.TFile.Open(path)
    t = f.Get("limit") if f else None
    if not t:
        return None
    points = sorted((getattr(e, param), 2. * e.deltaNLL) for i, e in enumerate(t) if i > 0)
    f.Close()
    x, y = np.array([p[0] for p in points]), np.array([p[1] for p in points])
    return x, y - y.min()


def read_scan_2d(path, px, py):
    """(x, y, -2 Delta ln L) points of a 2D MultiDimFit grid scan."""
    import ROOT
    if not os.path.exists(path):
        return None
    f = ROOT.TFile.Open(path)
    t = f.Get("limit") if f else None
    if not t:
        return None
    pts = np.array([(getattr(e, px), getattr(e, py), 2. * e.deltaNLL) for i, e in enumerate(t) if i > 0])
    f.Close()
    pts[:, 2] -= pts[:, 2].min()
    return pts


def interval(scan, level=1.):
    """(best, lo, hi) of a 1D scan: the minimum and the crossings of the level, by linear interpolation."""
    if scan is None:
        return None
    x, y = scan
    i0 = int(np.argmin(y))
    lo = hi = None
    for i in range(i0, 0, -1):
        if y[i - 1] >= level > y[i]:
            lo = x[i - 1] + (level - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])
            break
    for i in range(i0, len(x) - 1):
        if y[i] < level <= y[i + 1]:
            hi = x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
            break
    return x[i0], lo, hi


def rmc_variation_decomposition(model_dir, channel, variation="RMCPower"):
    """Per RMC process with RMCPower in a channel, and for their sum: nominal/up/down per bin, the edges, the yield
    ratios, and the shape ratios (each normalized to its own yield) to the nominal."""
    import ROOT
    info = json.load(open(os.path.join(model_dir, "model.json")))["channels"][channel]
    edges = np.array(info["observable"]["edges"])
    f = ROOT.TFile.Open(os.path.join(model_dir, "shapes", f"{channel}.root"))
    w = f.Get("w")
    out = {}
    for p, v in info["processes"].items():
        if variation not in v["nuisances"]:
            continue
        curves = {}
        for key, suffix in (("nominal", ""), ("up", f"_{variation}Up"), ("down", f"_{variation}Down")):
            dh = w.data(f"{channel}_{p}{suffix}")
            curves[key] = np.array([(dh.get(i), dh.weight())[1] for i in range(dh.numEntries())])
        out[p] = curves
    f.Close()
    out["total"] = {k: sum(out[p][k] for p in list(out)) for k in ("nominal", "up", "down")}
    summary = {}
    for p, c in out.items():
        n = c["nominal"]
        ok = n > 1e-3 * n.max()
        entry = {"yield": float(n.sum())}
        for key in ("up", "down"):
            ratio = (c[key] / c[key].sum()) / (n / n.sum())
            entry[f"{key}_yield_ratio"] = float(c[key].sum() / n.sum())
            entry[f"{key}_shape_ratio_range"] = [float(ratio[ok].min()), float(ratio[ok].max())]
        summary[p] = entry
    return edges, out, summary


def plot_cr_floating(power, norm, scan2d, edges, curves, path, cr_label):
    """Control region alone: the RMCPower and rmcNorm likelihood curves, their 2D contours, and the RMC shape (each
    template normalized) for n = 1.5, 2, 2.5 relative to n = 2."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.2))
    styles = {"floating": ("#1F1F1F", "-"), "fixed": ("#2a78d6", "--")}
    ax = axes[0]
    for key, label in (("floating", "RMC rate floating (rmcNorm free)"),
                       ("fixed", "RMC rate fixed (rmcNorm = 1)")):
        if power.get(key) is not None:
            ax.plot(*power[key], color=styles[key][0], linestyle=styles[key][1], linewidth=2, label=label)
    ax.set_xlabel("RMCPower $\\theta$ (unconstrained)", fontsize=12)
    top = ax.secondary_xaxis("top", functions=(lambda t: 2. + POWER_SHIFT * t, lambda n: (n - 2.) / POWER_SHIFT))
    top.set_xlabel("RMC spectrum power n", fontsize=11)
    ax = axes[1]
    for key, label in (("floating", "RMCPower free"), ("fixed", "RMCPower fixed (n = 2)")):
        if norm.get(key) is not None:
            ax.plot(*norm[key], color=styles[key][0], linestyle=styles[key][1], linewidth=2, label=label)
    ax.set_xlabel("RMC rate scale factor rmcNorm (unconstrained)", fontsize=12)
    for ax in axes[:2]:
        for level, text in ((1., "68% CL"), (3.84, "95% CL")):
            ax.axhline(level, color="#5C5C5C", linestyle=":", linewidth=1)
        ax.set_ylabel(r"$-2\Delta\ln L$", fontsize=12)
        ax.set_ylim(0., 10.)
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False, fontsize=9, loc="upper center")
    ax = axes[2]
    if scan2d is not None and len(scan2d) > 10:
        ax.tricontour(scan2d[:, 0], scan2d[:, 1], scan2d[:, 2], levels=[2.30, 5.99], colors=["#1F1F1F", "#2a78d6"],
                      linestyles=["-", "--"])
        ax.plot([0.], [1.], marker="+", color="#eb6834", markersize=12)
        ax.plot([], [], color="#1F1F1F", label="68% CL")
        ax.plot([], [], color="#2a78d6", linestyle="--", label="95% CL")
        ax.legend(frameon=False, fontsize=9)
    ax.set_xlabel("RMCPower $\\theta$", fontsize=12)
    ax.set_ylabel("rmcNorm", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax = axes[3]
    n = curves["nominal"]
    for key, color, label in (("down", "#eb6834", f"n = {2 - POWER_SHIFT:g}"), ("nominal", "#1F1F1F", "n = 2"),
                              ("up", "#2a78d6", f"n = {2 + POWER_SHIFT:g}")):
        c = curves[key]
        ax.stairs((c / c.sum()) / (n / n.sum()), edges, color=color, linewidth=2, label=label)
    ax.set_xlabel("p [MeV/c]", fontsize=12)
    ax.set_ylabel("RMC shape / nominal (each normalized)", fontsize=12)
    ax.set_xlim(edges[0], edges[0] + 0.75 * (edges[-1] - edges[0]))
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=9)
    fig.suptitle(f"{cr_label}: RMC rate and spectrum power", fontsize=13)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_cr_endpoint(endpoint, scan2d, edges, curves, path, cr_label, labels):
    """The RMCEndpoint likelihood curves (one-sided: theta <= 0), its 2D contours with the power in the control region,
    and the RMC shape with the endpoint lowered by 0.5 MeV relative to the nominal."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    ax = axes[0]
    styles = {"sr": ("#eb6834", "-"), "cr": ("#2a78d6", "-"), "sr_cr": ("#1F1F1F", "-"), "cr_power_fixed": ("#2a78d6", "--")}
    names = dict(labels, cr_power_fixed=labels["cr"] + ", power fixed")
    for key, sc in endpoint.items():
        if sc is not None:
            ax.plot(*sc, color=styles[key][0], linestyle=styles[key][1], linewidth=2, label=names[key])
    for level, text in ((1., "68% CL"), (2.71, "95% CL one-sided")):
        ax.axhline(level, color="#5C5C5C", linestyle=":", linewidth=1)
        ax.text(-4.9, level, f" {text}", va="bottom", fontsize=9, color="#5C5C5C")
    ax.set_xlabel("RMCEndpoint $\\theta$ (unconstrained, $\\theta \\leq 0$)", fontsize=12)
    ax.set_ylabel(r"$-2\Delta\ln L$", fontsize=12)
    ax.set_ylim(0., 10.)
    ax.set_xlim(-5., 0.)
    top = ax.secondary_xaxis("top", functions=(lambda t: ENDPOINT_SHIFT * t, lambda d: d / ENDPOINT_SHIFT))
    top.set_xlabel("Endpoint shift [MeV] (linear in $\\theta$)", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=8, loc="upper center")
    ax = axes[1]
    if scan2d is not None and len(scan2d) > 10:
        ax.tricontour(scan2d[:, 0], scan2d[:, 1], scan2d[:, 2], levels=[2.30, 5.99], colors=["#1F1F1F", "#2a78d6"],
                      linestyles=["-", "--"])
        ax.plot([0.], [0.], marker="+", color="#eb6834", markersize=12)
        ax.plot([], [], color="#1F1F1F", label="68% CL")
        ax.plot([], [], color="#2a78d6", linestyle="--", label="95% CL")
        ax.legend(frameon=False, fontsize=9)
    ax.set_xlabel("RMCPower $\\theta$", fontsize=12)
    ax.set_ylabel("RMCEndpoint $\\theta$", fontsize=12)
    ax.set_title("control region, rate floating", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax = axes[2]
    n = curves["nominal"]
    for key, color, label in (("nominal", "#1F1F1F", "kmax = 101.87 MeV"), ("down", "#eb6834", "kmax - 0.5 MeV")):
        c = curves[key]
        ax.stairs((c / c.sum()) / (n / n.sum()), edges, color=color, linewidth=2, label=label)
    ax.set_xlabel("p [MeV/c]", fontsize=12)
    ax.set_ylabel("RMC shape / nominal (each normalized)", fontsize=12)
    ax.set_xlim(edges[0], edges[0] + 0.75 * (edges[-1] - edges[0]))
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=9)
    fig.suptitle(f"{cr_label}: RMC spectrum endpoint", fontsize=13)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def postfit_yields(model_dir, fit_file, channel, edges, processes, n_samples=1000, seed=1):
    """Per-bin yields of the processes of a channel at the fit result, and sampled from its covariance (all floating
    parameters with their correlations): {process: nominal}, {process: samples (n_samples x bins)}."""
    import ROOT
    ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
    ROOT.RooRandom.randomGenerator().SetSeed(seed)
    f = ROOT.TFile.Open(os.path.join(model_dir, "model.root"))
    w = f.Get("w")
    g = ROOT.TFile.Open(fit_file)
    fr = g.Get("fit_mdf")
    centers, widths = 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)
    obs = None
    terms = {}
    for p in processes:
        pdf = None
        for kind in ("shapeBkg", "shapeSig"):
            pdf = w.pdf(f"{kind}_{channel}_{p}_morph") or w.pdf(f"{kind}_{p}_{channel}Pdf") or w.pdf(f"{kind}_{p}_{channel}")
            if pdf:
                break
        norm = w.function(f"n_exp_bin{channel}_proc_{p}") or w.function(f"n_exp_final_bin{channel}_proc_{p}")
        if not pdf or not norm:
            raise RuntimeError(f"{channel}/{p}: no shape or normalization in the workspace")
        if obs is None:
            obs = next(v for v in pdf.getObservables(w.set("observables")) if v.GetName() != "CMS_channel")
        terms[p] = (pdf, norm)
    nset = ROOT.RooArgSet(obs)

    def evaluate():
        out = {}
        for p, (pdf, norm) in terms.items():
            values = np.empty(len(centers))
            for i, c in enumerate(centers):
                obs.setVal(c)
                values[i] = pdf.getVal(nset)
            out[p] = values * widths * norm.getVal()
        return out

    def assign(pars):
        for par in pars:
            v = w.var(par.GetName())
            if v:
                v.setVal(par.getVal())

    assign(fr.floatParsFinal())
    assign(fr.constPars())
    nominal = evaluate()
    samples = {p: np.empty((n_samples, len(centers))) for p in processes}
    for k in range(n_samples):
        assign(fr.randomizePars())
        for p, y in evaluate().items():
            samples[p][k] = y
    assign(fr.floatParsFinal())
    f.Close()
    g.Close()
    return nominal, samples


def plot_cr_postfit(model_dir, fit_file, channel, cr_label, out):
    """The control region data with the fitted model and its 68% band, and the RMC model alone with its band. Returns
    the relative 68% uncertainty of the total and RMC yields."""
    import ROOT
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    info = json.load(open(os.path.join(model_dir, "model.json")))["channels"][channel]
    edges = np.array(info["observable"]["edges"])
    processes = list(info["processes"])
    rmc = [p for p in processes if p.startswith("rmc")]
    nominal, samples = postfit_yields(model_dir, fit_file, channel, edges, processes)
    f = ROOT.TFile.Open(os.path.join(model_dir, "shapes", f"{channel}.root"))
    dh = f.Get("w").data(f"{channel}_data_obs")
    data = np.array([(dh.get(i), dh.weight())[1] for i in range(dh.numEntries())])
    f.Close()
    total_s = sum(samples[p] for p in processes)
    rmc_s = sum(samples[p] for p in rmc)
    total, rmc_n = sum(nominal[p] for p in processes), sum(nominal[p] for p in rmc)
    band = lambda s: np.percentile(s, [16., 84.], axis=0)
    t_lo, t_hi = band(total_s)
    r_lo, r_hi = band(rmc_s)
    centers = 0.5 * (edges[1:] + edges[:-1])
    asimov = info["data"]["source"] == "asimov"
    data_label = "Asimov data" if asimov else "Data"
    colors = {"rmc_int": "#e87ba4", "rmc_ext": "#eb6834", "cosmic": "#2a78d6", "other": "#1baf7a"}

    def step(ax, y, **kw):
        ax.stairs(y, edges, **kw)

    def band_fill(ax, lo, hi, **kw):
        ax.fill_between(edges, np.append(lo, lo[-1]), np.append(hi, hi[-1]), step="post", **kw)

    # 1) data with the model
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    others = [p for p in processes if p not in rmc and p != "cosmic"]
    stack = [("other", sum((nominal[p] for p in others), np.zeros(len(centers))), "Other (RPC, pbar)"),
             ("cosmic", nominal.get("cosmic", np.zeros(len(centers))), "Cosmic rays")] + \
            [(p, nominal[p], {"rmc_ext": "RMC (external)", "rmc_int": "RMC (internal)"}.get(p, p)) for p in rmc]
    cumulative = np.zeros(len(centers))
    for key, y, label in stack:
        ax.stairs(cumulative + y, edges, baseline=cumulative, fill=True, color=colors.get(key, "#94a4a2"), alpha=0.8,
                  label=label)
        cumulative = cumulative + y
    band_fill(ax, t_lo, t_hi, facecolor="none", edgecolor="#1F1F1F", hatch="////", linewidth=0,
              label="Fit uncertainty (68%)")
    step(ax, total, color="#1F1F1F", linewidth=1.5)
    ax.errorbar(centers, data, yerr=np.sqrt(np.clip(data, 0., None)), fmt="o", color="black", markersize=3.5,
                label=f"{data_label} ({data.sum():.0f})")
    ax.set_ylabel("Events / 0.2 MeV/c", fontsize=12)
    ax.set_title(f"{cr_label}\nfit with the RMC rate and spectrum power free", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    ax.set_ylim(0., 1.25 * max(data.max(), t_hi.max()))
    ok = total > 0
    band_fill(rax, np.where(ok, t_lo / np.where(ok, total, 1.), 1.), np.where(ok, t_hi / np.where(ok, total, 1.), 1.),
              color="#94a4a2", alpha=0.5)
    rax.errorbar(centers[ok], data[ok] / total[ok], yerr=np.sqrt(np.clip(data[ok], 0., None)) / total[ok], fmt="o",
                 color="black", markersize=3)
    rax.axhline(1., color="#1F1F1F", linewidth=1)
    rax.set_ylim(0.5, 1.5)
    rax.set_ylabel("Data / model")
    rax.set_xlabel("p [MeV/c]", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "rmc_cr_data_model.png"))
    plt.close(fig)

    # 2) the RMC model alone
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    band_fill(ax, r_lo, r_hi, color="#94a4a2", alpha=0.5, label="RMC fit uncertainty (68%)")
    step(ax, rmc_n, color="#1F1F1F", linewidth=2, label=f"RMC total ({rmc_n.sum():.1f})")
    for p in rmc:
        lo, hi = band(samples[p])
        band_fill(ax, lo, hi, color=colors.get(p), alpha=0.25)
        step(ax, nominal[p], color=colors.get(p), linewidth=1.5, linestyle="--",
             label=f"{ {'rmc_ext': 'RMC (external)', 'rmc_int': 'RMC (internal)'}.get(p, p)} ({nominal[p].sum():.1f})")
    ax.set_ylabel("Events / 0.2 MeV/c", fontsize=12)
    ax.set_title(f"{cr_label}\nRMC model after the control region fit", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    ok = rmc_n > 0
    rel_lo = np.where(ok, r_lo / np.where(ok, rmc_n, 1.) - 1., 0.)
    rel_hi = np.where(ok, r_hi / np.where(ok, rmc_n, 1.) - 1., 0.)
    band_fill(rax, rel_lo, rel_hi, color="#94a4a2", alpha=0.5)
    rax.axhline(0., color="#1F1F1F", linewidth=1)
    rax.set_ylabel("Rel. uncertainty")
    rax.set_xlabel("p [MeV/c]", fontsize=12)
    rax.set_ylim(-0.6, 1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "rmc_cr_rmc_model.png"))
    plt.close(fig)
    rel = lambda s, n: (np.percentile(s.sum(axis=1), [16., 84.]) / n.sum() - 1.).tolist()
    return {"total": rel(total_s, total), "rmc": rel(rmc_s, rmc_n),
            "rmc_bins_rel": [rel_lo.tolist(), rel_hi.tolist()], "edges": edges.tolist()}


def read_fit(path):
    """{parameter: (value, error)} of the fit result saved by MultiDimFit --saveFitResult."""
    import ROOT
    if not os.path.exists(path):
        return None
    f = ROOT.TFile.Open(path)
    fr = f.Get("fit_mdf")
    out = {p.GetName(): (p.getVal(), p.getError()) for p in fr.floatParsFinal()} if fr else None
    f.Close()
    return out


def plot_scans(scans, path, labels):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"sr": "#eb6834", "cr": "#2a78d6", "sr_cr": "#1F1F1F"}
    fig, ax = plt.subplots(figsize=(8, 6))
    for name, (x, y) in scans.items():
        ax.plot(x, y, color=colors.get(name), linewidth=2, label=labels[name])
    for level, text in ((1., "68% CL"), (3.84, "95% CL")):
        ax.axhline(level, color="#5C5C5C", linestyle=":", linewidth=1)
        ax.text(min(x.min() for x, _ in scans.values()), level, f" {text}",
                va="bottom", fontsize=9, color="#5C5C5C")
    ax.set_xlabel("RMCPower $\\theta$ (unconstrained)", fontsize=13)
    ax.set_ylabel(r"$-2\Delta\ln L$", fontsize=13)
    ax.set_ylim(0., 10.)
    top = ax.secondary_xaxis("top", functions=(lambda t: 2. + POWER_SHIFT * t, lambda n: (n - 2.) / POWER_SHIFT))
    top.set_xlabel("RMC spectrum power n (linear in $\\theta$, x(1-x)$^n$)", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)
    ax.set_title("RMC spectrum power: Asimov likelihood scans", fontsize=12)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_templates(model_dir, path, variation="RMCPower", labels=("n = 2", f"n = {2 + POWER_SHIFT:g}",
                                                                   f"n = {2 - POWER_SHIFT:g}")):
    """Nominal and up/down templates of a variation of the 0n RMC processes in each channel of a model."""
    import ROOT
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    info = json.load(open(os.path.join(model_dir, "model.json")))["channels"]
    panels = []
    for ch, c in info.items():
        f = ROOT.TFile.Open(os.path.join(model_dir, "shapes", f"{ch}.root"))
        w = f.Get("w")
        edges = np.array(c["observable"]["edges"])
        for p in c["processes"]:
            if not w.data(f"{ch}_{p}_{variation}Up"):
                continue
            curves = {}
            for key, suffix in (("nominal", ""), ("up", f"_{variation}Up"), ("down", f"_{variation}Down")):
                dh = w.data(f"{ch}_{p}{suffix}")
                curves[key] = np.array([(dh.get(i), dh.weight())[1] for i in range(dh.numEntries())])
            panels.append((ch, p, edges, curves))
        f.Close()
    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 4.5), squeeze=False)
    for ax, (ch, p, edges, curves) in zip(axes[0], panels):
        for key, color, label in (("nominal", "#1F1F1F", labels[0]), ("up", "#2a78d6", labels[1]),
                                  ("down", "#eb6834", labels[2])):
            if key == "up" and np.allclose(curves["up"], curves["nominal"]):
                continue
            ax.stairs(curves[key], edges, color=color, linewidth=2,
                      label=f"{label} ({curves[key].sum():.3g})")
        ax.set_title(f"{ch}: {p}", fontsize=11)
        ax.set_xlabel("p [MeV/c]")
        ax.set_ylabel("Events / bin")
        ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-s", "--selection", type=int, default=75, help="mu- --> e- selection set (default: %(default)s)")
    parser.add_argument("--cr-range", type=float, nargs=2, default=(94., 105.), metavar=("LO", "HI"),
                        help="momentum range of the e+ control region (default: 94 105)")
    parser.add_argument("-t", "--tag", default="", help="output tag")
    parser.add_argument("--out-dir", default=None, help="output directory (default: rmc_studies/mumem_<set>[_<tag>])")
    parser.add_argument("--model-args", default="", help="additional build_model.py arguments, quoted (e.g. --2bb)")
    parser.add_argument("--fit-args", default="", help="additional Combine arguments for every fit, quoted")
    parser.add_argument("--jobs", type=int, default=4, help="Combine fits in parallel, one core each (default: %(default)s)")
    parser.add_argument("--skip-models", action="store_true", help="use the existing models")
    parser.add_argument("--dry-run", action="store_true", help="print the commands without running them")
    args = parser.parse_args()
    if not 1 <= args.jobs <= 20:
        parser.error("--jobs must be between 1 and 20")

    name = f"mumem_{args.selection}" + (f"_{args.tag}" if args.tag else "")
    out = os.path.abspath(args.out_dir or os.path.join("rmc_studies", name))
    logs = os.path.join(out, "logs")
    os.makedirs(logs, exist_ok=True)
    sr, cr = args.selection, args.selection + 4000
    sr_ch, cr_ch = f"set_{sr}", f"cr4000_{sr}"
    build = (f"{shlex.quote(sys.executable)} {shlex.quote(os.path.join(HERE, 'build_model.py'))} -p mumem --no-cc --plots"
             f" -s {sr} -c 4000 --charge-inverted-range {args.cr_range[0]:g} {args.cr_range[1]:g}")
    models = {"free": "--rmc-endpoint free", "prior": "--rmc-power constrained --rmc-endpoint constrained"}
    configs = {"sr": f"mask_{cr_ch}=1", "cr": f"mask_{sr_ch}=1", "sr_cr": ""}
    labels = {"sr": f"signal region (set {sr}, 100-110 MeV/c) only",
              "cr": f"e+ control region (set {cr}, {args.cr_range[0]:g}-{args.cr_range[1]:g} MeV/c) only",
              "sr_cr": "signal + e+ control region"}

    os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    if not args.skip_models:
        for m, opts in models.items():
            status, _ = run(f"{build} {opts} {args.model_args} -o {shlex.quote(os.path.join(out, m))}", HERE,
                            os.path.join(logs, f"model_{m}.log"), args.dry_run)
            if status:
                sys.exit(f"rmc_constraint_study.py: building {m} failed, see {logs}/model_{m}.log")
            status, _ = run("text2workspace.py combined.txt -o model.root --channel-masks", os.path.join(out, m),
                            os.path.join(logs, f"t2w_{m}.log"), args.dry_run)
            if status:
                sys.exit(f"rmc_constraint_study.py: text2workspace of {m} failed, see {logs}/t2w_{m}.log")

    extra = f" {args.fit_args}" if args.fit_args else ""
    ranges = "RMCPower=-5,5:RMCEndpoint=-5,0:rmcNorm=0.2,5:r=-20,20"

    def fit(config, freeze=(), ranges=ranges):
        """Asimov (r = 0) options of a configuration: its channel mask, r frozen without the signal region."""
        sets = ",".join(x for x in ("r=0", configs[config]) if x)
        frozen = (["r"] if config == "cr" else []) + list(freeze)
        return (f"combine -M MultiDimFit model.root -t -1 --setParameters {sets}"
                + (f" --freezeParameters {','.join(frozen)}" if frozen else "")
                + f" --setParameterRanges {ranges} --cminDefaultMinimizerStrategy=1")

    # Only unconstrained parameters become POIs: in this Combine a redefined POI loses its constraint term, so intN
    # stays a nuisance and its post-fit uncertainty comes from the fit covariance.
    jobs = []
    for config in configs:
        jobs.append(("free", f"{config}_singles", f"{fit(config)} --redefineSignalPOIs RMCPower,RMCEndpoint,rmcNorm"
                                                  " --algo singles"
                                                  f" --saveFitResult{extra} -n .{config}_singles"))
        jobs.append(("free", f"{config}_scan_endpoint", f"{fit(config)} --redefineSignalPOIs RMCEndpoint"
                                                        " -P RMCEndpoint --floatOtherPOIs 1 --algo grid --points 51"
                                                        f"{extra} -n .{config}_scan_endpoint"))
        jobs.append(("free", f"{config}_scan", f"{fit(config, ranges='RMCPower=-4,4:RMCEndpoint=-5,0:rmcNorm=0.2,5:r=-20,20')}"
                                               " --redefineSignalPOIs RMCPower -P RMCPower --floatOtherPOIs 1"
                                               f" --algo grid --points 81{extra} -n .{config}_scan"))
    norm_ranges = "RMCPower=-4,4:RMCEndpoint=-5,0:rmcNorm=0.5,1.5:r=-20,20"
    jobs += [("free", "cr_scan_endpoint_power_fixed", f"{fit('cr', ['RMCPower'])} --redefineSignalPOIs RMCEndpoint"
                                                      " -P RMCEndpoint --floatOtherPOIs 1 --algo grid --points 51"
                                                      f"{extra} -n .cr_scan_endpoint_power_fixed"),
             ("free", "cr_scan_2d_spectrum", f"{fit('cr', ranges='RMCPower=-2.5,2.5:RMCEndpoint=-5,0:rmcNorm=0.2,5:r=-20,20')}"
                                             " --redefineSignalPOIs RMCPower,RMCEndpoint -P RMCPower -P RMCEndpoint"
                                             f" --floatOtherPOIs 1 --algo grid --points 1600{extra} -n .cr_scan_2d_spectrum"),
             ("free", "cr_scan_power_rate_fixed", f"{fit('cr', ['rmcNorm'], 'RMCPower=-4,4:RMCEndpoint=-5,0:r=-20,20')}"
                                                  " --redefineSignalPOIs RMCPower -P RMCPower --floatOtherPOIs 1"
                                                  f" --algo grid --points 81{extra} -n .cr_scan_power_rate_fixed"),
             ("free", "cr_scan_norm", f"{fit('cr', ranges=norm_ranges)} --redefineSignalPOIs rmcNorm -P rmcNorm"
                                      f" --floatOtherPOIs 1 --algo grid --points 81{extra} -n .cr_scan_norm"),
             ("free", "cr_scan_norm_power_fixed", f"{fit('cr', ['RMCPower'], norm_ranges)} --redefineSignalPOIs rmcNorm"
                                                  " -P rmcNorm --floatOtherPOIs 1 --algo grid --points 81"
                                                  f"{extra} -n .cr_scan_norm_power_fixed"),
             ("free", "cr_scan_2d", f"{fit('cr', ranges='RMCPower=-2.5,2.5:RMCEndpoint=-5,0:rmcNorm=0.7,1.3:r=-20,20')}"
                                    " --redefineSignalPOIs RMCPower,rmcNorm -P RMCPower -P rmcNorm --floatOtherPOIs 1"
                                    f" --algo grid --points 1600{extra} -n .cr_scan_2d"),
             ("prior", "sr_cr_singles", f"{fit('sr_cr')} --redefineSignalPOIs rmcNorm --algo singles --saveFitResult"
                                        f"{extra} -n .sr_cr_singles")]
    limit = ("combine -M AsymptoticLimits model.root -t -1 --rMin 0 --rMax 100 --cl 0.9 --cminDefaultMinimizerStrategy=1"
             " --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rAbsAcc 0.001 --rRelAcc 0.001 --noFitAsimov"
             f" --setParameterRanges {ranges}" + extra)
    mask_cr = f" --setParameters {configs['sr']}"
    jobs += [("free", "sr_limit", f"{limit}{mask_cr} -n .sr_limit"),
             ("free", "sr_cr_limit", f"{limit} -n .sr_cr_limit"),
             ("prior", "sr_limit", f"{limit}{mask_cr} -n .sr_limit"),
             ("prior", "sr_cr_limit", f"{limit} -n .sr_cr_limit"),
             ("prior", "sr_cr_limit_frozen", f"{limit} --freezeParameters RMCPower,RMCEndpoint -n .sr_cr_limit_frozen")]

    def execute(job):
        m, kind, command = job
        return job, run(command, os.path.join(out, m), os.path.join(logs, f"{m}_{kind}.log"), args.dry_run)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        done = list(pool.map(execute, jobs))
    if args.dry_run:
        return

    results = {m: {} for m in models}
    for (m, kind, _), (status, text) in done:
        if status:
            print(f"WARNING: {m} {kind} failed (status {status}), see {logs}/{m}_{kind}.log")
        if kind.endswith("singles"):
            results[m][kind] = parse_singles(text)
            results[m][kind.replace("singles", "fit")] = read_fit(os.path.join(out, m, f"multidimfit.{kind}.root"))
        elif "limit" in kind:
            results[m][kind] = parse_limits(text)
    scan = lambda kind, param="RMCPower": read_scan(
        os.path.join(out, "free", f"higgsCombine.{kind}.MultiDimFit.mH120.root"), param)
    scans = {c: s for c in configs if (s := scan(f"{c}_scan")) is not None}
    if scans:
        plot_scans(scans, os.path.join(out, "rmc_power_scan.png"), labels)
    plot_templates(os.path.join(out, "free"), os.path.join(out, "rmc_power_templates.png"))
    power = {"floating": scans.get("cr"), "fixed": scan("cr_scan_power_rate_fixed")}
    norm = {"floating": scan("cr_scan_norm", "rmcNorm"), "fixed": scan("cr_scan_norm_power_fixed", "rmcNorm")}
    scan2d = read_scan_2d(os.path.join(out, "free", "higgsCombine.cr_scan_2d.MultiDimFit.mH120.root"), "RMCPower",
                          "rmcNorm")
    edges, curves, decomposition = rmc_variation_decomposition(os.path.join(out, "free"), cr_ch)
    plot_cr_floating(power, norm, scan2d, edges, curves["total"], os.path.join(out, "rmc_cr_fit.png"), labels["cr"])
    results["free"]["cr_decomposition"] = decomposition
    endpoint = {c: scan(f"{c}_scan_endpoint", "RMCEndpoint") for c in configs}
    endpoint["cr_power_fixed"] = scan("cr_scan_endpoint_power_fixed", "RMCEndpoint")
    spectrum2d = read_scan_2d(os.path.join(out, "free", "higgsCombine.cr_scan_2d_spectrum.MultiDimFit.mH120.root"),
                              "RMCPower", "RMCEndpoint")
    e_edges, e_curves, e_decomposition = rmc_variation_decomposition(os.path.join(out, "free"), cr_ch, "RMCEndpoint")
    plot_cr_endpoint(endpoint, spectrum2d, e_edges, e_curves["total"], os.path.join(out, "rmc_cr_endpoint.png"),
                     labels["cr"], labels)
    plot_templates(os.path.join(out, "free"), os.path.join(out, "rmc_endpoint_templates.png"), "RMCEndpoint",
                   ("kmax = 101.87 MeV", "kmax (up = nominal)", "kmax - 0.5 MeV"))
    results["free"]["cr_endpoint_decomposition"] = e_decomposition
    bands = plot_cr_postfit(os.path.join(out, "free"), os.path.join(out, "free", "multidimfit.cr_singles.root"), cr_ch,
                            labels["cr"], out)
    results["free"]["cr_postfit_bands"] = bands

    fmt = lambda x, w=8: f"{x:{w}.3f}" if x is not None else f"{'-':>{w}s}"
    lines = [f"RMC constraint study: mu- --> e- set {sr}, e+ control region set {cr} at "
             f"{args.cr_range[0]:g}-{args.cr_range[1]:g} MeV/c (background-only Asimov)",
             f"build_model.py options: {args.model_args or '(defaults)'}",
             "RMC rate: rmcNorm (flat); spectrum power: RMCPower (shape only; theta = +-1 is n = 2 +- 0.5); spectrum "
             "endpoint: RMCEndpoint (shape only; theta = -1 lowers kmax by 0.5 MeV, theta <= 0)", "",
             "Profile-likelihood 68% CL intervals (both parameters free):",
             f"{'fit':8s} {'RMCPower':>9s} {'-68%':>8s} {'+68%':>8s} {'sigma(n)':>9s}   "
             f"{'rmcNorm':>8s} {'-68%':>8s} {'+68%':>8s}"]
    for config in configs:
        v = results["free"].get(f"{config}_singles", {})
        p, nrm = v.get("RMCPower"), v.get("rmcNorm")
        row = f"{config:8s} "
        row += (f"{p[0]:9.3f} {p[1]:8.3f} {p[2]:8.3f} {POWER_SHIFT * 0.5 * (abs(p[1]) + abs(p[2])):9.3f}   " if p
                else f"{'-':>9s} {'':8s} {'':8s} {'':9s}   ")
        row += f"{nrm[0]:8.3f} {nrm[1]:8.3f} {nrm[2]:8.3f}" if nrm else f"{'-':>8s}"
        lines.append(row)
    lines += ["", "RMCEndpoint, one-sided (theta <= 0): lower limits from the likelihood scans (others profiled)",
              f"{'fit':26s} {'68% CL theta':>13s} {'shift [MeV]':>12s} {'95% CL theta':>13s} {'shift [MeV]':>12s}"]
    for key, label in (("sr", "sr"), ("cr", "cr"), ("sr_cr", "sr_cr"), ("cr_power_fixed", "cr, power fixed")):
        sc = endpoint.get(key)
        lims = [interval(sc, level)[1] if sc is not None else None for level in (1., 2.71)]
        row = f"{label:26s} "
        for lim in lims:
            row += (f"{lim:13.2f} {ENDPOINT_SHIFT * lim:12.2f} " if lim is not None
                    else f"{'< -5':>13s} {'< -2.5':>12s} ")
        lines.append(row.rstrip())
    lines += ["", "Fit covariance uncertainties (RMCPower, RMCEndpoint, rmcNorm absolute; intN in units of its prior "
                  "sigma; RMCEndpoint sits at its boundary theta = 0, so its Hesse error is approximate):",
              f"{'fit':26s} {'RMCPower':>9s} {'RMCEndp.':>9s} {'rmcNorm':>8s} {'intN':>8s}"]
    for m, kind, label in [("free", f"{c}_fit", f"{c}") for c in configs] + \
                          [("prior", "sr_cr_fit", "sr_cr (unit prior on power)")]:
        f = results[m].get(kind) or {}
        lines.append(f"{label:26s} {fmt(f.get('RMCPower', (None, None))[1], 9)} "
                     f"{fmt(f.get('RMCEndpoint', (None, None))[1], 9)} "
                     f"{fmt(f.get('rmcNorm', (None, None))[1])} {fmt(f.get('intN', (None, None))[1])}")
    lines += ["", "Control region alone, likelihood scans (68% CL):", f"{'scan':34s} {'best fit':>9s} {'-68%':>8s} {'+68%':>8s}"]
    for label, sc in (("RMCPower, rate floating", power["floating"]), ("RMCPower, rate fixed (rmcNorm = 1)", power["fixed"]),
                      ("rmcNorm, power floating", norm["floating"]), ("rmcNorm, power fixed (n = 2)", norm["fixed"])):
        iv = interval(sc)
        if iv is None:
            lines.append(f"{label:34s} {'-':>9s}")
            continue
        best, lo, hi = iv
        lines.append(f"{label:34s} {best:9.3f} {fmt(lo - best if lo is not None else None)} "
                     f"{fmt(hi - best if hi is not None else None)}")
    lines += ["", f"RMCPower +-1 in the control region: yield ratio (1 = shape only) and the range of the shape ratio",
              f"{'process':10s} {'yield':>9s} {'up yield':>9s} {'down yield':>10s} {'up shape':>15s} {'down shape':>15s}"]
    for p, e in decomposition.items():
        lines.append(f"{p:10s} {e['yield']:9.2f} {e['up_yield_ratio']:9.3f} {e['down_yield_ratio']:10.3f} "
                     f"{e['up_shape_ratio_range'][0]:6.3f}-{e['up_shape_ratio_range'][1]:6.3f} "
                     f"{e['down_shape_ratio_range'][0]:6.3f}-{e['down_shape_ratio_range'][1]:6.3f}")
    lines += ["", "Control region fit (rmc_cr_data_model.png, rmc_cr_rmc_model.png): 68% uncertainty of the yields in "
                  f"{args.cr_range[0]:g}-{args.cr_range[1]:g} MeV/c, from sampling the fit covariance",
              f"  total model: {100 * bands['total'][0]:+.1f}% / {100 * bands['total'][1]:+.1f}%",
              f"  RMC:         {100 * bands['rmc'][0]:+.1f}% / {100 * bands['rmc'][1]:+.1f}%"]
    lines += ["", "RMCEndpoint -1 (kmax - 0.5 MeV) in the control region: the range of the shape ratio",
              f"{'process':10s} {'yield':>9s} {'down shape':>15s}"]
    for p, e in e_decomposition.items():
        lines.append(f"{p:10s} {e['yield']:9.2f} {e['down_shape_ratio_range'][0]:6.3f}-{e['down_shape_ratio_range'][1]:6.3f}")
    lines += ["", "Expected 90% CL upper limits on r (Asimov median, -1 sigma, +1 sigma):"]
    for m, kind, label in (("free", "sr_limit", "signal region only, spectrum free"),
                           ("prior", "sr_limit", "signal region only, unit priors"),
                           ("free", "sr_cr_limit", "signal + e+ control region, spectrum free"),
                           ("prior", "sr_cr_limit", "signal + e+ control region, unit priors"),
                           ("prior", "sr_cr_limit_frozen", "signal + e+ control region, spectrum frozen")):
        lim = results[m].get(kind, {})
        lines.append(f"  {label:46s} {fmt(lim.get('median'))} {fmt(lim.get('16.0'))} {fmt(lim.get('84.0'))}")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(out, "summary.txt"), "w") as f:
        f.write(text + "\n")
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump({"results": results, "labels": labels}, f, indent=2, default=str)
    print(f"Outputs in {out}")


if __name__ == "__main__":
    main()
