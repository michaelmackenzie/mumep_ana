#! /usr/bin/env python3
"""RMC normalization and shape uncertainties for the mu- --> e- signal region (e- tracks, p > 100 MeV/c), motivated by
an e+ control region measurement that is not fit simultaneously: for analyses without the e+ channel.

    ./rmc_sr_prescription.py [-s 75] [--sr-range 100 110] [--cr-range 94 105] [--hist-path DIR] [--endpoint-prior 0.5]

Rerun it whenever the selection of the sets (e.g. a cut in set 75), the momentum ranges, or the inputs change: after
the histograms are remade, the control region study is rerun automatically (its outputs are older than the inputs).
See README.md (RMC prescription for the signal region) for what it needs and how to read it.

Steps:
 1. Checks: both 0n RMC samples have the photon spectrum variations (sys_<set> histograms, RMCPower and RMCEndpoint)
    for the signal-region set S and the e+ control-region set S + 4000, made after the Mu2eEvtAna RMCEndpoint fix,
    and the signal region has MC events.
 2. The e+ control region fit (rmc_constraint_study.py, run if its outputs are missing or with --run-study): the RMC
    rate and the likelihood of (RMCPower, RMCEndpoint) with the rate profiled.
 3. Normalization: the signal-region external RMC yield predicted from the control region is
        N_SR = N_CR(measured) x f_ext x [SR / CR](spectrum)
    with the rate measured in the control region (rmcNorm), the external fraction (intN, 4.5% on the internal/external
    ratio), and the signal/control region yield ratio. The ratio depends on the photon spectrum: e- above 100 MeV/c need
    photons within ~2 MeV of the endpoint, while the control region sees k ~ 94-98 MeV. The spectrum is sampled from the
    control region likelihood, with a flat power prior and a flat prior on an endpoint lowered by up to
    --endpoint-prior MeV (the control region alone cannot separate the two), and the variation factors of the yields
    come from the MC (their arbitrary Mu2eEvtAna normalization cancels in the ratio). The recommended lnN is the
    envelope of the 16% and 84% quantiles of the endpoint-fixed and endpoint-prior cases around the nominal yield.
 4. Shape: a generator-level model of the signal-region e- spectrum (the photon spectrum x(1-x)^n up to kmax, a
    pair-conversion lepton taking a uniform fraction of the photon energy, smeared with the DIO e- momentum response
    shifted by an extra momentum loss fit to the MC shape) is sampled the same way. The shape uncertainty is expressed
    with the RMCPower variation (shape only, n = 2 +- 0.5 at theta = +-1): the prior width on theta that reproduces the
    sampled shape spread (event-weighted), checked over the 68% CL range of the momentum loss. The fraction of the
    spread it captures is compared with the leading principal component, which is also written out as ratio curves.

Outputs in <study-dir>/sr_prescription/: summary.txt (recommendation and checks), prescription.json,
rmc_sr_shape.root (ratio histograms of the leading shape mode), rmc_sr_prescription.png.
"""

import argparse
import json
import os
import shlex
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
KMAX = 101.866
POWER, POWER_SHIFT, ENDPOINT_SHIFT = 2., 0.5, 0.5
INT_RATIO_UNC = 0.045   # intN
MIN_SR_MC = 5           # fewer MC events in the signal region: the shape fit and the MC factors are unreliable


def import_mc():
    import build_model  # noqa: F401  (locates mu2ecombine)
    return build_model.mc


def integral(h, lo, hi):
    return h.Integral(h.FindBin(lo + 1e-4), h.FindBin(hi - 1e-4))


def mc_factors(samples, process, selection, lo, hi):
    """Yield (weighted) and the raw variation factors (power up/down, endpoint down) of an RMC sample in [lo, hi]."""
    nominal = integral(samples[process].hist(f"trk_{selection}/obs"), lo, hi)
    up, down = samples[process].variation("RMCPower", selection)
    _, e_down = samples[process].variation("RMCEndpoint", selection)
    if nominal <= 0:
        return 0., (np.nan, np.nan, np.nan)
    return nominal, (integral(up, lo, hi) / nominal, integral(down, lo, hi) / nominal, integral(e_down, lo, hi) / nominal)


def spectrum(k, n, kmax):
    x = k / kmax
    return np.where((x > 0) & (x < 1), x * (1. - x) ** n, 0.)


def tail(kmin, n, kmax):
    k = np.linspace(kmin, kmax, 20001)
    return np.trapz(spectrum(k, n, kmax), k)


def model_shape(mc, n, kmax, response, edges, step=0.01):
    """e- momentum distribution of the generator-level model in the channel binning: the photon spectrum normalized
    above k = 97 MeV (as the Mu2eEvtAna variations), a lepton momentum uniform in [0, k], smeared with the response."""
    p = np.arange(80. + step / 2, KMAX + 0.5, step)
    k = np.linspace(80., KMAX, 4001)
    s = spectrum(k, n, kmax) / tail(97., n, kmax)
    # density of p: integral over k > p of s(k) / k
    dens = np.array([np.trapz(np.where(k > x, s / k, 0.), k) for x in p])
    truth = mc.hist.make(np.append(p - step / 2, p[-1] + step / 2), dens * step, None, "truth")
    full = np.round(np.arange(80., edges[-1] + 1e-6, edges[1] - edges[0]), 6)
    out = mc.convolved(truth, response, full)
    _, values, _ = mc.hist.arrays(out)
    i0 = int(np.argmin(np.abs(full - edges[0])))
    return values[i0:i0 + len(edges) - 1]


def shifted_response(mc, response, delta):
    edges, values, variances = mc.hist.arrays(response)
    return mc.hist.make(edges - delta, values, variances, "response")


def fit_delta(mc, response, edges, counts, grid=np.arange(0., 1.51, 0.05)):
    """The extra momentum loss of the model fit to the MC counts (Poisson likelihood of the normalized shape): the best
    value and the 68% CL range (Delta -ln L = 0.5)."""
    nll = []
    for d in grid:
        m = model_shape(mc, POWER, KMAX, shifted_response(mc, response, d), edges)
        mu = m / m.sum() * counts.sum()
        nll.append(np.sum(mu - counts * np.log(np.clip(mu, 1e-12, None))))
    nll = np.array(nll) - min(nll)
    inside = grid[nll < 0.5]
    return float(grid[np.argmin(nll)]), (float(inside.min()), float(inside.max()))


def read_scan_2d(path):
    import ROOT
    f = ROOT.TFile.Open(path)
    t = f.Get("limit")
    pts = np.array([(e.RMCPower, e.RMCEndpoint, 2. * e.deltaNLL) for i, e in enumerate(t) if i > 0])
    f.Close()
    pts[:, 2] -= pts[:, 2].min()
    return pts


def interp(theta, up, down):
    """Exponential (log-normal-like) interpolation of a factor: up^theta for theta > 0, down^-theta below."""
    return np.where(theta >= 0, np.power(up, theta), np.power(down, -theta))


def check_inputs(build_model, samples, sr, cr):
    """The 0n RMC samples have the spectrum variations of the sets, made after the RMCEndpoint fix."""
    for process, sets in (("rmc_ext", (sr, cr)), ("rmc_int", (cr,))):
        for s in sets:
            try:
                build_model.rmc_spectrum_variation(samples[process], "RMCEndpoint", s)
                samples[process].variation("RMCPower", s)
            except Exception as err:  # noqa: BLE001  (an InputError: the histograms are missing)
                sys.exit(f"rmc_sr_prescription.py: {process}, set {s}: no RMC spectrum variations ({err}). The sets "
                         "need sys_<set> histograms from Mu2eEvtAna (see README.md, RMC prescription)")


def sample_spectrum(scan, prior, n, rng):
    """(theta_power, theta_endpoint) sampled from the control region likelihood scan, with a flat prior on the endpoint
    lowered by up to `prior` MeV (0: fixed), and the posterior weight at the edge of the scanned power range."""
    theta_min = -prior / ENDPOINT_SHIFT
    sel = scan[:, 1] >= theta_min - 1e-9 if prior > 0 else np.isclose(scan[:, 1], scan[:, 1].max())
    pts = scan[sel]
    w = np.exp(-0.5 * pts[:, 2])
    w /= w.sum()
    step_p = np.min(np.diff(np.unique(scan[:, 0])))
    step_e = np.min(np.diff(np.unique(scan[:, 1])))
    idx = rng.choice(len(pts), size=n, p=w)
    tp = pts[idx, 0] + rng.uniform(-step_p / 2, step_p / 2, n)
    te = (np.clip(pts[idx, 1] + rng.uniform(-step_e / 2, step_e / 2, n), theta_min, 0.) if prior > 0 else np.zeros(n))
    edge = w[np.abs(pts[:, 0]) > np.abs(scan[:, 0]).max() - step_p].sum()
    return tp, te, float(edge)


def shape_analysis(model, tp, te):
    """The sampled signal-region shapes (per-bin exponential interpolation of the model templates), event-weighted:
    the rms shape change, the RMCPower prior width (theta) that reproduces the spread and the fraction of the spread it
    captures, the same for the endpoint, and the leading principal component (ratio curve, fraction)."""
    nom = model["nominal"]
    f = nom / nom.sum()
    ok = f > 1e-4
    rat = {k: np.where(ok, model[k] / np.where(ok, nom, 1.), 1.) for k in ("power_up", "power_down", "endpoint_down")}
    t = nom * interp(tp[:, None], rat["power_up"], rat["power_down"]) * np.power(rat["endpoint_down"], -te[:, None])
    ratio = (t / t.sum(axis=1, keepdims=True)) / f
    dev = (ratio - ratio.mean(axis=0))[:, ok] * np.sqrt(f[ok])
    total = np.sum(dev ** 2)
    shape_only = lambda k: ((model[k] / model[k].sum()) / f)[ok] - 1.
    directions = {"power": 0.5 * (shape_only("power_up") - shape_only("power_down")),
                  "endpoint": shape_only("endpoint_down")}
    out = {"rms": float(np.sqrt(total / len(dev)))}
    for name, v in directions.items():
        u = v * np.sqrt(f[ok])
        proj = dev @ (u / np.linalg.norm(u))
        out[name] = {"theta_sigma": float(np.std(proj) / np.linalg.norm(u)),
                     "captured": float(np.sum(proj ** 2) / total)}
    _, sv, vt = np.linalg.svd(dev, full_matrices=False)
    pc = vt[0] * sv[0] / np.sqrt(len(dev)) / np.sqrt(f[ok])
    if pc[0] < 0:
        pc = -pc
    pc_full = np.zeros(len(nom))
    pc_full[ok] = pc
    out["leading"] = {"captured": float(sv[0] ** 2 / np.sum(sv ** 2)), "ratio_up": (1. + pc_full).tolist(),
                      "ratio_down": (1. - pc_full).tolist(),
                      "band": np.percentile(ratio, [16., 84.], axis=0).tolist()}
    return out


def input_files(samples):
    """{path: modification time} of the histogram files the model reads."""
    out = {}
    for sample in samples.values():
        path = sample.file.path
        if os.path.exists(path):
            out[path] = os.path.getmtime(path)
    return out


def run_study(args, study_dir):
    command = [sys.executable, os.path.join(HERE, "rmc_constraint_study.py"), "-s", str(args.selection),
               "--cr-range", f"{args.cr_range[0]:g}", f"{args.cr_range[1]:g}", "--out-dir", study_dir,
               "--jobs", str(args.jobs)]
    if args.hist_path:
        command += [f"--model-args=--hist-path {args.hist_path}"]
    print("Running the control region study:", " ".join(shlex.quote(c) for c in command))
    if subprocess.run(command, stdout=subprocess.DEVNULL).returncode != 0:
        sys.exit(f"rmc_sr_prescription.py: rmc_constraint_study.py failed, see {study_dir}/logs")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-s", "--selection", type=int, default=75, help="mu- --> e- selection set (default: %(default)s)")
    parser.add_argument("--sr-range", type=float, nargs=2, default=(100., 110.), metavar=("LO", "HI"),
                        help="signal region momentum range (default: 100 110)")
    parser.add_argument("--cr-range", type=float, nargs=2, default=(94., 105.), metavar=("LO", "HI"),
                        help="e+ control region momentum range (default: 94 105)")
    parser.add_argument("--hist-path", default=None,
                        help="histogram directory (default: CONVANA_HIST_PATH, or the analysis histogram directory)")
    parser.add_argument("--endpoint-prior", type=float, default=0.5,
                        help="largest lowering of the photon spectrum endpoint, in MeV, with a flat prior (default: "
                             "%(default)s, the Mu2eEvtAna example shift); the endpoint-fixed case is always evaluated")
    parser.add_argument("-t", "--tag", default="", help="output tag")
    parser.add_argument("--study-dir", default=None,
                        help="rmc_constraint_study.py output directory (default: rmc_studies/mumem_<set>[_<tag>])")
    parser.add_argument("--run-study", action="store_true", help="rerun the control region study even if it exists")
    parser.add_argument("--jobs", type=int, default=4, help="Combine jobs of the control region study (default: 4)")
    parser.add_argument("--samples", type=int, default=200000, help="posterior samples (default: %(default)s)")
    args = parser.parse_args()

    mc = import_mc()
    import ROOT
    import build_model
    from samples import Exposure, make_samples
    ROOT.gROOT.SetBatch(True)
    study_dir = args.study_dir or os.path.join("rmc_studies", f"mumem_{args.selection}" + (f"_{args.tag}" if args.tag else ""))
    out = os.path.join(study_dir, "sr_prescription")
    samples = make_samples(Exposure("run1a"), args.hist_path)
    sr, cr = args.selection, args.selection + 4000
    (slo, shi), (clo, chi) = args.sr_range, args.cr_range

    # 1. inputs
    check_inputs(build_model, samples, sr, cr)
    scan_path = os.path.join(study_dir, "free", "higgsCombine.cr_scan_2d_spectrum.MultiDimFit.mH120.root")
    summary_path = os.path.join(study_dir, "summary.json")
    files = input_files(samples)
    # the control region study is redone when it is missing or older than an input histogram file (a changed selection)
    newer = [f for f, t in files.items() if os.path.exists(summary_path) and t > os.path.getmtime(summary_path)]
    if newer:
        print(f"Input histograms are newer than the control region study ({', '.join(map(os.path.basename, newer))}): "
              "rerunning it")
    if args.run_study or newer or not os.path.exists(scan_path) or not os.path.exists(summary_path):
        run_study(args, study_dir)
    os.makedirs(out, exist_ok=True)
    checks = []

    # 2. MC variation factors: signal region external RMC, control region external + internal RMC
    n_sr, f_sr = mc_factors(samples, "rmc_ext", sr, slo, shi)
    n_ext, f_ext = mc_factors(samples, "rmc_ext", cr, clo, chi)
    n_int, f_int = mc_factors(samples, "rmc_int", cr, clo, chi)
    f_cr = tuple((n_ext * a + n_int * b) / (n_ext + n_int) for a, b in zip(f_ext, f_int))
    ratio = tuple(a / b for a, b in zip(f_sr, f_cr))   # SR / CR yield ratio: power up, down, endpoint down
    int_fraction = n_int / (n_ext + n_int)
    edges = np.round(np.arange(slo, shi + 1e-6, 0.2), 6)
    _, mc_values, mc_variances = mc.hist.arrays(mc.rebin_to(samples["rmc_ext"].hist(f"trk_{sr}/obs"), edges))
    _, mc_counts, _ = mc.hist.arrays(mc.rebin_to(samples["rmc_ext"].file.hist(f"trk_{sr}/obs"), edges))
    if mc_counts.sum() < MIN_SR_MC:
        sys.exit(f"rmc_sr_prescription.py: only {mc_counts.sum():.0f} MC events of rmc_ext in set {sr} at "
                 f"{slo:g}-{shi:g} MeV/c: too few for the factors and the shape fit")
    checks.append(f"signal region MC: {mc_counts.sum():.0f} external RMC events ({n_sr:.3f} expected); the MC "
                  f"statistical uncertainty of this yield (+-{100 * np.sqrt(mc_variances.sum()) / mc_values.sum():.0f}%) "
                  "is not part of this prescription")

    # 3. generator-level model of the signal region e- spectrum, with its momentum loss fit to the MC shape
    response = samples["dio"].hist("trk_10/dP")
    delta, delta_range = fit_delta(mc, response, edges, mc_counts)

    def model_set(d):
        resp = shifted_response(mc, response, d)
        return {"nominal": model_shape(mc, POWER, KMAX, resp, edges),
                "power_up": model_shape(mc, POWER + POWER_SHIFT, KMAX, resp, edges),
                "power_down": model_shape(mc, POWER - POWER_SHIFT, KMAX, resp, edges),
                "endpoint_down": model_shape(mc, POWER, KMAX - ENDPOINT_SHIFT, resp, edges)}
    models = {d: model_set(d) for d in dict.fromkeys((delta,) + delta_range)}
    f_model = tuple(models[delta][k].sum() / models[delta]["nominal"].sum()
                    for k in ("power_up", "power_down", "endpoint_down"))
    worst = max(abs(np.log(a / b)) for a, b in zip(f_model, f_sr))
    checks.append(f"model vs MC signal region variation factors: largest difference {100 * (np.exp(worst) - 1):.0f}% "
                  f"(MC from {mc_counts.sum():.0f} events)" + ("  <-- CHECK" if worst > np.log(1.4) else ""))
    if delta_range[0] <= 0. or delta_range[1] >= 1.5:
        checks.append(f"momentum loss {delta:.2f} MeV reaches the edge of its fit range ({delta_range[0]:.2f}-"
                      f"{delta_range[1]:.2f}): the model shape may not describe the MC  <-- CHECK")

    # 4. the control region measurement, the posterior samples, and the signal region yield and shape
    scan = read_scan_2d(scan_path)
    fit = json.load(open(os.path.join(study_dir, "summary.json")))["results"]["free"]["cr_fit"]
    sigma_norm, sigma_ext = fit["rmcNorm"][1], int_fraction * INT_RATIO_UNC
    rng = np.random.default_rng(1)
    q = lambda x: [float(v) for v in np.percentile(x, [16., 50., 84.])]
    scenarios = {}
    for prior in dict.fromkeys((0., args.endpoint_prior)):
        tp, te, edge = sample_spectrum(scan, prior, args.samples, rng)
        if edge > 0.01:
            checks.append(f"endpoint prior {prior:g}: {100 * edge:.1f}% of the posterior at the edge of the power scan "
                          "<-- CHECK (widen the scan in rmc_constraint_study.py)")
        spectrum_factor = interp(tp, ratio[0], ratio[1]) * np.power(ratio[2], -te)
        rate_factor = np.exp(rng.normal(0., sigma_norm, args.samples) + rng.normal(0., sigma_ext, args.samples))
        sub = slice(0, min(20000, args.samples))
        scenarios[prior] = {"rate": q(rate_factor), "spectrum": q(spectrum_factor),
                            "total": q(spectrum_factor * rate_factor),
                            "shape": {f"{d:.2f}": shape_analysis(m, tp[sub], te[sub]) for d, m in models.items()},
                            "samples": (tp, spectrum_factor * rate_factor)}

    # recommendation: lnN envelope around the nominal; RMCPower shape width (largest over the momentum loss range)
    kappa_down = min(s["total"][0] for s in scenarios.values())
    kappa_up = max(s["total"][2] for s in scenarios.values())
    main = scenarios[args.endpoint_prior]
    widths = {d: a["power"]["theta_sigma"] for d, a in main["shape"].items()}
    theta = max(widths.values())
    captured = min(a["power"]["captured"] for a in main["shape"].values())
    theta_fixed = max(a["power"]["theta_sigma"] for a in scenarios[0.]["shape"].values())

    # outputs
    leading = main["shape"][f"{delta:.2f}"]["leading"]
    nominal_frac = models[delta]["nominal"] / models[delta]["nominal"].sum()
    f = ROOT.TFile(os.path.join(out, "rmc_sr_shape.root"), "RECREATE")
    for name, values in (("rmc_shape_up", leading["ratio_up"]), ("rmc_shape_down", leading["ratio_down"]),
                         ("rmc_shape_nominal", nominal_frac)):
        h = ROOT.TH1D(name, f"{name};p [MeV/c];" + ("fraction" if name.endswith("nominal") else "ratio to nominal"),
                      len(edges) - 1, float(edges[0]), float(edges[-1]))
        for i, v in enumerate(values):
            h.SetBinContent(i + 1, v)
        h.Write()
    f.Close()

    fmt3 = lambda v: " / ".join(f"{x:.3f}" for x in v)
    lines = [f"RMC prescription for the mu- --> e- signal region: set {sr}, e-, {slo:g}-{shi:g} MeV/c, from the e+ "
             f"control region set {cr} at {clo:g}-{chi:g} MeV/c (Run 1A Asimov; not a simultaneous fit)",
             f"histograms: {args.hist_path or os.environ.get('CONVANA_HIST_PATH') or 'default directory'}; "
             f"control region study: {study_dir}", "",
             "RECOMMENDATION",
             f"  normalization: lnN {kappa_down:.2f}/{kappa_up:.2f} on the signal-region RMC (envelope of the endpoint "
             f"fixed and lowered by up to {args.endpoint_prior:g} MeV, around the nominal yield)",
             f"  shape: the Mu2eEvtAna RMCPower variation, shape only (normalized to the nominal yield in the channel), "
             f"with a Gaussian prior of width {theta:.2f} on theta, i.e. n = 2 +- {POWER_SHIFT * theta:.2f}",
             f"         ({theta_fixed:.2f}, n = 2 +- {POWER_SHIFT * theta_fixed:.2f}, if the endpoint is taken as known); "
             f"it captures >= {100 * captured:.0f}% of the sampled shape spread",
             "  smooth the RMCPower templates (e.g. a fit to the variations, as build_model.py does): the signal region "
             "MC is sparse", "",
             "Inputs",
             f"  control region RMC: {n_ext + n_int:.1f} events (internal fraction {int_fraction:.2f}); rate measured to "
             f"+-{100 * sigma_norm:.1f}%; external fraction +-{100 * sigma_ext:.1f}%",
             f"  variation factors (power n = 2.5 / 1.5, endpoint -0.5 MeV; raw Mu2eEvtAna normalization, it cancels in "
             "the ratio):",
             f"    signal region (MC):         {f_sr[0]:.2f} / {f_sr[1]:.2f}, {f_sr[2]:.2f}",
             f"    signal region (model):      {f_model[0]:.2f} / {f_model[1]:.2f}, {f_model[2]:.2f}   (momentum loss "
             f"{delta:.2f} MeV, 68% CL {delta_range[0]:.2f}-{delta_range[1]:.2f}, fit to the MC shape)",
             f"    control region (MC):        {f_cr[0]:.2f} / {f_cr[1]:.2f}, {f_cr[2]:.2f}",
             f"    signal / control region:    {ratio[0]:.2f} / {ratio[1]:.2f}, {ratio[2]:.2f}", "",
             "Signal-region RMC yield factor (16% / 50% / 84%)"]
    for prior, s in scenarios.items():
        lines.append(f"  {'endpoint fixed' if prior == 0 else f'endpoint lowered by up to {prior:g} MeV':38s} "
                     f"rate {fmt3(s['rate'])} | spectrum {fmt3(s['spectrum'])} | total {fmt3(s['total'])}")
    lines += ["", "Signal-region shape (event-weighted), per momentum loss: rms change; RMCPower width (theta) and "
                  "captured fraction; endpoint width and fraction; leading mode fraction"]
    for prior, s in scenarios.items():
        for d, a in s["shape"].items():
            lines.append(f"  {'endpoint fixed' if prior == 0 else f'endpoint <= {prior:g} MeV':20s} loss {d} MeV: rms "
                         f"{100 * a['rms']:4.1f}% | power {a['power']['theta_sigma']:.2f} ({100 * a['power']['captured']:.0f}%)"
                         f" | endpoint {a['endpoint']['theta_sigma']:.2f} ({100 * a['endpoint']['captured']:.0f}%)"
                         f" | leading {100 * a['leading']['captured']:.0f}%")
    lines += ["", "Checks"] + [f"  {c}" for c in checks]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(out, "summary.txt"), "w") as fo:
        fo.write(text + "\n")
    with open(os.path.join(out, "prescription.json"), "w") as fo:
        json.dump({"selection": sr, "control_region_selection": cr, "sr_range": [slo, shi], "cr_range": [clo, chi],
                   "endpoint_prior_MeV": args.endpoint_prior,
                   "recommendation": {"lnN": [kappa_down, kappa_up], "rmc_power_theta_sigma": theta,
                                      "rmc_power_n_sigma": POWER_SHIFT * theta,
                                      "rmc_power_theta_sigma_endpoint_fixed": theta_fixed},
                   "edges": edges.tolist(), "nominal_shape": nominal_frac.tolist(),
                   "scenarios": {str(p): {k: v for k, v in s.items() if k != "samples"} for p, s in scenarios.items()},
                   "inputs": {"sr_mc_factors": f_sr, "sr_model_factors": f_model, "cr_mc_factors": f_cr,
                              "sr_cr_ratio": ratio, "rate_sigma": sigma_norm, "external_fraction_sigma": sigma_ext,
                              "momentum_loss_MeV": delta, "momentum_loss_68cl": delta_range},
                   "checks": checks,
                   "input_files": {f: time.strftime("%Y-%m-%d %H:%M", time.localtime(t)) for f, t in files.items()}}, fo, indent=2)

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    ax = axes[0]
    bins = np.linspace(0., 3., 121)
    for prior, s in scenarios.items():
        ax.hist(s["samples"][1], bins=bins, histtype="step", linewidth=2, density=True,
                label="endpoint fixed" if prior == 0 else f"endpoint $\\leq$ {prior:g} MeV lower")
    for k in (kappa_down, kappa_up):
        ax.axvline(k, color="#5C5C5C", linestyle="--")
    ax.axvline(1., color="#1F1F1F", linestyle=":")
    ax.set_xlabel("signal-region RMC yield / nominal")
    ax.set_title(f"normalization: lnN {kappa_down:.2f}/{kappa_up:.2f}", fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax = axes[1]
    centers = 0.5 * (edges[1:] + edges[:-1])
    ax.stairs(nominal_frac, edges, color="#1F1F1F", linewidth=2, label=f"model (momentum loss {delta:.2f} MeV)")
    ax.errorbar(centers, mc_values / mc_values.sum(), yerr=np.sqrt(mc_variances) / mc_values.sum(), fmt="o",
                color="#eb6834", markersize=4, label=f"MC ({mc_counts.sum():.0f} events)")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1.)
    ax.set_xlim(slo, slo + 3.)
    ax.set_xlabel("p [MeV/c]")
    ax.set_ylabel("fraction per 0.2 MeV/c")
    ax.set_title("signal-region external RMC shape", fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax = axes[2]
    m = models[delta]
    nom = m["nominal"]
    ok = nom > 1e-4 * nom.max()
    band = np.array(leading["band"])
    ax.fill_between(edges, np.append(band[0], band[0][-1]), np.append(band[1], band[1][-1]), step="post",
                    color="#94a4a2", alpha=0.5, label="sampled shapes (68%)")
    for key, color, sign in (("power_up", "#2a78d6", 1.), ("power_down", "#eb6834", -1.)):
        r = np.where(ok, (m[key] / m[key].sum()) / (nom / nom.sum()), 1.)
        scaled = np.exp(theta * np.log(np.clip(r, 1e-9, None)))   # the RMCPower variation at theta = +-width
        ax.stairs(scaled / np.sum(scaled * nom / nom.sum()), edges, color=color, linewidth=2,
                  label=f"RMCPower shape, n = 2 {'+' if sign > 0 else '-'} {POWER_SHIFT * theta:.2f}")
    ax.axhline(1., color="#1F1F1F", linewidth=1)
    ax.set_xlim(slo, slo + 3.)
    ax.set_ylim(0., 2.5)
    ax.set_xlabel("p [MeV/c]")
    ax.set_ylabel("ratio to the nominal shape")
    ax.set_title("shape: sampled spread and the RMCPower width", fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "rmc_sr_prescription.png"))
    plt.close(fig)
    print(f"Outputs in {out}")


if __name__ == "__main__":
    main()
