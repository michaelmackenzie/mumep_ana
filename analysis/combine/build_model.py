#! /usr/bin/env python3
"""Build the mu- --> e- / mu- --> e+ Combine model of a list of selection sets with mu2ecombine.

    python3 build_model.py -p mumep -s "40" [-c "1000 4000"] [-t <tag>] [options]

Writes models/<process>_<sets>[_cr<offset>...][_<tag>]/ (or --output): the shape model (combined.txt, a card per
channel, shapes/, model.json, commands.sh, figures/) and the cut-and-count model in cc/. Run Combine on the cards
directly, e.g. with the commands in commands.sh. Every channel is written as a RooWorkspace on its observable
obs_<set>, so any model of the same channels (templates, functions, envelope) can generate the data another fits.

--only PROCESS builds and fits only that process (to develop or debug its shape), draws it to <output>/debug/, and
exits without writing a model.

--range LO HI restricts the signal-region channels to a momentum range, keeping the full-range shapes of the MC
processes bin by bin (an envelope is fit to the data in the range); see sensitivity_vs_fit_range.py.

Channels: set_<S> for each selection set S, and cr<O>_<S> for each control region offset O (set S + O). The offsets
are the region flags of Mu2eEvtAna's ConvAna: 1000 inverted CRV veto, 2000 inverted time cut, 4000 charge-inverted.
They combine, e.g. 3000. A charge-inverted region has no signal, gets DIO only for e- tracks, and models both RMC
sources. Process and nuisance names are shared by all channels, so the control regions constrain them.

Process shapes (the default PDF type 'auto', see --functions and the FIT_* environment overrides below):
  signal      the MC histogram (--functions: a double-sided Crystal Ball fit, up to 107 / 97 MeV/c)
  dio         the MC histogram, with the tail above 102 MeV/c replaced by the generated spectrum convolved with
              the momentum response (set 10)
  cosmic      Chebychev (2 coefficients) fit to the set plus its CRV-tagged version (set + 1000); rate from the set
  pbar, rpc   Chebychev fits (pbar, rpc_ext: 1 coefficient; rpc_int: 2)
  rmc         mu- --> e+: exp(quadratic), i.e. a broad Gaussian (1n), and a Gaussian + linear polynomial (0n); mu- --> e-: a turn-on x
              power law x exponential shape
  A non-RMC background filling less than 75% of the bins is modeled as uniform.

RMC with a charge-inverted control region (offset 4000): the control region measures the RMC, so its rate floats
(rmcNorm, a free scale of every RMC process in every channel, in place of rmcN and lumi), and its photon spectrum power
(RMCPower, shape only) has a flat prior. intN keeps the internal/external ratio. Without one, the RMC rate has the
rmcN prior and no power uncertainty (--rmc-power sets the power treatment explicitly). Fitted shapes are used as binned
  templates; with --functions they enter the workspace as the fitted PDFs.

Environment (defaults of the options, as used by the macro-based workflow):
  CONVANA_HIST_PATH  histogram directory          CONVANA_HIST_PDFS=0  --functions
  CONVANA_USE_ENV=1  --envelope                   CONVANA_USE_2BB=1    --2bb
  FIT_PDF_TYPE[_<COMP>]      auto, analytic (auto without the uniform fallback), hist, uniform, polyN, cb
  FIT_TAIL_MODEL[_DIO]       convolution (default) or none
  FIT_SHAPE_SETS[_<COMP>]    comma-separated sets whose summed histograms give the shape (the rate stays the set's)
  FIT_CONTROL_SETS[_COSMIC]  sets added to the cosmic shape fit (default: set + 1000)
  <COMP>: SIGNAL DIO COSMIC RPC_EXT RPC_INT PBAR RMC_EXT RMC_INT, and RMC_EXT_0N ... RMC_INT_1N for mu- --> e+
"""

import argparse
import copy
import glob
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def import_mu2ecombine():
    """mu2ecombine from PYTHONPATH, $MU2ECOMBINE_DIR, or a Mu2eEvtAna clone next to this repository."""
    try:
        import mu2ecombine
        return mu2ecombine
    except ImportError:
        pass
    candidates = [os.environ.get("MU2ECOMBINE_DIR", ""),
                  os.path.normpath(os.path.join(HERE, "..", "..", "..", "Mu2eEvtAna", "combine"))]
    for path in candidates:
        if path and os.path.isdir(os.path.join(path, "mu2ecombine")):
            sys.path.insert(0, path)
            import mu2ecombine
            return mu2ecombine
    sys.exit("build_model.py: mu2ecombine not found. Add <Mu2eEvtAna>/combine to PYTHONPATH or set MU2ECOMBINE_DIR")


mc = import_mu2ecombine()
import ROOT  # noqa: E402  (after mu2ecombine, which may extend sys.path)

from samples import Exposure, make_samples, rate_uncertainties  # noqa: E402

ROOT.gROOT.SetBatch(True)

#----------------------------------------------------------------------------------------------------------------------
# Regions

CRV_INVERTED, TIME_INVERTED, CHARGE_INVERTED = 1000, 2000, 4000
BIN_WIDTH = 0.2
FIT_RANGE = {"mumem": (100., 110.), "mumep": (87., 97.)}
FIT_RANGE_CHARGE_INVERTED = {"mumem": (95., 105.), "mumep": (87., 97.)}
BLIND = {"mumem": (102., 105.5), "mumep": (90., 93.)}   # excluded from the envelope fits
SIGNAL_FIT_MAX = {"mumem": 107., "mumep": 97.}
DIO_TAIL_START = 102.
SPARSE_FRACTION = 0.75
COSMIC_CR_RATIO = {True: 754., False: 471.}  # CRV-tagged / vetoed cosmic yield, e- / e+ tracks
# RMC samples with the Mu2eEvtAna RMCPower variations (photon spectrum x(1-x)^n, n = 2 +- 0.5): the 0 neutron
# knockout samples. The power is a shape-only uncertainty: its up/down templates keep each process's nominal yield in
# every channel, so it is separate from the RMC rate
RMC_POWER_SAMPLES = {"mumem": ("rmc_ext", "rmc_int"), "mumep": ("rmc_ext_0n", "rmc_int_0n")}


class Region:
    def __init__(self, process, selection):
        self.selection = selection
        self.base = selection % 1000
        flags = selection // 1000
        self.crv_inverted = bool(flags & (CRV_INVERTED // 1000))
        self.time_inverted = bool(flags & (TIME_INVERTED // 1000))
        self.charge_inverted = bool(flags & (CHARGE_INVERTED // 1000))
        self.negative_tracks = (process == "mumem") != self.charge_inverted
        self.has_signal = not self.charge_inverted
        ranges = FIT_RANGE_CHARGE_INVERTED if self.charge_inverted else FIT_RANGE
        self.xmin, self.xmax = ranges[process]

    @property
    def is_control(self):
        return self.crv_inverted or self.time_inverted or self.charge_inverted


TITLES = {"signal": {"mumem": "#mu^{-}#rightarrowe^{-}", "mumep": "#mu^{-}#rightarrowe^{+}"}, "dio": "DIO",
          "cosmic": "Cosmic ray", "pbar": "Antiproton", "rpc_ext": "RPC (external)", "rpc_int": "RPC (internal)",
          "rmc_ext": "RMC (external)", "rmc_int": "RMC (internal)", "rmc_ext_0n": "RMC (0n-ext)",
          "rmc_ext_1n": "RMC (1n-ext)", "rmc_int_0n": "RMC (0n-int)", "rmc_int_1n": "RMC (1n-int)",
          "env": "Background envelope"}


def title(process, component):
    t = TITLES.get(component, component)
    return t[process] if isinstance(t, dict) else t


def background_components(process, region):
    """Background processes of a region; the backgrounds follow the charge of the selected tracks."""
    comps = ["dio"] if region.negative_tracks else []
    comps += ["cosmic", "rpc_ext", "rpc_int", "pbar"]
    if process == "mumep":
        comps += ["rmc_ext_0n", "rmc_ext_1n", "rmc_int_0n", "rmc_int_1n"]
    else:
        comps += ["rmc_ext"] + (["rmc_int"] if region.charge_inverted else [])
    return comps


#----------------------------------------------------------------------------------------------------------------------
# Configuration from the environment

def env_flag(name, default):
    value = os.environ.get(name)
    if value is None:
        return default
    if value.lower() in ("1", "true", "yes"):
        return True
    if value.lower() in ("0", "false", "no"):
        return False
    print(f"build_model.py: unrecognized value {name}={value}, using {default}")
    return default


def env_setting(prefix, component, default=""):
    """FIT_<prefix>_<COMPONENT>, else FIT_<prefix>, else the default."""
    return (os.environ.get(f"FIT_{prefix}_{component.upper()}") or os.environ.get(f"FIT_{prefix}") or default).strip()


def env_sets(prefix, component):
    value = env_setting(prefix, component)
    return [int(s) for s in value.replace(",", " ").split()] if value else []


#----------------------------------------------------------------------------------------------------------------------
# Shapes

KEEP = []  # RooFit objects built here must outlive the model


def filled_fraction(h, edges):
    _, values, _ = mc.hist.arrays(mc.rebin_to(h, edges))
    return np.count_nonzero(values > 0.) / len(values)


def analytic_pdf(process, component, region, obs, name):
    """The analytic shape of a background: a RooAbsPdf on obs, or (predefined mu2ecombine function, order)."""
    if component == "cosmic":
        return ("chebychev", 2)
    if component in ("pbar", "rpc_ext"):
        return ("chebychev", 1)
    if component == "rpc_int":
        return ("chebychev", 2)
    if component.startswith("rmc"):
        if process == "mumep":
            if component.endswith("_1n"):
                # a broad Gaussian peaking below the range, as exp(quadratic): the fit lands far in the tail, where
                # Gaussian bin integrals (erf differences) lose all precision
                return ("exp_poly", 2)
            # a broad Gaussian plus a linear polynomial (0n)
            mean = ROOT.RooRealVar(f"{name}_mean", "Gaussian mean", 75., 30., 90.)
            sigma = ROOT.RooRealVar(f"{name}_sigma", "Gaussian sigma", 12., 5., 20.)
            KEEP.extend([mean, sigma])
            p0 = ROOT.RooRealVar(f"{name}_p0", "constant", -0.001, -1., 1.)
            p1 = ROOT.RooRealVar(f"{name}_p1", "slope", 0., -0.1, 0.1)
            o = obs.GetName()
            pdf = ROOT.RooGenericPdf(name, "RMC shape", f"max(0., TMath::Gaus({o}, {mean.GetName()}, {sigma.GetName()})"
                                     f" + {p0.GetName()} + {p1.GetName()}*{o}/90)",
                                     ROOT.RooArgList(obs, mean, sigma, p0, p1))
            KEEP.extend([p0, p1, pdf])
            return pdf
        x0 = ROOT.RooRealVar(f"{name}_x0", "Low threshold edge", 94. if region.charge_inverted else 96.)
        x0.setConstant(True)
        alpha = ROOT.RooRealVar(f"{name}_alpha", "Low-edge turn-on power", 1.915, 0.1, 4.0)
        beta = ROOT.RooRealVar(f"{name}_beta", "Bulk curvature power", 0.457, 0.01, 20.0)
        lam = ROOT.RooRealVar(f"{name}_lambda", "Tail exponential decay", 1.932, 0.01, 5.0)
        o = obs.GetName()
        pdf = ROOT.RooGenericPdf(name, "RMC shape",
                                 f"pow(max(0., {o} - {x0.GetName()}), {alpha.GetName()}) * pow(120.0 - {o}, "
                                 f"{beta.GetName()}) * exp(-{lam.GetName()} * {o})",
                                 ROOT.RooArgList(obs, x0, alpha, beta, lam))
        KEEP.extend([x0, alpha, beta, lam, pdf])
        return pdf
    return None


def fit_shape(function, fit_hists, rate, binned, fit_range=None):
    if isinstance(function, tuple):
        name, order = function
        return mc.Fit(name, fit_hists, rate=rate, order=order, binned=binned, fit_range=fit_range)
    return mc.Fit(function, fit_hists, rate=rate, binned=binned, fit_range=fit_range)


def component_shape(args, samples, region, component, edges, obs, channel):
    """(shape, rate histogram) of one process in one channel, or None if it has no expected events."""
    process, sel = args.process, region.selection
    key = args.process if component == "signal" else component
    hist = lambda s, var="obs": samples[key].hist(f"trk_{s}/{var}")
    in_channel = lambda s: mc.rebin_to(hist(s), edges)
    h = in_channel(sel)
    rate_hist = h
    label = f"{channel}/{component}"

    shape_sets = env_sets("SHAPE_SETS", component)
    shape_hists = [in_channel(s) for s in shape_sets] if shape_sets else [h]
    if component == "cosmic":
        control_sets = env_sets("CONTROL_SETS", component)
        if not control_sets and not region.crv_inverted:
            control_sets = [sel + CRV_INVERTED]
        shape_hists = shape_hists + [in_channel(s) for s in control_sets]

    rate = mc.hist.integral(h)
    if component == "cosmic" and rate <= 0.:
        # no simulated cosmic rays pass the veto: the rate from the CRV-tagged region
        rate = sum(mc.hist.integral(x) for x in shape_hists[1:]) / COSMIC_CR_RATIO[
            region.negative_tracks]
        print(f"  {label}: no events in the set, rate {rate:.4g} from the CRV-tagged region")
    if rate <= 0.:
        print(f"  {label}: no expected events in {region.xmin:g}-{region.xmax:g}, not modeled")
        return None

    pdf_type = env_setting("PDF_TYPE", component, "auto").lower()
    binned = not args.functions

    if component == "signal":
        if not args.functions and pdf_type in ("auto", "analytic"):
            return mc.Template(h, shape=shape_hists), rate_hist
        if pdf_type in ("auto", "analytic", "cb"):
            return mc.Fit("dscb", shape_hists, rate=h, binned=binned,
                          fit_range=(region.xmin, min(region.xmax, SIGNAL_FIT_MAX[process]))), rate_hist
    if component == "dio":
        if env_setting("TAIL_MODEL", component, "convolution").lower() == "convolution" and \
           region.xmax > DIO_TAIL_START:
            reference = mc.convolved(hist(10, "MC_GenE"), hist(10, "dP"), edges)
            h = mc.replace_tail(h, DIO_TAIL_START, reference, rel_error=0.2)
            shape_hists = [h] if not shape_sets else shape_hists
            rate_hist = h
        if pdf_type in ("auto", "analytic"):
            return mc.Template(h, shape=shape_hists), rate_hist

    rate_arg = h if component != "cosmic" or mc.hist.integral(h) > 0. else rate
    if pdf_type == "hist":
        if not isinstance(rate_arg, float):
            return mc.Template(h, shape=shape_hists), rate_hist
        return mc.Fit("uniform", shape_hists, rate=rate_arg, binned=binned), rate_hist
    if pdf_type == "uniform":
        return mc.Fit("uniform", shape_hists, rate=rate_arg, binned=binned), rate_hist
    if pdf_type.startswith("poly"):
        order = int(pdf_type[4:] or 1) + 1
        return mc.Fit("chebychev", shape_hists, rate=rate_arg, order=order, binned=binned), rate_hist
    if pdf_type == "cb":
        return mc.Fit("dscb", shape_hists, rate=rate_arg, binned=binned), rate_hist
    if pdf_type not in ("auto", "analytic"):
        sys.exit(f"build_model.py: unknown PDF type '{pdf_type}' for {component}")

    if pdf_type == "auto" and not component.startswith("rmc") and filled_fraction(h, edges) < SPARSE_FRACTION:
        print(f"  {label}: only {filled_fraction(h, edges):.0%} of the bins are filled, using a uniform shape")
        return mc.Fit("uniform", shape_hists, rate=rate_arg, binned=binned), rate_hist
    function = analytic_pdf(process, component, region, obs, f"{channel}_{component}_shape")
    if function is None:
        return mc.Template(h, shape=shape_hists), rate_hist
    return fit_shape(function, shape_hists, rate_arg, binned), rate_hist


#----------------------------------------------------------------------------------------------------------------------
# Data

def data_hists(args, samples, regions):
    """{channel: TH1} for the --data / --data-model choice; channels missing from the dict use the Asimov data."""
    out = {}
    if args.data_model:
        info = json.load(open(os.path.join(args.data_model, "model.json")))
        for ch, region in regions.items():
            if ch not in info["channels"]:
                sys.exit(f"build_model.py: channel {ch} is not in the model {args.data_model}")
            out[ch] = model_data(args.data_model, ch, info["channels"][ch])
        return out
    if args.data in ("asimov",) or args.data.startswith("toy"):
        return out
    if args.data in ("mds3c", "mds3d"):
        for ch, region in regions.items():
            out[ch] = samples[f"data_{args.data}"].hist(f"trk_{region.selection}/obs")
        return out
    if ":" not in args.data:
        sys.exit(f"build_model.py: unknown --data '{args.data}' (asimov, toy[:seed], mds3c, mds3d, or "
                 "<file>:<histogram path>, where {set} is replaced by the channel's set)")
    path, name = args.data.split(":", 1)
    f = ROOT.TFile.Open(path)
    for ch, region in regions.items():
        h = f.Get(name.format(set=region.selection)) if f else None
        if not h:
            sys.exit(f"build_model.py: no histogram {name.format(set=region.selection)} in {path}")
        out[ch] = mc.hist.copy(h)
    return out


def model_data(model_dir, channel, info):
    """The data (data_obs) of a channel of a model written by mu2ecombine, as a TH1D."""
    path = os.path.join(model_dir, "shapes", f"{channel}.root")
    f = ROOT.TFile.Open(path)
    edges = np.array(info["observable"]["edges"])
    if info["mode"] == "binned":
        return mc.hist.copy(f.Get("data_obs"))
    data = f.Get("w").data(f"{channel}_data_obs")
    values = np.array([(data.get(i), data.weight())[1] for i in range(data.numEntries())])
    return mc.hist.make(edges, values, values, "data_obs")


#----------------------------------------------------------------------------------------------------------------------
# Model

def channel_name(region):
    offset = region.selection - region.base
    return f"set_{region.base}" if offset == 0 else f"cr{offset}_{region.base}"


def new_model(args, name):
    # workspaces on the channel observables even for template channels, so that every model of a channel (templates,
    # functions, envelope) can generate data for the others
    model = mc.Model(name, poi_range=tuple(args.poi_range), workspace=True)
    model.mc_stats(enabled=args.mc_stat)
    return model


def common_terms(args, model, channels, processes):
    """Rate uncertainties, the mcstat group, and frozen yield scale factors (set them with --setParameters, e.g.
    yieldScale=2 or scale_dio=0.5)."""
    rate_uncertainties(model, rmc_floating=args.rmc_floating)
    for ch in channels:
        for p in processes:
            model.group("mcstat", f"mcstat_{ch}_{p}")
    model.rate_param("yieldScale", "*", range=(0., 100.))
    for p in sorted(processes):
        model.rate_param(f"scale_{p}", p, range=(0., 100.))


def channel_edges(region):
    return np.round(np.arange(region.xmin, region.xmax + BIN_WIDTH / 2, BIN_WIDTH), 6)


def define_model(args, samples, regions, data=None, envelope=True):
    """The mu2ecombine Model, and the rate histogram of each process, {channel: {process: TH1}}."""
    model = new_model(args, args.name)
    data = data or {}
    inputs = {}
    processes = set()
    for ch, region in regions.items():
        edges = channel_edges(region)
        obs = mc.observable(f"obs_{region.selection}", edges=edges, title="p [MeV/c]")
        KEEP.append(obs)
        use_envelope = envelope and args.envelope and not region.is_control
        channel = model.channel(ch, obs, data=data.get(ch), blind=BLIND[args.process] if region.has_signal else None)
        print(f"Channel {ch}: set {region.selection}, {region.xmin:g}-{region.xmax:g} MeV/c")
        inputs[ch] = {}
        components = (["signal"] if region.has_signal else []) + background_components(args.process, region)
        for comp in components:
            result = component_shape(args, samples, region, comp, edges, obs, ch)
            if result is None:
                continue
            shape, inputs[ch][comp] = result
            if args.rmc_power != "off" and comp in RMC_POWER_SAMPLES[args.process]:
                up, down = samples[comp].variation("RMCPower", region.selection)
                shape.variations = mc.shapes.as_variations({"RMCPower": (mc.rebin_to(up, edges),
                                                                         mc.rebin_to(down, edges), "shape_only")})
            channel.add(comp, shape, signal=comp == "signal", title=title(args.process, comp))
            processes.add(comp)
        if use_envelope:
            # the envelope absorbs every MC background; it is only used where the signal is fit
            channel.add("env", mc.Envelope(absorbs="*"), title=title(args.process, "env"))
            processes.add("env")
    common_terms(args, model, regions, processes)
    if args.rmc_power == "free":
        model.unconstrained("RMCPower")
    return model, inputs


def narrowed_regions(args, regions):
    """The regions with the signal-region channels restricted to --range, moved to the nearest bin edges."""
    out = {}
    for ch, region in regions.items():
        narrow = copy.copy(region)
        if not region.is_control:
            edges = channel_edges(region)
            lo, hi = (float(edges[np.argmin(np.abs(edges - x))]) for x in args.range)
            if hi - lo < BIN_WIDTH - 1e-6:
                sys.exit(f"build_model.py: the range {args.range[0]:g}-{args.range[1]:g} is empty in channel {ch}")
            narrow.xmin, narrow.xmax = lo, hi
        out[ch] = narrow
    return out


def narrowed_model(args, full, inputs, regions):
    """The model on the narrowed regions, with the full-range model's shapes: each process is the full-range model's
    expectation in the bins inside the range, so only the range changes. The MC uncertainty of each bin is the
    relative uncertainty of the process's input histogram in that bin. An envelope is fit to the data in the range."""
    model = new_model(args, args.name)
    processes = set()
    narrow_inputs = {}
    for ch, r in full.items():
        region = regions[ch]
        edges = r.channel.edges
        lo = int(np.argmin(np.abs(edges - region.xmin)))
        hi = int(np.argmin(np.abs(edges - region.xmax)))
        sub = edges[lo:hi + 1]
        obs = mc.observable(f"obs_{region.selection}", edges=sub, title="p [MeV/c]")
        KEEP.append(obs)
        data = None
        if r.data_source == "data":
            _, values, variances = mc.hist.arrays(r.data)
            data = mc.hist.make(sub, values[lo:hi], variances[lo:hi], "data_obs")
        channel = model.channel(ch, obs, data=data, blind=r.channel.blind)
        print(f"Channel {ch}: set {region.selection}, {region.xmin:g}-{region.xmax:g} MeV/c (full-range shapes)")
        narrow_inputs[ch] = {}
        for p in r.processes:
            values = r.expected[p][lo:hi]
            rel2 = np.zeros_like(values)
            if p in inputs[ch]:
                _, in_values, in_variances = mc.hist.arrays(inputs[ch][p])
                ok = in_values[lo:hi] > 0.
                rel2[ok] = in_variances[lo:hi][ok] / in_values[lo:hi][ok] ** 2
            h = mc.hist.make(sub, values, values ** 2 * rel2, p)
            if h.Integral() <= 0.:
                print(f"  {ch}/{p}: no expected events in {region.xmin:g}-{region.xmax:g}, not modeled")
                continue
            channel.add(p, h, signal=r.channel.processes[p].signal, title=r.channel.processes[p].title)
            narrow_inputs[ch][p] = h
            processes.add(p)
        if args.envelope and not region.is_control:
            channel.add("env", mc.Envelope(absorbs="*"), title=title(args.process, "env"))
            processes.add("env")
    common_terms(args, model, full, processes)
    return model, narrow_inputs


def toy_data(args, samples, regions):
    """Background-only pseudo-data: the expected backgrounds (MC models, no envelope) Poisson-fluctuated per bin."""
    seed = int(args.data.split(":", 1)[1]) if ":" in args.data else 123456
    rng = np.random.default_rng(seed)
    model, _ = define_model(args, samples, regions, envelope=False)
    out = {}
    for ch, r in model.build().items():
        bkg = sum(r.expected[p] for p in r.expected if not r.channel.processes[p].signal)
        out[ch] = mc.hist.make(r.channel.edges, rng.poisson(bkg).astype(float), None, "data_obs")
        print(f"  {ch}: toy data (seed {seed}), {out[ch].Integral():.0f} events for {bkg.sum():.1f} expected")
    return out


def window(signal, background, optimize):
    """The bin range [lo, hi) maximizing S / sqrt(S + B), or the full range."""
    n = len(background)
    if not optimize:
        return 0, n
    best, best_value = (0, n), -1.
    s_cum = np.concatenate([[0.], np.cumsum(signal)])
    b_cum = np.concatenate([[0.], np.cumsum(background)])
    for lo in range(n):
        s = s_cum[lo + 1:] - s_cum[lo]
        b = b_cum[lo + 1:] - b_cum[lo]
        ok = (s > 0.) & (b > 0.)
        if not ok.any():
            continue
        value = np.where(ok, s / np.sqrt(np.where(ok, s + b, 1.)), -1.)
        i = int(np.argmax(value))
        if value[i] > best_value:
            best_value, best = value[i], (lo, lo + i + 1)
    return best


def counting_model(args, results, inputs, regions):
    """Cut-and-count model: one bin per channel, the momentum window maximizing S / sqrt(S + B) in a signal region
    and the full range in a control region, with the shape model's expected yields in it."""
    model = new_model(args, f"{args.name}_cc")
    windows = {}
    processes = set()
    for ch, r in results.items():
        region = regions[ch]
        procs = r.processes
        signal = sum((r.expected[p] for p in procs if r.channel.processes[p].signal), np.zeros(len(r.channel.edges) - 1))
        bkg = sum((r.expected[p] for p in procs if not r.channel.processes[p].signal), np.zeros_like(signal))
        lo, hi = window(signal, bkg, region.has_signal and not region.is_control)
        edges = r.channel.edges
        xlo, xhi = float(edges[lo]), float(edges[hi])
        windows[ch] = (xlo, xhi)
        obs = mc.observable(f"obs_{region.selection}_cc", edges=[xlo, xhi], title="p [MeV/c]")
        KEEP.append(obs)
        data = None
        if r.data_source == "data":
            _, values, _ = mc.hist.arrays(r.data)
            data = mc.hist.make([xlo, xhi], [values[lo:hi].sum()], None, "data_obs")
        channel = model.channel(ch, obs, data=data)
        for p in procs:
            yield_ = float(r.expected[p][lo:hi].sum())
            rel = 0.
            if p in inputs[ch]:
                _, values, variances = mc.hist.arrays(inputs[ch][p])
                total = values[lo:hi].sum()
                rel = np.sqrt(variances[lo:hi].sum()) / total if total > 0. else 0.
            channel.add(p, mc.hist.make([xlo, xhi], [yield_], [(rel * yield_) ** 2], p),
                        signal=r.channel.processes[p].signal, title=r.channel.processes[p].title)
            processes.add(p)
        print(f"  {ch}: counting window {xlo:.2f}-{xhi:.2f} MeV/c, N(signal) = {signal[lo:hi].sum():.4f}, "
              f"N(background) = {bkg[lo:hi].sum():.4f}")
    common_terms(args, model, results, processes)
    return model, windows


def write(model, outdir, plots, results=None):
    """Write the model, reusing already built results (Model.write builds the model itself otherwise)."""
    if results is not None:
        model.build = lambda: results
    return model.write(outdir, plots=plots)


def finish_cards(outdir, header):
    """Header comments (read by sensitivity_vs_livetime.py) and frozen yield scale factors in every card."""
    for card in glob.glob(os.path.join(outdir, "*.txt")):
        text = open(card).read()
        # only the yield scale factors are frozen; a fitted rate (rmcNorm) floats
        params = re.findall(r"^(yieldScale|scale_\S+) rateParam ", text, flags=re.M)
        freeze = [f"nuisance edit freeze {p}" for p in dict.fromkeys(params)]
        with open(card, "w") as f:
            f.write("\n".join(header) + "\n" + text.rstrip("\n") + "\n" + "\n".join(freeze) + "\n")


def make_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--process", default="mumem", choices=("mumem", "mumep"))
    parser.add_argument("-s", "--selections", default="20", help="selection sets, quoted (default: %(default)s)")
    parser.add_argument("-c", "--control-regions", default="",
                        help="control region offsets added to each selection set, quoted, e.g. \"1000 4000\"")
    parser.add_argument("-t", "--tag", default="", help="output name tag")
    parser.add_argument("-o", "--output", default=None,
                        help="output directory (default: models/<process>_<sets>[_cr<offset>...][_<tag>])")
    parser.add_argument("--exposure", default=None, choices=Exposure.NAMES,
                        help="exposure (default: run1a, or the data ensemble of --data mds3c/mds3d)")
    parser.add_argument("--2bb", dest="two_bb", action="store_true", default=env_flag("CONVANA_USE_2BB", False),
                        help="2BB exposure, with the inputs re-weighted to 2BB (ConvAna.*.m2.root)")
    parser.add_argument("--reference", type=float, default=None,
                        help="branching fraction the signal is normalized to, the unit of r "
                             "(default: 1e-15, 1e-13 for the MDS ensembles)")
    parser.add_argument("--functions", action="store_true", default=not env_flag("CONVANA_HIST_PDFS", True),
                        help="fitted PDFs in the workspace instead of binned templates")
    parser.add_argument("--envelope", action="store_true", default=env_flag("CONVANA_USE_ENV", False),
                        help="a data-driven envelope (discrete profiling) absorbs the backgrounds of the signal "
                             "regions")
    parser.add_argument("--data", default="asimov",
                        help="asimov (default), toy[:seed] (Poisson-fluctuated backgrounds), mds3c, mds3d, or "
                             "<file>:<histogram>, with {set} replaced by each channel's set")
    parser.add_argument("--data-model", default=None,
                        help="use the data of another model directory, e.g. the Asimov data of a template model "
                             "for an envelope model")
    parser.add_argument("--mc-stat", action="store_true", help="MC statistical uncertainties of the yields")
    parser.add_argument("--poi-range", type=float, nargs=2, default=(-50., 50.), metavar=("LO", "HI"))
    parser.add_argument("--charge-inverted-range", type=float, nargs=2, default=None, metavar=("LO", "HI"),
                        help="momentum range of the charge-inverted control regions (default: 95-105 MeV/c for "
                             "mumem, 87-97 for mumep)")
    parser.add_argument("--rmc-power", default="auto", choices=("auto", "off", "constrained", "free"),
                        help="the RMC photon spectrum power uncertainty (Mu2eEvtAna RMCPower, n = 2 +- 0.5, shape "
                             "only) on the 0n RMC processes of every channel: none, a unit-Gaussian (constrained) or "
                             "a flat (free) prior; auto (default): free with a charge-inverted control region, which "
                             "measures it, and off otherwise")
    parser.add_argument("--range", type=float, nargs=2, default=None, metavar=("LO", "HI"),
                        help="restrict the signal-region channels to this momentum range (moved to the nearest bin "
                             "edges), keeping the full-range shapes of the MC processes")
    parser.add_argument("--hist-path", default=None, help="histogram directory (default: CONVANA_HIST_PATH, or the "
                                                          "analysis histogram directory)")
    parser.add_argument("--plots", action="store_true", help="expected stacks and fit figures in <output>/figures")
    parser.add_argument("--no-cc", action="store_true", help="skip the cut-and-count model")
    parser.add_argument("--only", default=None, metavar="PROCESS",
                        help="develop/debug one process: build and fit only it (e.g. signal, cosmic, rmc_ext_0n, env), "
                             "print the result, draw it to <output>/debug/, and exit without writing a model")
    parser.add_argument("--channel", default=None, help="with --only: the channel (default: every channel with the "
                                                       "process)")
    return parser


def parse_args(argv=None):
    """The options, with the model name (args.name) and output directory (args.output) filled in."""
    parser = make_parser()
    args = parser.parse_args(argv)

    selections = [int(s) for s in args.selections.split()]
    offsets = [int(o) for o in args.control_regions.split()]
    for o in offsets:
        if o % 1000 or not 1 <= o // 1000 <= 7:
            parser.error(f"control region offset {o} is not a sum of the flags 1000, 2000, and 4000")
    name = "_".join([args.process] + [str(s) for s in selections] + [f"cr{o}" for o in offsets]
                    + ([args.tag] if args.tag else []))
    args.name = name
    args.output = args.output or os.path.join("models", name)
    args.selection_list, args.offsets = selections, offsets
    # a charge-inverted (RMC) control region measures the RMC: its rate (rmcNorm) and spectrum power float
    args.rmc_floating = any((x // 1000) & (CHARGE_INVERTED // 1000) for x in offsets + selections)
    if args.rmc_power == "auto":
        args.rmc_power = "free" if args.rmc_floating else "off"
    return args


def run(args):
    """Build and write the shape and cut-and-count models of the parsed options. Returns the output directory."""
    outdir = args.output
    exposure_name = args.exposure or (args.data if args.data in ("mds3c", "mds3d") else "run1a")
    exposure = Exposure(exposure_name, args.two_bb)
    if args.reference:
        exposure.reference = args.reference
    if args.data in ("mds3c", "mds3d") and exposure_name != args.data:
        print(f"build_model.py: WARNING: the {args.data} data with the {exposure_name} exposure")
    samples = make_samples(exposure, args.hist_path)
    print(exposure, f"signal reference {exposure.reference:.1e}")

    regions = {}
    for s in args.selection_list:
        for sel in [s] + [s + o for o in args.offsets]:
            region = Region(args.process, sel)
            if args.charge_inverted_range and region.charge_inverted:
                region.xmin, region.xmax = args.charge_inverted_range
            regions[channel_name(region)] = region

    data = toy_data(args, samples, regions) if args.data.startswith("toy") else data_hists(args, samples, regions)
    if args.only:
        model, _ = define_model(args, samples, regions, data)
        channels = [args.channel] if args.channel else [ch for ch, c in model.channels.items()
                                                         if args.only in c.processes]
        if not channels:
            sys.exit(f"build_model.py: no channel has the process '{args.only}'")
        for ch in channels:
            model.fit_process(ch, args.only, figdir=os.path.join(outdir, "debug"))
        return outdir
    if args.range:
        full, full_inputs = define_model(args, samples, regions, data, envelope=False)
        regions = narrowed_regions(args, regions)
        model, inputs = narrowed_model(args, full.build(), full_inputs, regions)
    else:
        model, inputs = define_model(args, samples, regions, data)
    results = model.build()
    write(model, outdir, args.plots, results)

    n_signal = exposure.nmuons * mc.mu2e_physics.muon_capture_fraction_al * exposure.reference
    header = ["# Combine model written by analysis/combine/build_model.py",
              f"# R_mue used for signal: {exposure.reference:.3e}",
              f"# N(POT): {exposure.npot:.3e}",
              f"# Livetime: {exposure.livetime:.3e}",
              f"# N(muons): {exposure.nmuons:.3e}"]
    for ch, r in results.items():
        if "signal" in r.processes:
            header.append(f"# Signal efficiency ({ch}): {r.card_rates['signal'] / n_signal:.3e}")
    finish_cards(outdir, header)

    if not args.no_cc:
        print("Cut-and-count model:")
        cc_model, windows = counting_model(args, results, inputs, regions)
        cc_dir = os.path.join(outdir, "cc")
        write(cc_model, cc_dir, False)
        finish_cards(cc_dir, header + [f"# Window ({ch}): {lo:.2f} - {hi:.2f} MeV/c" for ch, (lo, hi)
                                       in windows.items()])
    print(f"Model written to {outdir}: combined.txt" + ("" if args.no_cc else ", cc/combined.txt (cut-and-count)"))
    return outdir


if __name__ == "__main__":
    run(parse_args())
