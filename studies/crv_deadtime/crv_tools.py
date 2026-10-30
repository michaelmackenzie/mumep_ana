"""Shared tools for the CRV veto timing window study.

Reads the ConvAna track trees (trs_<set>) of the e+ selection with and without the ConvAna CRV veto (set and
set + 1000, the kCRVVetoOffset control region), normalizes each sample with analysis/optimize_selection.C, and scans
veto windows lo < dt < hi on the track - CRV cluster time differences:
  trk_crv_dt:     track t0 - CRV cluster time propagated to the tracker (via the calorimeter front or the stopping
                  target, whichever is closer in time), for the cluster closest in time
  trk_crv_raw_dt: track t0 - raw CRV cluster time, for the cluster closest in time
  trk_crv_cal_dt: track t0 - CRV cluster time propagated via the calorimeter front, for the cluster closest in time
  trk_crv_st_dt:  track t0 - CRV cluster time propagated via the stopping target, for the cluster closest in time
  trk_upstream_dt: track t0 - t0 of the upstream-going track closest in time among those > 25 ns earlier (cosmic
                  rays reflected in the solenoid field)
Tracks without a CRV cluster (upstream track) have -10000 and are never vetoed. trk_crv_dt is the one of trk_crv_cal_dt
and trk_crv_st_dt with the smaller |dt|.

The upstream window can only be studied with ConvAna set 43, which has no CRV, upstream, or multi-track veto; the
vetoes are then emulated from the stored variables. Set 40 (and 1040) already has the upstream veto applied.
"""

import os
from pathlib import Path

import numpy as np
import uproot

STUDY_DIR    = Path(__file__).resolve().parent
REPO_DIR     = STUDY_DIR.parents[1]                     # mumep_ana/studies/crv_deadtime
ANALYSIS_DIR = REPO_DIR / "analysis"
FIG_DIR      = STUDY_DIR / "figures"
OUT_DIR      = STUDY_DIR / "outputs"

CRV_VETO_OFFSET = 1000 # ConvAna kCRVVetoOffset: set + 1000 holds the tracks that fail the CRV veto
NO_VETO_SET     = 43     # ConvAna set 40 without the CRV, upstream, and multi-track vetoes
NO_CRV          = -9999. # values below this have no CRV cluster (or no upstream track)
MULTITRK_DT     = 150.   # ConvAna multi-track veto: another downstream track within this time

# Track - CRV time variables: axis labels, scan range, and the ConvAna veto window (Mu2eEvtAna SetTrackIDs)
VARIABLES = {
    "trk_crv_dt"     : dict(label=r"$t_\mathrm{track} - t_\mathrm{CRV}$ (propagated) [ns]"  , scan=(-150., 150.),
                            current=(-50., 60.)),
    "trk_crv_raw_dt" : dict(label=r"$t_\mathrm{track} - t_\mathrm{CRV}$ (raw) [ns]"         , scan=(-100., 300.),
                            current=(-25., 0.)),
    "trk_crv_cal_dt" : dict(label=r"$t_\mathrm{track} - t_\mathrm{CRV}$ (via calo) [ns]"    , scan=(-150., 150.),
                            current=(-50., 60.)),
    "trk_crv_st_dt"  : dict(label=r"$t_\mathrm{track} - t_\mathrm{CRV}$ (via target) [ns]"  , scan=(-150., 150.),
                            current=(-50., 60.)),
    "trk_upstream_dt": dict(label=r"$t_\mathrm{track} - t_\mathrm{upstream\ track}$ [ns]"  , scan=(25., 300.),
                            current=(40., 110.)),
}
UPSTREAM = "trk_upstream_dt"
# The ConvAna veto in terms of the stored variables: the propagated window applies to whichever hypothesis is closer,
# for the matched cluster (closest in propagated time), and the raw window to that same cluster
CONVANA_VETO = {"trk_crv_dt": (-50., 60.), "trk_crv_raw_dt": (-25., 0.)}
# The ConvAna upstream (reflection) veto: any upstream-going track in this window
CONVANA_UPSTREAM_VETO = {UPSTREAM: (40., 110.)}
# Pairs of variables to compare in 2D
PAIRS = [("trk_crv_dt", "trk_crv_raw_dt"), ("trk_crv_cal_dt", "trk_crv_st_dt")]

# Samples (analysis/datasets.C names): mu- --> e+ signal and the backgrounds in the e+ momentum window
SIGNAL      = "mumep"
BACKGROUNDS = ["cosmic", "rmc_ext_0n", "rmc_ext_1n", "rmc_int_0n", "rmc_int_1n", "rpc_ext", "rpc_int"]

# Colors: the analysis process colors (analysis/datasets.C set_style, as in studies/eminus_rmc)
BLACK, BLUE, ORANGE, GRAY = "#1F1F1F", "#2a78d6", "#eb6834", "#94a4a2"
COSMIC, RMC, RPC = "#6699ff", "#ff6666", "#66cc66"


#------------------------------------------------------------------------------------------------------
# Inputs
#------------------------------------------------------------------------------------------------------

_root = None
def analysis_macros():
    """ROOT with analysis/optimize_selection.C loaded (dataset info and normalization), set up for EventNtuple inputs."""
    global _root
    if _root is not None: return _root
    import ROOT
    ROOT.gROOT.SetBatch(True)
    if ROOT.gROOT.LoadMacro(str(ANALYSIS_DIR / "optimize_selection.C")) != 0:
        raise RuntimeError("Failed to load analysis/optimize_selection.C")
    ROOT.init_physics("run1a")
    ROOT.use_evtana_ = True
    ROOT.set_evtana_defaults()
    ROOT.init_dataset_info()
    _root = ROOT
    return ROOT


def hist_file(name):
    """ConvAna histogram file of a dataset (the path respects CONVANA_HIST_PATH, as in analysis/defaults.C)."""
    ROOT = analysis_macros()
    info = ROOT.get_dataset_info(name)
    if info.name_ == "": raise ValueError(f"Dataset {name} is not defined in analysis/datasets.C")
    path = f"{ROOT.hist_path_}ConvAna.{ROOT.hist_func_}.{info.name_}.m{ROOT.hist_mode_}.{ROOT.file_type_}"
    return path, info


def sample_norm(name, path, info):
    """Expected yield per weighted entry, as in optimize_selection.C."""
    ROOT = analysis_macros()
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie(): raise OSError(f"Can't open {path}")
    norm = ROOT.evaluate_norm(name, f, info) * info.theory_
    f.Close()
    return norm


def load_sample(name, track_set=40, pmin=90., pmax=93., charge=1):
    """Tracks in the momentum window.

    For track_set = NO_VETO_SET, the tracks of that set that pass the multi-track veto, with the ConvAna CRV and
    upstream vetoes emulated. Otherwise the tracks of track_set (pass the CRV veto) and track_set + 1000 (fail it),
    which have the upstream veto applied.

    Returns a dict of arrays: the time variables, the expected yield per track (w), convana_crv_veto and
    convana_upstream_veto (the track fails that ConvAna veto), and convana_veto (either).
    """
    path, info = hist_file(name)
    norm = sample_norm(name, path, info)
    no_veto = track_set == NO_VETO_SET
    branches = ["trk_p", "trk_charge", "weight", "trk_multitrk_dt"] + list(VARIABLES)
    sets = ((None, track_set),) if no_veto else ((False, track_set), (True, track_set + CRV_VETO_OFFSET))
    parts = []
    with uproot.open(path) as f:
        for crv_veto, s in sets:
            key = f"Ana/Hist/trs_{s}/tree"
            if key not in f: raise KeyError(f"{path} has no {key}")
            a = f[key].arrays(branches, library="np")
            sel = (a["trk_p"] > pmin) & (a["trk_p"] < pmax) & (charge*a["trk_charge"] > 0) & (a["weight"] > 0)
            if no_veto: sel &= np.abs(a["trk_multitrk_dt"]) > MULTITRK_DT # no other downstream track in time
            part = {v: a[v][sel] for v in VARIABLES}
            part["w"] = a["weight"][sel] * norm
            if no_veto:
                part["convana_crv_veto"] = veto_flags(part, CONVANA_VETO)
                part["convana_upstream_veto"] = veto_flags(part, CONVANA_UPSTREAM_VETO)
            else:
                part["convana_crv_veto"] = np.full(sel.sum(), crv_veto)
                part["convana_upstream_veto"] = np.zeros(sel.sum(), bool)
            part["convana_veto"] = part["convana_crv_veto"] | part["convana_upstream_veto"]
            parts.append(part)
    sample = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    sample["name"] = name
    return sample


def load_samples(**kwargs):
    """Signal and background samples, keyed by dataset name."""
    samples = {}
    for name in [SIGNAL] + BACKGROUNDS:
        samples[name] = load_sample(name, **kwargs)
        s = samples[name]
        print(f"{name:>11}: N = {len(s['w']):8d}, N(CRV vetoed) = {s['convana_crv_veto'].sum():5d}, "
              f"N(upstream vetoed) = {s['convana_upstream_veto'].sum():5d}, yield = {s['w'].sum():.4g}")
    return samples


def merge(samples, names):
    """Concatenate samples."""
    keys = list(VARIABLES) + ["w", "convana_veto", "convana_crv_veto", "convana_upstream_veto"]
    return {k: np.concatenate([samples[n][k] for n in names]) for k in keys}


#------------------------------------------------------------------------------------------------------
# Veto window scans
#------------------------------------------------------------------------------------------------------

class WindowYield:
    """Weighted yield of a sample in veto windows lo < x < hi, from the cumulative weighted distribution of x."""
    def __init__(self, x, w):
        order = np.argsort(x)
        self.x = x[order]
        self.cw = np.concatenate([[0.], np.cumsum(w[order])])
        self.cw2 = np.concatenate([[0.], np.cumsum(w[order]**2)])
        self.n = np.arange(len(x) + 1)
        self.total = self.cw[-1]

    def _in(self, c, lo, hi):
        lo, hi = np.asarray(lo, float), np.asarray(hi, float)
        i_lo = np.searchsorted(self.x, lo, side="right") # entries <= lo
        i_hi = np.searchsorted(self.x, hi, side="left" ) # entries <  hi
        return np.where(hi > lo, c[i_hi] - c[np.minimum(i_lo, i_hi)], 0.)

    def vetoed(self, lo, hi):    return self._in(self.cw , lo, hi)
    def vetoed_w2(self, lo, hi): return self._in(self.cw2, lo, hi)
    def vetoed_n(self, lo, hi):  return self._in(self.n  , lo, hi)


def significance(s, b):
    """S / sqrt(B)."""
    return np.where(b > 0., s / np.sqrt(np.where(b > 0., b, 1.)), 0.)


class WindowScan:
    """Signal and background yields after a veto lo < var < hi, with S / sqrt(B) relative to no veto, for all
    backgrounds (rel_sig) and for the cosmic-ray background alone (rel_sig_cosmic)."""
    def __init__(self, samples, var):
        self.var = var
        self.sig = WindowYield(samples[SIGNAL][var], samples[SIGNAL]["w"])
        self.cosmic = WindowYield(samples["cosmic"][var], samples["cosmic"]["w"])
        bkg = merge(samples, BACKGROUNDS)
        self.bkg = WindowYield(bkg[var], bkg["w"])

    def evaluate(self, lo, hi):
        """Efficiencies and S / sqrt(B) relative to no veto. The MC statistical uncertainties (_err) are those of the
        change from no veto, from the vetoed background yield."""
        s = self.sig.total - self.sig.vetoed(lo, hi)
        res = dict(eff_sig=s / self.sig.total, n_cosmic_mc=self.cosmic.vetoed_n(lo, hi))
        for tag, bkg in (("", self.bkg), ("_cosmic", self.cosmic)):
            b = bkg.total - bkg.vetoed(lo, hi)
            r = significance(s, b) / significance(self.sig.total, bkg.total)
            res[f"eff{tag or '_bkg'}"] = b / bkg.total
            res[f"rel_sig{tag}"] = r
            res[f"rel_sig{tag}_err"] = 0.5*r*np.sqrt(bkg.vetoed_w2(lo, hi))/np.where(b > 0., b, 1.)
        return res

    def grid(self, step):
        """(lo, hi) grid over the scan range, and the evaluated quantities on it (lo along axis 1)."""
        r_lo, r_hi = VARIABLES[self.var]["scan"]
        edges = np.arange(r_lo, r_hi + 0.5*step, step)
        lo, hi = np.meshgrid(edges, edges)
        return edges, lo, hi, self.evaluate(lo, hi)

    def best_window(self, step, margin=0., fom="rel_sig"):
        """Window maximizing the figure of merit (the narrowest if several do), widened by margin on each side."""
        edges, lo, hi, res = self.grid(step)
        r = np.where(hi > lo, res[fom], -np.inf)
        best = r.max()
        tied = np.isclose(r, best, rtol=1.e-9, atol=0.)
        width = np.where(tied, hi - lo, np.inf)
        i = np.unravel_index(np.argmin(width), width.shape)
        return float(lo[i] - margin), float(hi[i] + margin)


def best_split_windows(samples, var_a, var_b, step, margin=0.):
    """Separate windows on two variables, vetoing a track if either contains it, that veto every cosmic-ray track
    with a value of either variable at the least signal loss.

    Each window on var_a (on the scan grid, or none) fixes the var_b window: the narrowest grid window containing the
    cosmic-ray tracks the var_a window misses. Returns the windows (widened by margin) and the optimum's signal loss.
    """
    sig, cos = samples[SIGNAL], samples["cosmic"]
    s_has = (sig[var_a] > NO_CRV) | (sig[var_b] > NO_CRV)
    c_has = (cos[var_a] > NO_CRV) | (cos[var_b] > NO_CRV)
    s_a, s_b, s_w = sig[var_a][s_has], sig[var_b][s_has], sig["w"][s_has]
    c_a, c_b = cos[var_a][c_has], cos[var_b][c_has]
    r_lo, r_hi = VARIABLES[var_a]["scan"]
    edges = np.arange(r_lo, r_hi + 0.5*step, step)
    candidates = [None] + [(lo, hi) for lo in edges for hi in edges if hi > lo]
    best = (np.inf, None)
    for win_a in candidates:
        in_a = np.zeros(len(c_a), bool) if win_a is None else (c_a > win_a[0]) & (c_a < win_a[1])
        if np.any(~in_a & (c_b <= NO_CRV)): continue # a track without var_b must be in the var_a window
        rest = c_b[~in_a]
        win_b = None if len(rest) == 0 else (np.floor(rest.min()/step)*step - step, np.ceil(rest.max()/step)*step + step)
        veto = np.zeros(len(s_w), bool)
        for x, win in ((s_a, win_a), (s_b, win_b)):
            if win is not None: veto |= (x > win[0]) & (x < win[1])
        loss = s_w[veto].sum() / sig["w"].sum()
        width = sum(w[1] - w[0] for w in (win_a, win_b) if w is not None)
        if (loss, width) < (best[0], best[2] if best[1] is not None else np.inf):
            best = (loss, {var_a: win_a, var_b: win_b}, width)
    windows = {v: (float(w[0] - margin), float(w[1] + margin)) for v, w in best[1].items() if w is not None}
    return windows, float(best[0])


def veto_flags(sample, windows):
    """Tracks vetoed by any of the windows {var: (lo, hi)}."""
    flag = np.zeros(len(sample["w"]), bool)
    for var, (lo, hi) in windows.items():
        x = sample[var]
        flag |= (x > lo) & (x < hi)
    return flag


def evaluate_veto(samples, flags_fn):
    """Signal/background efficiencies and S / sqrt(B) relative to no veto for a per-sample veto flag function."""
    s0 = samples[SIGNAL]["w"].sum()
    s = samples[SIGNAL]["w"][~flags_fn(samples[SIGNAL])].sum()
    b0 = sum(samples[n]["w"].sum() for n in BACKGROUNDS)
    b = sum(samples[n]["w"][~flags_fn(samples[n])].sum() for n in BACKGROUNDS)
    c0 = samples["cosmic"]["w"].sum()
    c = samples["cosmic"]["w"][~flags_fn(samples["cosmic"])].sum()
    return dict(eff_sig=s/s0, eff_bkg=b/b0, eff_cosmic=c/c0,
                rel_sig=float(significance(s, b) / significance(s0, b0)),
                rel_sig_cosmic=float(significance(s, c) / significance(s0, c0)), n_sig=s, n_bkg=b, n_cosmic=c)


#------------------------------------------------------------------------------------------------------
# Plotting helpers (as in studies/pid)
#------------------------------------------------------------------------------------------------------

def setup_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        # Times-like text and math, as with ROOT font 132 in mumep_ana/analysis (fonts bundled with matplotlib)
        "font.family": "STIXGeneral", "mathtext.fontset": "stix",
        "font.size": 14, "axes.labelsize": 16, "legend.fontsize": 12, "legend.frameon": False,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True, "figure.figsize": (6.4, 4.8),
        "axes.formatter.use_mathtext": True, "savefig.bbox": "tight",
    })
    return plt


def mu2e_label(ax, status="Simulation"):
    """Mu2e label as in mumep_ana/analysis/tools/utilities.C: bold sans "Mu2e" (ROOT font 61) followed by a
    smaller serif italic status (font 52), aligned with the left axis edge and just above the frame."""
    from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
    size = ax.xaxis.label.get_fontsize()
    mu2e = TextArea("Mu2e", textprops=dict(family=["Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
                                           weight="bold", size=1.1*size))
    text = TextArea(status, textprops=dict(style="italic", size=0.88*size))
    box = HPacker(children=[mu2e, text], align="baseline", pad=0, sep=0.3*size)
    ax.add_artist(AnchoredOffsetbox(loc="lower left", child=box, pad=0., borderpad=0., frameon=False,
                                    bbox_to_anchor=(0., 1.01), bbox_transform=ax.transAxes))


def label_figure(fig):
    """Add the Mu2e label to each axes with no other axes above it (the top row of a multi-panel figure)."""
    for ax in fig.axes:
        if ax.get_label() == "<colorbar>": continue
        p = ax.get_position()
        covered = any(o is not ax and o.get_label() != "<colorbar>" and o.get_position().y0 >= p.y1 - 1e-3
                      and o.get_position().x0 < p.x1 and o.get_position().x1 > p.x0 for o in fig.axes)
        if not covered: mu2e_label(ax)


def save(fig, name, figdir=FIG_DIR, formats=("pdf", "png")):
    """Save the figure with the Mu2e label on the top row of axes."""
    label_figure(fig)
    Path(figdir).mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(Path(figdir) / f"{name}.{ext}", dpi=150)
    print(f">>> Saved {Path(figdir) / name}.{{{','.join(formats)}}}")
