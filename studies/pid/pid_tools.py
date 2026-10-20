"""Shared tools for the PID paper plots.

The PID (v1) model is a small dense network trained in MLTrain/TrkPID/TrackPIDTrain.py on flat-spectrum
e- (signal) and mu- (background) tracks. This module:
  - reproduces that training's track selection, feature definitions, and train/test split on the skimmed
    training ntuples, and evaluates the deployed model from the TMVA SOFIE weight file used by Offline
    (ArtAnalysis/TrkDiag/data/TrackPID_v1.dat). The training ntuples predate v1, so their trkpid branch
    holds the v0 score and the v1 score is evaluated here instead.
  - reads the ConvAna histograms of the MDC2025au-001 ntuples, whose trkpid branch holds the Offline v1 score
"""

import os
from pathlib import Path

import awkward as ak
import numpy as np
import pandas as pd
import uproot

# Default inputs
DEFAULT_DATA_DIR = "/exp/mu2e/data/users/mmackenz/trkpid/data/"
DEFAULT_SIGNAL   = "nts.flate.root"
DEFAULT_BKG      = "nts.flatmu.root"
DEFAULT_HIST_DIR = "/exp/mu2e/data/projects/run1a/mumep_ana/histograms/"
WORKING_POINT    = 0.574 # Default selection (ConvAna set 79): PID() > 0.574; the selection is still being optimized

STUDY_DIR = Path(__file__).resolve().parent

# Model input features (v1), in the order the network expects them
FEATURES = ["EoverP", "dt", "fitcon", "dtdt_slope"]

# Axis labels and plotting ranges for the features and other variables
LABELS = {
    "EoverP"     : r"$E_\mathrm{cluster} / p_\mathrm{track}$",
    "dt"         : r"$t_\mathrm{track} - t_\mathrm{cluster}$ [ns]",
    "fitcon"     : r"Track fit consistency",
    "logfitcon"  : r"$\log_{10}$(track fit consistency)",
    "dtdt_slope" : r"Straw hit time slope $dt_\mathrm{hit}/dt_\mathrm{track}$",
    "mom"        : r"$p_\mathrm{track}$ [MeV/c]",
    "prediction" : r"PID score",
}
RANGES = {  # chosen to contain ~99.9% of each training-sample distribution
    "EoverP"     : (0.0, 1.1),
    "dt"         : (-6., 10.),
    "fitcon"     : (0.0, 1.0),
    "logfitcon"  : (-30., 0.),
    "dtdt_slope" : (-0.5, 3.0),
    "mom"        : (80., 130.),
    "prediction" : (0.0, 1.0),
}


def area_dir():
    """Top of the Muse work area (contains ArtAnalysis/), from $MUSE_WORK_DIR or the repository location."""
    env = os.environ.get("MUSE_WORK_DIR")
    if env: return Path(env)
    return STUDY_DIR.parents[2]  # <area>/mumep_ana/studies/pid


def default_weights():
    return str(area_dir() / "ArtAnalysis" / "TrkDiag" / "data" / "TrackPID_v1.dat")


def hist_dir():
    """ConvAna histogram directory, overridden by CONVANA_HIST_PATH as in mumep_ana/analysis/defaults.C."""
    return os.environ.get("CONVANA_HIST_PATH", DEFAULT_HIST_DIR)


#------------------------------------------------------------------------------------------------------
# Model evaluation
#------------------------------------------------------------------------------------------------------

class SofieMLP:
    """Dense ReLU network with a sigmoid output, read from a TMVA SOFIE weight file.

    The file lists (name, size, values) per tensor. Each Gemm kernel is stored row-major as (n_in, n_out),
    i.e. the Keras kernel layout, so a layer is relu(x @ W + b). The layer tensor names and order follow
    the generated inference code (ArtAnalysis/TrkDiag/inc/TrackPID_v1.hxx).
    """

    LAYERS = [  # (kernel, bias) tensor names, input to output
        ("tensor_sequential1dense1CastReadVariableOp0" , "tensor_sequential1dense1BiasAddReadVariableOp0" ),
        ("tensor_sequential1dense12CastReadVariableOp0", "tensor_sequential1dense12BiasAddReadVariableOp0"),
        ("tensor_sequential1dense21CastReadVariableOp0", "tensor_sequential1dense21BiasAddReadVariableOp0"),
        ("tensor_sequential1dense31CastReadVariableOp0", "tensor_sequential1dense31AddReadVariableOp0"    ),
    ]

    def __init__(self, path):
        tensors = {}
        with open(path) as f:
            tokens = f.read().split()
        i = 0
        while i < len(tokens):
            name, size = tokens[i], int(tokens[i+1])
            tensors[name] = np.array(tokens[i+2:i+2+size], dtype=np.float32)
            i += 2 + size
        self.layers = []
        n_in = len(FEATURES)
        for kname, bname in self.LAYERS:
            bias = tensors[bname]
            self.layers.append((tensors[kname].reshape(n_in, bias.size), bias))
            n_in = bias.size

    def __call__(self, x):
        x = np.asarray(x, dtype=np.float32)
        for i, (w, b) in enumerate(self.layers):
            x = x @ w + b
            if i < len(self.layers) - 1: x = np.maximum(x, 0.)
        return (1. / (1. + np.exp(-x[:, 0]))).astype(np.float32)


#------------------------------------------------------------------------------------------------------
# Training sample
#------------------------------------------------------------------------------------------------------

def _mag(v):
    c = v["fCoordinates"]
    return np.sqrt(c["fX"]**2 + c["fY"]**2 + c["fZ"]**2)


def select_tracks(path, mc_pdg, label, step_size="200 MB"):
    """Vectorized copy of TrackPIDTrain.apply_cut: one row per (track, tracker-entrance segment) with
    an e- fit hypothesis, the requested true pdg, an active calo hit, TrkQual >= 0.2, and a downstream
    tracker-entrance segment with 80 < p < 130 MeV/c. Rows keep the event/track/segment order, matching
    the training csv files row by row."""
    branches = ["trk/trk.pdg", "trk/trk.fitcon", "trkqual/trkqual.result",
                "trkcalohit/trkcalohit.active", "trkcalohit/trkcalohit.edep", "trkcalohit/trkcalohit.dt",
                "trksegs", "trkmcsim", "trkdtdt_slope"]
    frames = []
    tree = uproot.open(path)["EventNtuple/ntuple"]
    for arr in tree.iterate(branches, step_size=step_size, library="ak"):
        pdg_mc = ak.firsts(arr["trkmcsim"]["pdg"], axis=2)
        trk_mask = ((arr["trk/trk.pdg"] == 11) & (pdg_mc == mc_pdg) & (arr["trkcalohit/trkcalohit.active"] != 0)
                    & (arr["trkqual/trkqual.result"] >= 0.2))
        trk_mask = ak.fill_none(trk_mask, False)

        segs = arr["trksegs"]
        mom  = _mag(segs["mom"])
        seg_mask = (segs["sid"] == 1) & (segs["mom"]["fCoordinates"]["fZ"] >= 0.) & (mom >= 80.) & (mom <= 130.)
        seg_mask = seg_mask & trk_mask  # broadcast the track mask over its segments

        def per_seg(x):  # broadcast a per-track quantity over the selected segments, then flatten
            return ak.to_numpy(ak.flatten(ak.broadcast_arrays(x, segs["sid"])[0][seg_mask], axis=None))

        frames.append(pd.DataFrame({
            "mom"        : ak.to_numpy(ak.flatten(mom[seg_mask], axis=None)),
            "fitcon"     : per_seg(arr["trk/trk.fitcon"]),
            "dtdt_slope" : per_seg(arr["trkdtdt_slope"]),
            "edep"       : per_seg(arr["trkcalohit/trkcalohit.edep"]),
            "dt"         : per_seg(arr["trkcalohit/trkcalohit.dt"]),
        }))
    df = pd.concat(frames, ignore_index=True)
    df["EoverP"]    = df["edep"] / df["mom"]
    df["label"]     = label
    return df


def split_train_test(n_sig, n_bkg, n_train=70000, frac_train=0.7):
    """Number of leading rows per class used in training, as in TrackPIDTrain.py (defaults match it)."""
    max_train = int(2*int(min(int(2.*frac_train*min(n_sig, n_bkg)), n_train)/2))
    return int(max_train/2)


def load_training_sample(data_dir=DEFAULT_DATA_DIR, signal=DEFAULT_SIGNAL, background=DEFAULT_BKG,
                         weights=None, cache=None, n_train=70000, frac_train=0.7):
    """Selected e- (label 1) and mu- (label 0) tracks with features, the PID v1 score ('prediction'),
    and an 'is_train' flag. With cache set, the selection is read from/written to that parquet file."""
    if cache and Path(cache).exists():
        print(f">>> Reading cached selection from {cache}")
        df = pd.read_parquet(cache)
    else:
        print(f">>> Selecting signal tracks from {os.path.join(data_dir, signal)}")
        df_sig = select_tracks(os.path.join(data_dir, signal), mc_pdg=11, label=1)
        print(f">>> Selecting background tracks from {os.path.join(data_dir, background)}")
        df_bkg = select_tracks(os.path.join(data_dir, background), mc_pdg=13, label=0)
        half = split_train_test(len(df_sig), len(df_bkg), n_train, frac_train)
        df_sig["is_train"] = np.arange(len(df_sig)) < half
        df_bkg["is_train"] = np.arange(len(df_bkg)) < half
        df = pd.concat([df_sig, df_bkg], ignore_index=True)
        if cache:
            Path(cache).parent.mkdir(parents=True, exist_ok=True)
            df.to_parquet(cache)
            print(f">>> Wrote the selection to {cache}")
    # log10 of the fit consistency, undefined (NaN) for fits with a consistency of exactly 0
    df["logfitcon"] = np.log10(df["fitcon"].where(df["fitcon"] > 0.))
    model = SofieMLP(weights or default_weights())
    df["prediction"] = model(df[FEATURES].to_numpy())
    return df


#------------------------------------------------------------------------------------------------------
# ConvAna histograms
#------------------------------------------------------------------------------------------------------

def read_hist(dataset, set_id, name, directory=None, mode=1):
    """(counts, edges, underflow, overflow) of ConvAna histogram trk_<set>/<name> for a dataset ID, or None
    if the file or histogram is missing."""
    path = Path(directory or hist_dir()) / f"ConvAna.cnv_ana.{dataset}.m{mode}.root"
    if not path.exists():
        print(f"!!! Missing histogram file {path}, skipping {dataset}")
        return None
    with uproot.open(path) as f:
        key = f"Ana/Hist/trk_{set_id}/{name}"
        try:
            h = f[key]
        except KeyError:
            print(f"!!! Missing histogram {key} in {path}, skipping {dataset}")
            return None
        flow = h.values(flow=True)
        return h.values(), h.axis().edges(), flow[0], flow[-1]


#------------------------------------------------------------------------------------------------------
# Statistics helpers
#------------------------------------------------------------------------------------------------------

def efficiency(n_pass, n_total):
    """Efficiency with a Clopper-Pearson 68% interval, as (eff, err_low, err_high) arrays. Assumes unit weights."""
    from scipy.stats import beta
    n_pass, n_total = np.asarray(n_pass, dtype=float), np.asarray(n_total, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        eff = np.where(n_total > 0, n_pass/n_total, np.nan)
        alpha = 1. - 0.6827
        lo = np.where(n_pass > 0, beta.ppf(alpha/2, n_pass, n_total - n_pass + 1), 0.)
        hi = np.where(n_pass < n_total, beta.ppf(1 - alpha/2, n_pass + 1, n_total - n_pass), 1.)
    return eff, eff - lo, hi - eff


def roc_curve(score_sig, score_bkg, n_points=2001):
    """Signal and background efficiencies for score > cut, for cuts spanning [0, 1]."""
    cuts = np.linspace(0., 1., n_points)
    s = np.sort(score_sig)
    b = np.sort(score_bkg)
    eff_s = 1. - np.searchsorted(s, cuts, side="right") / len(s)
    eff_b = 1. - np.searchsorted(b, cuts, side="right") / len(b)
    return cuts, eff_s, eff_b


def auc(score_sig, score_bkg):
    """Area under the ROC curve: P(score_sig > score_bkg), from the rank-sum statistic."""
    from scipy.stats import rankdata
    ranks = rankdata(np.concatenate([score_sig, score_bkg]))
    n_s, n_b = len(score_sig), len(score_bkg)
    return (ranks[:n_s].sum() - n_s*(n_s + 1)/2.) / (n_s*n_b)


#------------------------------------------------------------------------------------------------------
# Plotting helpers
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


def cut_line(ax, x, **kwargs):
    """Vertical line at x, stopped just below the legend if the legend spans x. Call after ax.legend()."""
    ymax = 1.
    legend = ax.get_legend()
    if legend is not None:
        ax.figure.canvas.draw()
        box = legend.get_window_extent().transformed(ax.transAxes.inverted())
        x_axes = ax.transAxes.inverted().transform(ax.transData.transform((x, 0.)))[0]
        if box.x0 <= x_axes <= box.x1: ymax = box.y0 - 0.02
    ax.axvline(x, ymax=ymax, **kwargs)


def label_figure(fig):
    """Add the Mu2e label to each axes with no other axes above it (the top row of a multi-panel figure)."""
    for ax in fig.axes:
        p = ax.get_position()
        covered = any(o is not ax and o.get_position().y0 >= p.y1 - 1e-3
                      and o.get_position().x0 < p.x1 and o.get_position().x1 > p.x0 for o in fig.axes)
        if not covered: mu2e_label(ax)


def save(fig, figdir, name, formats=("pdf", "png"), label_axes=None):
    """Save the figure, adding the Mu2e label to label_axes (default: the top row of axes)."""
    if label_axes is None: label_figure(fig)
    else:
        for ax in label_axes: mu2e_label(ax)
    Path(figdir).mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(Path(figdir) / f"{name}.{ext}", dpi=150)
    print(f">>> Saved {Path(figdir) / name}.{{{','.join(formats)}}}")
