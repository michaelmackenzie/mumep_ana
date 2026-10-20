"""Shared tools for the RMC plots of the mu- --> e- sensitivity paper: paths, the analysis/combine modules, and the
plotting style (as in studies/pid)."""

import sys
from pathlib import Path

STUDY_DIR   = Path(__file__).resolve().parent
REPO_DIR    = STUDY_DIR.parents[1]                     # mumep_ana/studies/eminus_rmc
COMBINE_DIR = REPO_DIR / "analysis" / "combine"
FIG_DIR     = STUDY_DIR / "figures"
DEFAULT_STUDY = COMBINE_DIR / "rmc_studies" / "mumem_75"  # rmc_constraint_study.py / rmc_sr_prescription.py outputs

# Colors: the fixed order of the analysis RMC plots (validated as a categorical set)
BLACK, BLUE, ORANGE, GREEN, PINK, GRAY = "#1F1F1F", "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#94a4a2"
# Process colors of the nominal analysis plots (analysis/datasets.C set_style): RMC external/internal kRed-7/kRed-9,
# cosmic rays kAzure-4, RPC (combined) kGreen-6
RMC_EXT, RMC_INT, COSMIC, RPC = "#ff6666", "#ff9999", "#6699ff", "#66cc66"


def combine_module(name):
    """Import a module of analysis/combine (rmc_constraint_study, rmc_sr_prescription)."""
    if str(COMBINE_DIR) not in sys.path:
        sys.path.insert(0, str(COMBINE_DIR))
    return __import__(name)


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
    """Mu2e label as in mumep_ana/analysis/tools/utilities.C: bold sans "Mu2e" followed by a smaller italic status,
    aligned with the left axis edge and just above the frame."""
    from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
    size = ax.xaxis.label.get_fontsize()
    mu2e = TextArea("Mu2e", textprops=dict(family=["Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
                                           weight="bold", size=1.1*size))
    text = TextArea(status, textprops=dict(style="italic", size=0.88*size))
    box = HPacker(children=[mu2e, text], align="baseline", pad=0, sep=0.3*size)
    ax.add_artist(AnchoredOffsetbox(loc="lower left", child=box, pad=0., borderpad=0., frameon=False,
                                    bbox_to_anchor=(0., 1.01), bbox_transform=ax.transAxes))


def chi2_levels(ax, levels=((1., "68% CL"), (3.84, "95% CL")), x=None):
    """Horizontal -2 Delta ln L levels, labeled at the left edge (or at x)."""
    for level, text in levels:
        ax.axhline(level, color="#5C5C5C", linestyle=":", linewidth=1)
        ax.text(ax.get_xlim()[0] if x is None else x, level, f" {text}", va="bottom", fontsize=10, color="#5C5C5C")


def save(fig, name, label_axes=(), status="Simulation", formats=("pdf", "png")):
    """Save the figure in figures/, adding the Mu2e label to label_axes."""
    for ax in label_axes:
        mu2e_label(ax, status)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=150)
    print(f">>> Saved {FIG_DIR / name}.{{{','.join(formats)}}}")
