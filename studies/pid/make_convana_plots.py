#!/usr/bin/env python3
"""PID plots from ConvAna histograms of the MDC2025au-001 ntuples (Offline PID v1 score, trkpid branch).

Uses the track histogram sets filled before any PID requirement:
  set 4: downstream e- fit, negative charge, p > 95 MeV/c  (mu- -> e- analysis)
  set 5: downstream e- fit, positive charge, p > 75 MeV/c  (mu- -> e+ analysis)
Tracks without a PID score (e.g. no calorimeter cluster) are in the underflow and count as failing; saturated
scores of exactly 1 are in the overflow and are folded into the last bin.

Makes in figures/convana/, per analysis (prefix ele_ / pos_):
  - <prefix>pid:       normalized PID score distributions per sample
  - <prefix>pid_eff:   efficiency for score > cut vs. the cut
  - ele_roc:           flat mu- vs. flat e- efficiency (no pileup and with pileup e-)
and prints the efficiencies at the working point.

The histograms are produced with Mu2eEvtAna/scripts/make_histograms.C (function cnv_ana, mode 1) for the
dataset IDs below, e.g. root -q -b 'make_histograms.C(1, "fele0b0s5r0101", 1, "cnv_ana", 10, -1)'.
The flat-sample datasets are in Mu2eEvtAna/scripts/datasets.C but not processed by default.

Example:
  python make_convana_plots.py --hist-dir /path/to/histograms
"""

import argparse

import numpy as np

import pid_tools as pt

# (dataset ID, legend label, color, line style)
ELE_SAMPLES = [
    ("fele0b0s5r0101", r"Flat $e^-$"               , "#1f77b4", "-" ),
    ("fele0b1s5r0101", r"Flat $e^-$ (pileup)"      , "#1f77b4", "--"),
    ("cele1b1s5r0101", r"$\mu^- \to e^-$ (pileup)" , "#2ca02c", "-" ),
    ("fmum0b0s5r0101", r"Flat $\mu^-$"             , "#d62728", "-" ),
    ("cry4ab1s5r0101", r"Cosmics (pileup)"         , "#7f7f7f", "-" ),
]
POS_SAMPLES = [
    ("fpos0b0s5r0101", r"Flat $e^+$"               , "#1f77b4", "-" ),
    ("fpos0b1s5r0101", r"Flat $e^+$ (pileup)"      , "#1f77b4", "--"),
    ("cpos1b1s5r0101", r"$\mu^- \to e^+$ (pileup)" , "#2ca02c", "-" ),
    ("cry4ab1s5r0101", r"Cosmics (pileup)"         , "#7f7f7f", "-" ),
]
ANALYSES = [  # (prefix, track set, samples, description)
    ("ele_", 4, ELE_SAMPLES, r"$e^-$ tracks, $p > 95$ MeV/c"),
    ("pos_", 5, POS_SAMPLES, r"$e^+$ tracks, $p > 75$ MeV/c"),
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hist-dir", default=None, help="ConvAna histogram directory (default: $CONVANA_HIST_PATH or the run1a area)")
    parser.add_argument("--mode", type=int, default=1, help="Histogramming mode in the file names")
    parser.add_argument("--wp", type=float, default=pt.WORKING_POINT, help="Score cut working point")
    parser.add_argument("--figdir", default=str(pt.STUDY_DIR / "figures" / "convana"), help="Figure output directory")
    parser.add_argument("--formats", default="pdf,png", help="Comma-separated figure formats")
    return parser.parse_args()


def score_hist(dataset, set_id, args):
    """Score histogram as (counts, edges, total) with total including the under/overflow."""
    h = pt.read_hist(dataset, set_id, "pid", args.hist_dir, args.mode)
    if h is None: return None
    counts, edges, under, over = h
    counts = counts.copy()
    counts[-1] += over  # saturated scores of exactly 1 fall in the overflow
    total = counts.sum() + under
    if total <= 0:
        print(f"!!! No entries for {dataset} in set {set_id}, skipping")
        return None
    return counts, edges, total


def eff_vs_cut(counts, edges, total):
    """Fraction of all tracks with score above each bin's lower edge (cuts at the bin edges)."""
    above = np.cumsum(counts[::-1])[::-1]
    return edges[:-1], above / total


def main():
    args = parse_args()
    plt = pt.setup_style()
    formats = args.formats.split(",")
    curves = {}

    for prefix, set_id, samples, desc in ANALYSES:
        fig_pid, ax_pid = plt.subplots()
        fig_eff, ax_eff = plt.subplots()
        print(f">>> {desc} (track set {set_id}): efficiency for PID score > {args.wp:g}")
        n_drawn = 0
        for dataset, label, color, ls in samples:
            h = score_hist(dataset, set_id, args)
            if h is None: continue
            counts, edges, total = h
            cuts, eff = eff_vs_cut(counts, edges, total)
            curves[(set_id, dataset)] = (cuts, eff)
            i_wp = np.argmin(np.abs(cuts - args.wp))
            if abs(cuts[i_wp] - args.wp) > 1e-6:
                print(f"!!! Working point {args.wp} is not a bin edge, using {cuts[i_wp]:.4f}")
            scored = counts[edges[:-1] >= 0.].sum()  # tracks with a valid score in [0, 1]
            print(f"    {dataset}: {eff[i_wp]:.5f} of all tracks, {counts[i_wp:].sum()/scored if scored > 0 else np.nan:.5f}"
                  f" of tracks with a score (N(tracks) = {total:.0f})")

            # draw the [0, 1] score range, normalized to all tracks in the set
            i0 = np.searchsorted(edges, 0.)
            ax_pid.stairs(counts[i0:] / (total*np.diff(edges[i0:])), edges[i0:],
                          color=color, ls=ls, lw=1.5, label=label)
            ax_eff.plot(cuts[i0:], eff[i0:], color=color, ls=ls, lw=1.5, label=label)
            n_drawn += 1

        if n_drawn == 0:
            print(f"!!! No histograms found for {prefix}, skipping its figures")
            plt.close(fig_pid); plt.close(fig_eff)
            continue
        for ax in (ax_pid, ax_eff):
            ax.set_xlabel(pt.LABELS["prediction"])
            ax.set_yscale("log")
            ax.set_xlim(0., 1.)
        ax_pid.set_ylabel("Normalized entries")
        ax_pid.legend(loc="upper center", title=desc, ncol=2)
        ax_pid.set_ylim(top=ax_pid.get_ylim()[1]*50)
        ax_eff.set_ylabel("Efficiency")
        ax_eff.legend(loc="lower left", title=desc)
        for ax in (ax_pid, ax_eff):
            pt.cut_line(ax, args.wp, color="k", ls=":", lw=1)
        pt.save(fig_pid, args.figdir, prefix + "pid", formats)
        pt.save(fig_eff, args.figdir, prefix + "pid_eff", formats)
        plt.close(fig_pid); plt.close(fig_eff)

    # ROC: flat mu- vs. flat e- in the e- track set
    mu = curves.get((4, "fmum0b0s5r0101"))
    if mu is None:
        print("!!! No flat mu- histograms, skipping the ROC curve")
        return
    fig, ax = plt.subplots()
    for dataset, label, color, ls in ELE_SAMPLES:
        if not dataset.startswith("fele") or (4, dataset) not in curves: continue
        cuts, eff_e = curves[(4, dataset)]
        ok = mu[1] > 0
        ax.plot(eff_e[ok], mu[1][ok], color=color, ls=ls, lw=1.5, label=label)
        i_wp = np.argmin(np.abs(cuts - args.wp))
        ax.plot([eff_e[i_wp]], [mu[1][i_wp]], "*", color="#d62728", ms=10)
    ax.plot([], [], "*", color="#d62728", ms=10, label=f"Score > {args.wp:g}")
    ax.set_xlabel(r"$e^-$ efficiency")
    ax.set_ylabel(r"Flat $\mu^-$ efficiency (no pileup)")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 1.0)
    ax.legend(loc="upper left", title=r"$p > 95$ MeV/c")
    pt.save(fig, args.figdir, "ele_roc", formats)
    plt.close(fig)


if __name__ == "__main__":
    main()
