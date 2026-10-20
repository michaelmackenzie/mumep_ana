#!/usr/bin/env python3
"""PID performance plots from the flat e- / mu- training samples.

Reproduces the PID v1 training selection and train/test split (MLTrain/TrkPID/TrackPIDTrain.py), evaluates
the deployed model (ArtAnalysis/TrkDiag/data/TrackPID_v1.dat), and makes:
  - pid_features[_<feature>]: normalized input feature distributions for e- and mu- (test sample)
  - pid_score:                score distributions, test sample with the training sample overlaid
  - pid_roc:                  mu- efficiency vs. e- efficiency, with the working point marked
  - pid_eff_vs_cut:           e- and mu- efficiencies vs. the score cut
  - pid_eff_vs_p:             e- and mu- efficiencies at the working point vs. track momentum
A summary of the working-point numbers is written to <out-dir>/pid_summary.json.

Example:
  python make_training_plots.py                 # full selection (~1 min), cached in outputs/
  python make_training_plots.py --wp 0.67 --figdir figures/training
"""

import argparse
import json
from pathlib import Path

import numpy as np

import pid_tools as pt

ELE = dict(color="#1f77b4", label=r"$e^-$")
MU  = dict(color="#d62728", label=r"$\mu^-$")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", default=pt.DEFAULT_DATA_DIR, help="Directory with the skimmed training ntuples")
    parser.add_argument("--signal", default=pt.DEFAULT_SIGNAL, help="Flat e- ntuple")
    parser.add_argument("--background", default=pt.DEFAULT_BKG, help="Flat mu- ntuple")
    parser.add_argument("--weights", default=None, help="SOFIE weight file (default: ArtAnalysis/TrkDiag/data/TrackPID_v1.dat)")
    parser.add_argument("--wp", type=float, default=pt.WORKING_POINT, help="Score cut working point")
    parser.add_argument("--figdir", default=str(pt.STUDY_DIR / "figures" / "training"), help="Figure output directory")
    parser.add_argument("--out-dir", default=str(pt.STUDY_DIR / "outputs"), help="Cache and summary output directory")
    parser.add_argument("--no-cache", action="store_true", help="Redo the track selection, ignoring/without writing a cache")
    parser.add_argument("--formats", default="pdf,png", help="Comma-separated figure formats")
    return parser.parse_args()


def plot_features(plt, df, figdir, formats):
    test = df[~df.is_train]
    names = ["EoverP", "dt", "logfitcon", "dtdt_slope"]

    def draw(ax, var):
        lo, hi = pt.RANGES[var]
        bins = np.linspace(lo, hi, 61)
        zero_fit = []
        for label, style in ((1, ELE), (0, MU)):
            x = test.loc[test.label == label, var]
            # normalize to all tracks, so the few entries outside the range (or undefined) are not drawn
            w = np.full(len(x), 1. / (len(x) * (bins[1] - bins[0])))
            ax.hist(x, bins=bins, weights=w, histtype="step", lw=1.5, **style)
            zero_fit.append(f"{style['label']}: {100*np.mean(np.isnan(x)):.1f}%")
        if var == "logfitcon":
            ax.text(0.03, 0.97, "Fit consistency = 0:\n" + ", ".join(zero_fit), transform=ax.transAxes,
                    va="top", ha="left", fontsize=10)
        ax.set_xlabel(pt.LABELS[var])
        ax.set_ylabel("Normalized entries")
        ax.set_xlim(lo, hi)

    # Combined 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, var in zip(axes.flat, names):
        draw(ax, var)
        ax.set_yscale("log")
    axes[0, 0].legend(loc="upper left")
    fig.tight_layout()
    pt.save(fig, figdir, "pid_features", formats)
    plt.close(fig)

    # Individual figures, linear and log
    for var in names:
        for scale in ("linear", "log"):
            fig, ax = plt.subplots()
            draw(ax, var)
            ax.set_yscale(scale)
            ax.legend(loc="best")
            pt.save(fig, figdir, f"pid_feature_{var}" + ("_log" if scale == "log" else ""), formats)
            plt.close(fig)


def plot_score(plt, df, wp, figdir, formats):
    bins = np.linspace(0., 1., 51)
    centers = 0.5*(bins[1:] + bins[:-1])
    fig, ax = plt.subplots()
    for label, style in ((1, ELE), (0, MU)):
        test  = df.loc[(df.label == label) & ~df.is_train, "prediction"]
        train = df.loc[(df.label == label) &  df.is_train, "prediction"]
        ax.hist(test, bins=bins, density=True, histtype="stepfilled", alpha=0.3, color=style["color"],
                label=style["label"] + " (test)")
        ax.hist(test, bins=bins, density=True, histtype="step", lw=1.5, color=style["color"])
        n, _ = np.histogram(train, bins=bins)
        norm = 1. / (len(train) * np.diff(bins))
        ax.errorbar(centers, n*norm, yerr=np.sqrt(n)*norm, fmt="o", ms=3, color=style["color"],
                    label=style["label"] + " (training)")
    ax.set_xlabel(pt.LABELS["prediction"])
    ax.set_ylabel("Normalized entries")
    ax.set_yscale("log")
    ax.set_xlim(0., 1.)
    ax.legend(loc="upper center", ncol=2)
    ax.set_ylim(top=ax.get_ylim()[1]*20)
    pt.cut_line(ax, wp, color="k", ls="--", lw=1)
    pt.save(fig, figdir, "pid_score", formats)
    plt.close(fig)


def plot_roc(plt, df, wp, figdir, formats):
    test = df[~df.is_train]
    s = test.loc[test.label == 1, "prediction"].to_numpy()
    b = test.loc[test.label == 0, "prediction"].to_numpy()
    cuts, eff_s, eff_b = pt.roc_curve(s, b)
    area = pt.auc(s, b)

    fig, ax = plt.subplots()
    ok = eff_b > 0
    ax.plot(eff_s[ok], eff_b[ok], color="k", lw=1.5, label=f"PID (AUC = {area:.4f})")
    ax.plot([np.mean(s > wp)], [np.mean(b > wp)], "*", color="#d62728", ms=12, label=f"Score > {wp:g}")
    ax.set_xlabel(r"$e^-$ efficiency")
    ax.set_ylabel(r"$\mu^-$ efficiency")
    ax.set_yscale("log")
    ax.set_xlim(0.8, 1.0)
    ax.set_ylim(max(0.5/len(b), 1e-4), 1.)
    ax.legend(loc="upper left")
    pt.save(fig, figdir, "pid_roc", formats)
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(cuts, eff_s, color=ELE["color"], lw=1.5, label=ELE["label"])
    ax.plot(cuts[eff_b > 0], eff_b[eff_b > 0], color=MU["color"], lw=1.5, label=MU["label"])
    ax.set_xlabel("PID score cut")
    ax.set_ylabel("Efficiency")
    ax.set_yscale("log")
    ax.set_xlim(0., 1.)
    ax.set_ylim(max(0.5/len(b), 1e-4), 1.5)
    ax.legend(loc="lower left")
    pt.cut_line(ax, wp, color="k", ls="--", lw=1)
    pt.save(fig, figdir, "pid_eff_vs_cut", formats)
    plt.close(fig)
    return area


def plot_eff_vs_p(plt, df, wp, figdir, formats):
    test = df[~df.is_train]
    bins = np.linspace(80., 120., 9)  # the flat samples end near 120 MeV/c
    centers, half = 0.5*(bins[1:] + bins[:-1]), 0.5*np.diff(bins)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6.4, 6.4), gridspec_kw=dict(hspace=0.05))
    for label, style, ax in ((1, ELE, ax1), (0, MU, ax2)):
        sel = test[test.label == label]
        n_tot, _  = np.histogram(sel["mom"], bins=bins)
        n_pass, _ = np.histogram(sel.loc[sel.prediction > wp, "mom"], bins=bins)
        eff, lo, hi = pt.efficiency(n_pass, n_tot)
        ax.errorbar(centers, eff, xerr=half, yerr=[lo, hi], fmt="o", ms=4, color=style["color"],
                    label=style["label"] + f", score > {wp:g}")
        ax.set_ylabel(style["label"] + " efficiency")
        ax.legend(loc="lower right" if label == 1 else "upper left")
    ax1.set_ylim(0.9, 1.0)
    ax2.set_yscale("log")
    ax2.set_xlabel(pt.LABELS["mom"])
    ax2.set_xlim(bins[0], bins[-1])
    pt.save(fig, figdir, "pid_eff_vs_p", formats)
    plt.close(fig)


def summarize(df, wp, area):
    def numbers(sel):
        n, k = len(sel), int((sel.prediction > wp).sum())
        eff, lo, hi = pt.efficiency(k, n)
        return dict(n=n, n_pass=k, eff=float(eff), err_low=float(lo), err_high=float(hi))

    test = df[~df.is_train]
    summary = dict(working_point=wp, auc_test=area, n_train_per_class=int(df[df.label == 1].is_train.sum()))
    for name, label in (("electron", 1), ("muon", 0)):
        summary[name] = {
            "train": numbers(df[(df.label == label) & df.is_train]),
            "test" : numbers(test[test.label == label]),
            # momentum windows near the mu- -> e- (104.97 MeV/c) and mu- -> e+ (92.32 MeV/c) endpoints
            "test_p100_110": numbers(test[(test.label == label) & (test.mom > 100.) & (test.mom < 110.)]),
            "test_p87_97"  : numbers(test[(test.label == label) & (test.mom >  87.) & (test.mom <  97.)]),
        }
    return summary


def main():
    args = parse_args()
    plt = pt.setup_style()
    formats = args.formats.split(",")
    cache = None if args.no_cache else str(Path(args.out_dir) / "pid_training_tracks.parquet")
    df = pt.load_training_sample(args.data_dir, args.signal, args.background, args.weights, cache)
    print(f">>> {int((df.label == 1).sum())} e- and {int((df.label == 0).sum())} mu- tracks, "
          f"{int((df.label == 1).sum() - (df[df.label == 1].is_train).sum())} / "
          f"{int((df.label == 0).sum() - (df[df.label == 0].is_train).sum())} in the test sample")

    plot_features(plt, df, args.figdir, formats)
    plot_score(plt, df, args.wp, args.figdir, formats)
    area = plot_roc(plt, df, args.wp, args.figdir, formats)
    plot_eff_vs_p(plt, df, args.wp, args.figdir, formats)

    summary = summarize(df, args.wp, area)
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    out = Path(args.out_dir) / "pid_summary.json"
    with open(out, "w") as f:
        json.dump(summary, f, indent=2)
    print(f">>> Test-sample AUC = {area:.4f}")
    for name in ("electron", "muon"):
        for sample in ("train", "test", "test_p100_110", "test_p87_97"):
            r = summary[name][sample]
            print(f"    {name:8s} {sample:14s}: eff(score > {args.wp:g}) = {r['eff']:.5f} "
                  f"+{r['err_high']:.5f} -{r['err_low']:.5f} ({r['n_pass']}/{r['n']})")
    print(f">>> Wrote {out}")


if __name__ == "__main__":
    main()
