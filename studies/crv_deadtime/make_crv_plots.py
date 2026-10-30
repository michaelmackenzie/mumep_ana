#!/usr/bin/env python3
"""CRV veto timing window study for the mu- --> e+ selection (ConvAna set 40, 90 < p < 93 MeV/c by default).

A CRV veto window lo < t(track) - t(CRV) < hi rejects cosmic rays, but also signal tracks with an accidental
in-time CRV cluster (dead time). This compares the track - CRV time differences in the ConvAna track trees (see
crv_tools.py): the raw CRV time, the time propagated via the calorimeter or via the stopping target, and the closer of
the two propagated times (as in the ConvAna veto). It scans the veto window on each to maximize S / sqrt(B), and
jointly optimizes separate windows on the two propagation hypotheses (vetoing if either window contains the track).
The ConvAna CRV veto is removed by adding the vetoed tracks (set + 1000).

With --set 43 (ConvAna set 40 without the CRV, upstream, and multi-track vetoes; the multi-track veto is applied here
and the others emulated), it also scans the upstream (reflected) track time difference, trk_upstream_dt, and its
combination with the raw CRV time.

Makes in figures/:
  - crv_dt_<var>:      track - CRV time distributions (wide and zoomed), with the ConvAna and proposed windows
  - crv_dt_2d_<x>_<y>: one time difference vs. another for signal and cosmic-ray tracks
  - crv_scan1d_<var>:  efficiencies and S / sqrt(B) vs. one window edge, with the other at the proposal
  - crv_scan2d_<var>:  signal efficiency, background efficiency, and S / sqrt(B) vs. both window edges
and writes the windows and their efficiencies to outputs/crv_windows.json.

The proposed window maximizes S / sqrt(B) on the scan grid (the narrowest such window), widened by --margin on each
side: the cosmic-ray sample has ~100 tracks with a CRV cluster, so the optimum edges sit on the outermost of them.

Example:
  python make_crv_plots.py [--step 2.5] [--margin 10]
"""

import argparse
import json

import numpy as np

import crv_tools as ct


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--set", type=int, default=40,
                        help="ConvAna track set: 40 (with 1040, upstream veto applied) or 43 (no vetoes, adds the "
                             "upstream track time)")
    parser.add_argument("--pmin", type=float, default=90., help="Minimum track momentum [MeV/c]")
    parser.add_argument("--pmax", type=float, default=93., help="Maximum track momentum [MeV/c]")
    parser.add_argument("--step", type=float, default=2.5, help="Window edge scan step [ns]")
    parser.add_argument("--margin", type=float, default=10., help="Margin added to each edge of the optimum [ns]")
    parser.add_argument("--fom", choices=("all", "cosmic"), default="cosmic",
                        help="Background in the S / sqrt(B) that sets the proposal: all backgrounds or cosmic rays only")
    parser.add_argument("--figdir", default=None, help="Figure output directory (default: figures/, or figures/set<N> "
                                                        "for a set other than 40)")
    parser.add_argument("--formats", default="pdf,png", help="Comma-separated figure formats")
    return parser.parse_args()


def shade(ax, window, color, label, alpha=0.15):
    ax.axvspan(*window, color=color, alpha=alpha, lw=0, label=label)


def plot_distributions(plt, samples, var, windows, args):
    """Fraction of tracks per bin for the signal, cosmic rays, and the other backgrounds."""
    other = ct.merge(samples, [n for n in ct.BACKGROUNDS if n != "cosmic"])
    entries = [(samples[ct.SIGNAL], r"$\mu^- \to e^+$", ct.BLUE, "-"),
               (samples["cosmic"], "Cosmic rays", ct.COSMIC, "-"),
               (other, "RMC + RPC", ct.RMC, "--")]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8))
    for ax, (lo, hi, width) in zip(axes, ((-1000., 1000., 20.), (*ct.VARIABLES[var]["scan"], 5.))):
        bins = np.arange(lo, hi + 0.5*width, width)
        for s, label, color, ls in entries:
            x, w = s[var], s["w"]
            no_crv = w[x < ct.NO_CRV].sum() / w.sum()
            h, _ = np.histogram(x, bins=bins, weights=w / w.sum())
            ax.stairs(h, bins, color=color, ls=ls, lw=1.5, label=f"{label} ({100.*no_crv:.0f}% no CRV)")
        shade(ax, ct.VARIABLES[var]["current"], ct.GRAY, "ConvAna window")
        shade(ax, windows[var], ct.ORANGE, "Proposed window")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_xlabel(ct.VARIABLES[var]["label"])
        ax.set_ylabel(f"Fraction of tracks / {width:g} ns")
    axes[0].set_ylim(top=axes[0].get_ylim()[1]*30.)
    axes[1].set_ylim(top=axes[1].get_ylim()[1]*30.)
    axes[0].legend(loc="upper left", ncols=2, fontsize=10)
    ct.save(fig, f"crv_dt_{var}", args.figdir, args.formats)
    plt.close(fig)


def plot_2d(plt, samples, x_var, y_var, windows, args):
    """One track - CRV time difference vs. another."""
    from matplotlib.colors import LogNorm
    x_info, y_info = ct.VARIABLES[x_var], ct.VARIABLES[y_var]
    xr, yr = x_info["scan"], y_info["scan"]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), gridspec_kw=dict(wspace=0.45))
    for ax, name, title in zip(axes, (ct.SIGNAL, "cosmic"), (r"$\mu^- \to e^+$", "Cosmic rays")):
        s = samples[name]
        has = (s[x_var] > ct.NO_CRV) | (s[y_var] > ct.NO_CRV)
        if name == "cosmic": # few tracks: draw each one
            ax.scatter(s[x_var][has], s[y_var][has], s=12, color=ct.COSMIC,
                       edgecolors=np.where(s["convana_veto"][has], ct.BLACK, "none"))
            ax.scatter([], [], s=12, color=ct.COSMIC, edgecolors=ct.BLACK, label="ConvAna CRV vetoed")
        else:
            h = ax.hist2d(s[x_var][has], s[y_var][has], bins=(np.arange(*xr, 5.), np.arange(*yr, 5.)),
                          weights=s["w"][has] / s["w"].sum(), norm=LogNorm(), cmap="viridis")
            fig.colorbar(h[3], ax=ax, label="Fraction of tracks")
        ax.add_patch(plt.Rectangle((windows[x_var][0], yr[0]), np.diff(windows[x_var])[0], yr[1] - yr[0],
                                   fc="none", ec=ct.ORANGE, lw=1.5, label="Proposed windows"))
        ax.add_patch(plt.Rectangle((xr[0], windows[y_var][0]), xr[1] - xr[0], np.diff(windows[y_var])[0],
                                   fc="none", ec=ct.ORANGE, lw=1.5, ls="--"))
        ax.set_xlim(*xr)
        ax.set_ylim(*yr)
        ax.set_xlabel(x_info["label"])
        ax.set_ylabel(y_info["label"])
        ax.text(0.04, 0.93, title, transform=ax.transAxes)
    axes[1].legend(loc="upper right", fontsize=10)
    ct.save(fig, f"crv_dt_2d_{x_var}_{y_var}", args.figdir, args.formats)
    plt.close(fig)


def plot_scan_1d(plt, scan, window, args):
    """Efficiencies and S / sqrt(B) vs. each window edge, with the other edge at the proposal."""
    info = ct.VARIABLES[scan.var]
    edges = np.arange(*info["scan"], 0.5)
    fig, axes = plt.subplots(4, 2, figsize=(12.8, 12.), sharex="col", gridspec_kw=dict(hspace=0.05, wspace=0.25))
    for col, (moving, label) in enumerate((("lo", "Veto window lower edge [ns]"), ("hi", "Veto window upper edge [ns]"))):
        lo = edges if moving == "lo" else np.full_like(edges, window[0])
        hi = edges if moving == "hi" else np.full_like(edges, window[1])
        res = scan.evaluate(lo, hi)
        ax_b, ax_s, ax_r, ax_c = axes[:, col]
        ax_b.plot(edges, res["eff_cosmic"], color=ct.COSMIC, lw=1.5, label="Cosmic rays")
        ax_b.plot(edges, res["eff_bkg"], color=ct.BLACK, lw=1.5, ls="--", label="All backgrounds")
        ax_b.set_ylabel("Bkg. efficiency")
        ax_s.plot(edges, 100.*(1. - res["eff_sig"]), color=ct.BLUE, lw=1.5)
        ax_s.set_ylabel("Signal loss [%]")
        ax_s.set_ylim(bottom=0.)
        for ax, tag, ylabel in ((ax_r, "", r"$S/\sqrt{B}$ / no veto"),
                                (ax_c, "_cosmic", r"$S/\sqrt{B_\mathrm{cosmic}}$ / no veto")):
            r, e = res[f"rel_sig{tag}"], res[f"rel_sig{tag}_err"]
            ax.fill_between(edges, r - e, r + e, color=ct.BLACK, alpha=0.15, lw=0, label="MC stat.")
            ax.plot(edges, r, color=ct.BLACK, lw=1.5)
            ax.set_ylabel(ylabel)
        ax_c.set_xlabel(label)
        fixed = window[1] if moving == "lo" else window[0]
        ax_b.text(0.04, 0.08, f"{'Upper' if moving == 'lo' else 'Lower'} edge at {fixed:g} ns", transform=ax_b.transAxes)
        for ax in axes[:, col]:
            ax.axvline(window[0 if moving == "lo" else 1], color=ct.ORANGE, lw=1.2, ls=":")
            ax.axvline(info["current"][0 if moving == "lo" else 1], color=ct.GRAY, lw=1.2, ls=":")
    axes[0, 0].plot([], [], color=ct.ORANGE, ls=":", label="Proposed")
    axes[0, 0].plot([], [], color=ct.GRAY, ls=":", label="ConvAna")
    axes[0, 0].legend(loc="center right", fontsize=10)
    axes[2, 0].legend(loc="center right", fontsize=10)
    fig.suptitle(info["label"].replace(" [ns]", ""), y=0.95, fontsize=14)
    ct.save(fig, f"crv_scan1d_{scan.var}", args.figdir, args.formats)
    plt.close(fig)


def plot_scan_2d(plt, scan, window, args):
    """Signal efficiency, background efficiency, and S / sqrt(B) vs. both window edges."""
    info = ct.VARIABLES[scan.var]
    edges, lo, hi, res = scan.grid(args.step)
    valid = hi > lo
    panels = [(100.*(1. - res["eff_sig"]), "Signal loss [%]", "viridis"),
              (res["eff_bkg"], "Background efficiency", "viridis"),
              (res["rel_sig"], r"$S/\sqrt{B}$ / no veto", "viridis"),
              (res["rel_sig_cosmic"], r"$S/\sqrt{B_\mathrm{cosmic}}$ / no veto", "viridis")]
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 11.), gridspec_kw=dict(wspace=0.35, hspace=0.3))
    axes = axes.flatten()
    mesh = np.append(edges - 0.5*args.step, edges[-1] + 0.5*args.step)
    for ax, (z, label, cmap) in zip(axes, panels):
        z = np.where(valid, z, np.nan)
        pcm = ax.pcolormesh(mesh, mesh, z, cmap=cmap, shading="flat", rasterized=True)
        fig.colorbar(pcm, ax=ax, label=label)
        ax.plot(*window, marker="*", ms=14, color=ct.ORANGE, mec=ct.BLACK, ls="none", label="Proposed")
        ax.plot(*info["current"], marker="o", ms=8, color="white", mec=ct.BLACK, ls="none", label="ConvAna")
        ax.set_xlabel("Veto window lower edge [ns]")
        ax.set_ylabel("Veto window upper edge [ns]")
        ax.set_aspect("equal")
    axes[2].legend(loc="lower right", fontsize=10)
    fig.suptitle(info["label"].replace(" [ns]", ""), y=0.95, fontsize=14)
    ct.save(fig, f"crv_scan2d_{scan.var}", args.figdir, args.formats)
    plt.close(fig)


def main():
    args = parse_args()
    args.formats = args.formats.split(",")
    tag = "" if args.set == 40 else f"_set{args.set}"
    if args.figdir is None: args.figdir = str(ct.FIG_DIR / tag.lstrip("_"))
    plt = ct.setup_style()
    samples = ct.load_samples(track_set=args.set, pmin=args.pmin, pmax=args.pmax)
    no_veto = args.set == ct.NO_VETO_SET # upstream veto not applied: study it too
    variables = [v for v in ct.VARIABLES if no_veto or v != ct.UPSTREAM]
    pairs = ct.PAIRS + ([("trk_crv_raw_dt", ct.UPSTREAM)] if no_veto else [])

    if no_veto:
        # Check the emulated ConvAna CRV and upstream vetoes: the tracks passing them should match set 40 (pass) and
        # set 1040 (fail the CRV veto). The emulation uses the closest cluster or upstream track, which need not be
        # the one ConvAna vetoes on.
        print("\nEmulated ConvAna vetoes on set 43 vs. sets 40 and 1040 (weighted yields):")
        for name in [ct.SIGNAL, "cosmic"]:
            s, ref = samples[name], ct.load_sample(name, track_set=40, pmin=args.pmin, pmax=args.pmax)
            print(f"  {name:>8}: pass {s['w'][~s['convana_veto']].sum():.5g} vs. set 40 "
                  f"{ref['w'][~ref['convana_crv_veto']].sum():.5g}; CRV-vetoed, no upstream veto "
                  f"{s['w'][s['convana_crv_veto'] & ~s['convana_upstream_veto']].sum():.5g} vs. set 1040 "
                  f"{ref['w'][ref['convana_crv_veto']].sum():.5g}")
    else:
        # Check the emulated ConvAna veto (trk_crv_dt OR trk_crv_raw_dt window) against the ConvAna veto decision.
        # ConvAna applies the raw window to the matched cluster, which need not be the cluster closest in raw time.
        print("\nConvAna CRV veto vs. the emulation from the stored time differences (weighted fractions):")
        for name in [ct.SIGNAL, "cosmic"]:
            s = samples[name]
            emul = ct.veto_flags(s, ct.CONVANA_VETO)
            w = s["w"] / s["w"].sum()
            print(f"  {name:>8}: ConvAna {w[s['convana_veto']].sum():.5f}, emulated {w[emul].sum():.5f}, "
                  f"disagree {w[emul != s['convana_veto']].sum():.5f}")

    # Proposed window per variable, then the plots
    scans = {v: ct.WindowScan(samples, v) for v in variables}
    fom = "rel_sig" if args.fom == "all" else "rel_sig_cosmic"
    optimum = {v: scans[v].best_window(args.step, fom=fom) for v in variables}
    proposed = {v: scans[v].best_window(args.step, args.margin, fom=fom) for v in variables}
    alt_fom = "rel_sig_cosmic" if args.fom == "all" else "rel_sig"
    alt_optimum = {v: scans[v].best_window(args.step, fom=alt_fom) for v in variables}
    for var in variables:
        plot_distributions(plt, samples, var, proposed, args)
        plot_scan_1d(plt, scans[var], proposed[var], args)
        plot_scan_2d(plt, scans[var], proposed[var], args)
    for x_var, y_var in pairs:
        plot_2d(plt, samples, x_var, y_var, proposed, args)

    # Separate windows on two variables, vetoing if either contains the track, that veto all the tagged cosmic rays at
    # the least signal loss: the two propagation hypotheses, and the raw CRV time with the upstream track time
    joint = {}
    split_pairs = [("trk_crv_cal_dt", "trk_crv_st_dt"), ("trk_crv_st_dt", "trk_crv_cal_dt")]
    if no_veto: split_pairs += [("trk_crv_raw_dt", ct.UPSTREAM), (ct.UPSTREAM, "trk_crv_raw_dt")]
    for var_a, var_b in split_pairs:
        for margin in (0., args.margin):
            windows, _ = ct.best_split_windows(samples, var_a, var_b, args.step, margin)
            joint[f"{var_a} first, margin {margin:g}"] = windows

    # Summary of the veto options
    options = {"No veto": lambda s: np.zeros(len(s["w"]), bool),
               "ConvAna veto" + (" (CRV OR upstream)" if no_veto else ""): lambda s: s["convana_veto"]}
    if no_veto:
        options["ConvAna CRV veto only"] = lambda s: s["convana_crv_veto"]
        options["ConvAna upstream veto only"] = lambda s: s["convana_upstream_veto"]
        options[f"Raw proposed OR ConvAna upstream"] = lambda s: ct.veto_flags(
            s, {"trk_crv_raw_dt": proposed["trk_crv_raw_dt"], **ct.CONVANA_UPSTREAM_VETO})
        options[f"Raw proposed OR upstream proposed"] = lambda s: ct.veto_flags(
            s, {v: proposed[v] for v in ("trk_crv_raw_dt", ct.UPSTREAM)})
    for var in variables:
        options[f"{var} optimum {optimum[var]}"]   = lambda s, v=var: ct.veto_flags(s, {v: optimum[v]})
        options[f"{var} proposed {proposed[var]}"] = lambda s, v=var: ct.veto_flags(s, {v: proposed[v]})
        options[f"{var} {alt_fom} optimum {alt_optimum[var]}"] = lambda s, v=var: ct.veto_flags(s, {v: alt_optimum[v]})
    for label, windows in joint.items():
        text = ", ".join(f"{v.replace('trk_crv_', '').replace('trk_', '')} {w}" for v, w in windows.items())
        options[f"Split ({label.split(',')[1].strip()}): {text}"] = lambda s, ws=windows: ct.veto_flags(s, ws)
    summary = {"set": args.set, "p_window": [args.pmin, args.pmax], "step": args.step, "margin": args.margin,
               "fom": args.fom, "optimum": optimum, "proposed": proposed, "joint_or": joint, "options": {}}
    print(f"\n{'Veto':<78} {'eff(sig)':>9} {'eff(cosmic)':>11} {'eff(bkg)':>9} {'S/sqrt(B)':>10} {'S/sqrt(Bcos)':>12}")
    for label, fn in options.items():
        res = ct.evaluate_veto(samples, fn)
        summary["options"][label] = res
        print(f"{label:<78} {res['eff_sig']:9.5f} {res['eff_cosmic']:11.4f} {res['eff_bkg']:9.4f} {res['rel_sig']:10.4f}"
              f" {res['rel_sig_cosmic']:12.4f}")
    n_cosmic_crv = int((samples["cosmic"]["trk_crv_dt"] > ct.NO_CRV).sum())
    print(f"\nCosmic-ray MC tracks: {len(samples['cosmic']['w'])}, with a CRV cluster: {n_cosmic_crv}")
    summary["n_cosmic_mc"] = int(len(samples["cosmic"]["w"]))
    summary["n_cosmic_mc_with_crv"] = n_cosmic_crv
    ct.OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = ct.OUT_DIR / f"crv_windows{tag}.json"
    with open(out, "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f">>> Wrote {out}")


if __name__ == "__main__":
    main()
