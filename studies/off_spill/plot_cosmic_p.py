#!/usr/bin/env python3
"""Compare cosmic track momentum (trk_<set>/p_2) across cry4ab{1,0,n} samples.

Plots unit-area-normalized spectra in [pmin, pmax] with pol2 fits and 1-sigma
statistical bands, plus a ratio panel relative to b1.

Each sample's shape is compared with b1 (Mixed) by a likelihood-ratio test:
common pol2 shape vs separate shapes, normalizations profiled, Delta(-2lnL)
against chi2 with 2 dof. The p-value is shown in the legend.

Fits are Poisson-likelihood fits to the raw (rebinned) counts; the fit function
and its confidence band are then scaled by the same factor as the histogram.
"""
import argparse
import os
from array import array

import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLES = [("b1", "Mixed", ROOT.kBlack),
           ("b0", "On-spill", ROOT.kRed + 1),
           ("bn", "Off-spill", ROOT.kBlue + 1)]


def load(tag, hset, pmin, pmax, rebin, indir):
    path = os.path.join(indir, f"ConvAna.cnv_ana.cry4a{tag}s5r0101.m1.root")
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"cannot open {path}")
    hname = f"Ana/Hist/trk_{hset}/p_2"
    h = f.Get(hname)
    if not h:
        raise RuntimeError(f"{hname} not found in {path}")
    h = h.Clone(f"p_{tag}")
    h.SetDirectory(0)
    f.Close()
    h.Rebin(rebin)
    h.GetXaxis().SetRangeUser(pmin, pmax)
    return h


def fit_and_band(h, pmin, pmax, scale, color, npts=200):
    """pol2 likelihood fit to raw counts; return scaled TF1 and 1-sigma band."""
    xc = 0.5 * (pmin + pmax)  # centre the polynomial to decorrelate parameters
    form = f"[0] + [1]*(x-{xc}) + [2]*(x-{xc})^2"
    fn = ROOT.TF1(f"fit_{h.GetName()}", form, pmin, pmax)
    res = h.Fit(fn, "QLS0", "", pmin, pmax)
    xs = array("d", [pmin + (pmax - pmin) * i / (npts - 1) for i in range(npts)])
    ci = array("d", [0.0] * npts)
    res.GetConfidenceIntervals(npts, 1, 1, xs, ci, 0.682689, False)

    band = ROOT.TGraphErrors(npts)
    for i, x in enumerate(xs):
        band.SetPoint(i, x, fn.Eval(x) * scale)
        band.SetPointError(i, 0.0, ci[i] * scale)
    band.SetFillColorAlpha(color, 0.25)
    band.SetLineColor(color)

    fs = ROOT.TF1(f"fits_{h.GetName()}", form, pmin, pmax)
    for i in range(3):
        fs.SetParameter(i, fn.GetParameter(i) * scale)
    fs.SetLineColor(color)
    fs.SetLineWidth(2)
    print(f"  {h.GetName()}: entries={h.Integral(h.FindBin(pmin + 1e-6), h.FindBin(pmax - 1e-6)):.0f}"
          f"  status={int(res)}  params(unnormalized, p-pc)="
          + ", ".join(f"{fn.GetParameter(i):.3g}+-{fn.GetParError(i):.2g}" for i in range(3)))
    return fs, band


def _bin_moments(h, pmin, pmax):
    """Per-bin counts and integrals of 1, u, u^2 (u = p - centre) over the fit window."""
    xc = 0.5 * (pmin + pmax)
    rows = []
    for i in range(h.FindBin(pmin + 1e-6), h.FindBin(pmax - 1e-6) + 1):
        u1 = h.GetXaxis().GetBinLowEdge(i) - xc
        u2 = h.GetXaxis().GetBinUpEdge(i) - xc
        rows.append((h.GetBinContent(i), u2 - u1, (u2**2 - u1**2) / 2, (u2**3 - u1**3) / 3))
    return rows


def _shape_nll(rows, a, b):
    """-ln L of a binned pol2 shape (1 + a u + b u^2), normalization profiled out."""
    p = [m0 + a * m1 + b * m2 for _, m0, m1, m2 in rows]
    tot = sum(p)
    nll = 0.0
    for (n, *_), pi in zip(rows, p):
        if n <= 0:
            continue
        if pi <= 0 or tot <= 0:
            return 1e30
        nll -= n * ROOT.TMath.Log(pi / tot)
    return nll


def _minimize(nll):
    fcn = lambda x: nll(x[0], x[1])  # noqa: E731  (must outlive the Functor)
    functor = ROOT.Math.Functor(fcn, 2)
    m = ROOT.Math.Factory.CreateMinimizer("Minuit2", "Migrad")
    m.SetPrintLevel(-1)
    m.SetErrorDef(0.5)
    m.SetFunction(functor)
    m.SetVariable(0, "a", 0.0, 0.01)
    m.SetVariable(1, "b", 0.0, 0.01)
    m.Minimize()
    return m.MinValue(), m.Status()


def shape_compat(h, href, pmin, pmax):
    """Likelihood-ratio test that h and href share one pol2 shape.

    Separate shapes (4 shape params) vs a common shape (2), normalizations
    profiled. Delta(-2 ln L) ~ chi2(2 dof) asymptotically. Returns (q, p-value).
    """
    r, rref = _bin_moments(h, pmin, pmax), _bin_moments(href, pmin, pmax)
    n_sep = _minimize(lambda a, b: _shape_nll(r, a, b))[0] \
        + _minimize(lambda a, b: _shape_nll(rref, a, b))[0]
    n_joint, status = _minimize(lambda a, b: _shape_nll(r, a, b) + _shape_nll(rref, a, b))
    q = max(0.0, 2.0 * (n_joint - n_sep))
    return q, ROOT.TMath.Prob(q, 2), status


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--set", default=79, type=int, help="trk_<set> histogram set (default 79)")
    ap.add_argument("--pmin", default=100.0, type=float)
    ap.add_argument("--pmax", default=110.0, type=float)
    ap.add_argument("--rebin", default=10, type=int)
    ap.add_argument("--indir", default=HERE, help="directory with the ConvAna hist files")
    ap.add_argument("-o", "--out", default=None, help="output file (default cosmic_p_trk<set>.png)")
    args = ap.parse_args()
    out = args.out or f"cosmic_p_trk{args.set}.png"

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    hists, fits, bands, raw = {}, {}, {}, {}
    for label, tag, color in SAMPLES:
        h = load(label, args.set, args.pmin, args.pmax, args.rebin, args.indir)
        raw[label] = h.Clone(f"raw_{label}")
        raw[label].SetDirectory(0)
        b0, b1 = h.FindBin(args.pmin + 1e-6), h.FindBin(args.pmax - 1e-6)
        integ = h.Integral(b0, b1)
        if integ <= 0:
            raise RuntimeError(f"{tag}: no entries in [{args.pmin}, {args.pmax}]")
        scale = 1.0 / integ
        fits[label], bands[label] = fit_and_band(h, args.pmin, args.pmax, scale, color)
        h.Scale(scale)
        h.SetLineColor(color)
        h.SetMarkerColor(color)
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.8)
        hists[label] = h

    ref_tag = SAMPLES[0][1]
    pvals = {}
    for label, tag, _ in SAMPLES[1:]:
        q, pv, status = shape_compat(raw[label], raw["b1"], args.pmin, args.pmax)
        pvals[label] = pv
        print(f"  shape compatibility {tag} vs {ref_tag}: -2DeltalnL={q:.2f} (2 dof)"
              f"  p={pv:.3g}  (minimizer status {status})")

    c = ROOT.TCanvas("c", "c", 800, 800)
    top = ROOT.TPad("top", "", 0, 0.32, 1, 1)
    bot = ROOT.TPad("bot", "", 0, 0, 1, 0.32)
    top.SetBottomMargin(0.02)
    bot.SetTopMargin(0.03)
    bot.SetBottomMargin(0.32)
    for p in (top, bot):
        p.SetLeftMargin(0.13)
        p.Draw()

    top.cd()
    ymax = max(h.GetMaximum() + h.GetBinError(h.GetMaximumBin()) for h in hists.values())
    ref = hists["b1"]
    ref.SetMaximum(1.35 * ymax)
    ref.SetMinimum(0.0)
    ref.GetYaxis().SetTitle(f"Fraction / {ref.GetBinWidth(1):.2f} MeV/c")
    ref.GetYaxis().SetTitleSize(0.05)
    ref.GetYaxis().SetLabelSize(0.045)
    ref.GetXaxis().SetLabelSize(0)
    ref.Draw("E")
    for label, _, _ in SAMPLES:
        bands[label].Draw("3 same")
        fits[label].Draw("same")
    for label, _, _ in SAMPLES:
        hists[label].Draw("E same")
    leg = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg.SetBorderSize(0)
    for label, tag, _ in SAMPLES:
        ptxt = f", p(shape)={pvals[label]:.2g}" if label in pvals else ""
        leg.AddEntry(hists[label], f"{tag} ({hists[label].GetEntries():.0f} ent.{ptxt})", "lep")
    leg.AddEntry(bands["b1"], "pol2 fit #pm 1#sigma", "lf")
    leg.Draw()
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.04)
    lat.DrawLatex(0.16, 0.84, f"trk_{args.set}, unit area")

    # Ratio panel: hist / b1 hist (points) and fit / b1 fit (lines).
    bot.cd()
    keep = []
    frame = bot.DrawFrame(args.pmin, 0.0, args.pmax, 3.0)
    frame.GetXaxis().SetTitle("p [MeV/c]")
    frame.GetYaxis().SetTitle("Ratio to Mixed")
    for ax, s in ((frame.GetXaxis(), 0.11), (frame.GetYaxis(), 0.10)):
        ax.SetTitleSize(s)
        ax.SetLabelSize(0.09)
    frame.GetYaxis().SetTitleOffset(0.55)
    frame.GetYaxis().SetNdivisions(505)

    # b1 relative band around 1
    bref = bands["b1"]
    rb = ROOT.TGraphErrors(bref.GetN())
    for i in range(bref.GetN()):
        x, y = bref.GetPointX(i), bref.GetPointY(i)
        rb.SetPoint(i, x, 1.0)
        rb.SetPointError(i, 0.0, bref.GetErrorY(i) / y if y > 0 else 0.0)
    rb.SetFillColorAlpha(ROOT.kGray + 1, 0.35)
    rb.Draw("3 same")
    keep.append(rb)

    line = ROOT.TLine(args.pmin, 1.0, args.pmax, 1.0)
    line.SetLineStyle(2)
    line.Draw()
    for label, _, color in SAMPLES[1:]:
        fr = ROOT.TF1(f"rfit_{label}", lambda x, p, a=fits[label], b=fits["b1"]:
                      a.Eval(x[0]) / b.Eval(x[0]) if b.Eval(x[0]) != 0 else 0.0,
                      args.pmin, args.pmax, 0)
        fr.SetLineColor(color)
        fr.SetLineWidth(2)
        fr.SetNpx(200)
        fr.Draw("same")
        r = hists[label].Clone(f"ratio_{label}")
        r.Divide(hists["b1"])
        r.Draw("E same")
        keep += [fr, r]

    c.SaveAs(out)


if __name__ == "__main__":
    main()
