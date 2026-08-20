#!/usr/bin/env python3
"""Compare two TH1 inputs from different code bases.

Usage example:
  python compare_inputs.py \
	--file1 comp_a.root --hist1 workflow/raw/signal --rate-hist1 rates/signal \
	--file2 comp_b.root --hist2 workflow/raw/signal --rate2 2.29 \
	--output compare_signal.png
"""

from __future__ import annotations

import argparse
import math
import sys

try:
	import ROOT  # type: ignore[import-not-found]
except ImportError as exc:
	print(f"ERROR: Failed to import PyROOT (ROOT): {exc}")
	raise


def _open_hist(root_path: str, hist_path: str) -> ROOT.TH1:
	root_file = ROOT.TFile.Open(root_path, "READ")
	if not root_file or root_file.IsZombie():
		raise RuntimeError(f"Unable to open ROOT file: {root_path}")

	hist = root_file.Get(hist_path)
	if not hist:
		raise RuntimeError(f"Unable to retrieve histogram '{hist_path}' from '{root_path}'")
	if not hist.InheritsFrom("TH1"):
		raise RuntimeError(f"Object '{hist_path}' in '{root_path}' is not a TH1")

	hist_copy = hist.Clone(f"{hist.GetName()}_copy")
	hist_copy.SetDirectory(0)
	root_file.Close()
	return hist_copy


def _read_rate_from_hist(root_path: str, hist_path: str) -> float:
	root_file = ROOT.TFile.Open(root_path, "READ")
	if not root_file or root_file.IsZombie():
		raise RuntimeError(f"Unable to open ROOT file: {root_path}")

	rate_obj = root_file.Get(hist_path)
	if not rate_obj:
		raise RuntimeError(f"Unable to retrieve rate histogram '{hist_path}' from '{root_path}'")
	if not rate_obj.InheritsFrom("TH1"):
		raise RuntimeError(f"Object '{hist_path}' in '{root_path}' is not a TH1")

	rate_hist = rate_obj
	if rate_hist.GetNbinsX() == 1:
		rate = float(rate_hist.GetBinContent(1))
	else:
		rate = float(rate_hist.Integral())

	root_file.Close()
	return rate


def _resolve_rate(
	input_name: str,
	root_path: str,
	rate_value: float | None,
	rate_hist_path: str | None,
) -> float:
	if rate_value is not None:
		# print(f"Using explicit expected rate for {input_name}: {rate_value:.8g}")
		return float(rate_value)

	if rate_hist_path:
		rate = _read_rate_from_hist(root_path, rate_hist_path)
		# print(f"Using expected rate from histogram for {input_name}: {rate:.8g} ({rate_hist_path})")
		return rate

	raise RuntimeError(f"No expected rate source provided for {input_name}")


def _can_rebin_to_match(h_fine: ROOT.TH1, h_coarse: ROOT.TH1, rel_tol: float = 1e-6) -> int:
	fine_bw = h_fine.GetBinWidth(1)
	coarse_bw = h_coarse.GetBinWidth(1)
	if fine_bw <= 0.0 or coarse_bw <= 0.0:
		return 0

	ratio = coarse_bw / fine_bw
	factor = int(round(ratio))
	if factor < 1:
		return 0

	if abs(ratio - factor) > rel_tol:
		return 0

	x_min_fine = h_fine.GetXaxis().GetBinLowEdge(1)
	x_min_coarse = h_coarse.GetXaxis().GetBinLowEdge(1)
	if abs(x_min_fine - x_min_coarse) > rel_tol * max(1.0, abs(x_min_coarse)):
		return 0

	if h_fine.GetNbinsX() % factor != 0:
		return 0

	return factor


def _try_match_binning(h1: ROOT.TH1, h2: ROOT.TH1) -> tuple[ROOT.TH1, ROOT.TH1, str]:
	h1_out = h1.Clone(f"{h1.GetName()}_binmatch")
	h2_out = h2.Clone(f"{h2.GetName()}_binmatch")
	h1_out.SetDirectory(0)
	h2_out.SetDirectory(0)

	bw1 = h1_out.GetBinWidth(1)
	bw2 = h2_out.GetBinWidth(1)
	if bw1 <= 0.0 or bw2 <= 0.0:
		return h1_out, h2_out, "Warning: invalid bin width, skipping rebinning"

	if math.isclose(bw1, bw2, rel_tol=1e-9, abs_tol=1e-12):
		return h1_out, h2_out, "Binning already matched"

	if bw1 < bw2:
		factor = _can_rebin_to_match(h1_out, h2_out)
		if factor > 1:
			h1_out.Rebin(factor)
			return h1_out, h2_out, f"Rebinned input 1 by factor {factor} to match input 2"
	else:
		factor = _can_rebin_to_match(h2_out, h1_out)
		if factor > 1:
			h2_out.Rebin(factor)
			return h1_out, h2_out, f"Rebinned input 2 by factor {factor} to match input 1"

	return h1_out, h2_out, "Warning: could not match binning exactly; plotting with original bins"


def _fractional_difference(a: float, b: float) -> float:
	denom = max(abs(a), abs(b), 1e-30)
	return abs(a - b) / denom


def _print_mismatch_warnings(int1: float, int2: float, rate1: float, rate2: float, threshold: float) -> None:
	int_diff = _fractional_difference(int1, int2)
	rate_diff = _fractional_difference(rate1, rate2)

	print(f"Raw integral input 1: {int1:.8g}")
	print(f"Raw integral input 2: {int2:.8g}")
	print(f"Integral fractional difference: {100.0 * int_diff:.4f}%")
	if int_diff > threshold:
		print(f"WARNING: raw integrals differ by more than {100.0 * threshold:.2f}%")

	print(f"Expected rate input 1: {rate1:.8g}")
	print(f"Expected rate input 2: {rate2:.8g}")
	print(f"Rate fractional difference: {100.0 * rate_diff:.4f}%")
	if rate_diff > threshold:
		print(f"WARNING: expected rates differ by more than {100.0 * threshold:.2f}%")


def _make_plot(
	h1: ROOT.TH1,
	h2: ROOT.TH1,
	label1: str,
	label2: str,
	out_path: str,
	title: str,
) -> None:
	ROOT.gROOT.SetBatch(True)
	ROOT.gStyle.SetOptStat(0)

	c = ROOT.TCanvas("c_compare", "c_compare", 1000, 900)
	c.cd()

	pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.30, 1.0, 1.0)
	pad_bot = ROOT.TPad("pad_bot", "pad_bot", 0.0, 0.0, 1.0, 0.30)
	pad_top.SetBottomMargin(0.02)
	pad_bot.SetTopMargin(0.02)
	pad_bot.SetBottomMargin(0.32)
	pad_top.Draw()
	pad_bot.Draw()

	h1_draw = h1.Clone("h1_draw")
	h2_draw = h2.Clone("h2_draw")
	h1_draw.SetDirectory(0)
	h2_draw.SetDirectory(0)

	h1_draw.SetLineColor(ROOT.kBlue + 1)
	h1_draw.SetMarkerColor(ROOT.kBlue + 1)
	h1_draw.SetMarkerStyle(20)
	h1_draw.SetLineWidth(2)

	h2_draw.SetLineColor(ROOT.kRed + 1)
	h2_draw.SetMarkerColor(ROOT.kRed + 1)
	h2_draw.SetMarkerStyle(24)
	h2_draw.SetLineWidth(2)

	pad_top.cd()
	y_max = max(h1_draw.GetMaximum(), h2_draw.GetMaximum(), 1e-12)
	h1_draw.SetTitle(title)
	h1_draw.GetYaxis().SetTitle("Events")
	h1_draw.GetYaxis().SetRangeUser(0.0, 1.25 * y_max)
	h1_draw.GetXaxis().SetLabelSize(0.0)
	h1_draw.Draw("E1")
	h2_draw.Draw("E1 SAME")

	leg = ROOT.TLegend(0.60, 0.72, 0.88, 0.88)
	leg.SetBorderSize(0)
	leg.SetFillStyle(0)
	leg.AddEntry(h1_draw, label1, "lep")
	leg.AddEntry(h2_draw, label2, "lep")
	leg.Draw()

	ratio = h1_draw.Clone("ratio")
	ratio.SetDirectory(0)
	ratio.Divide(h2_draw)

	pad_bot.cd()
	ratio.SetTitle("")
	ratio.GetYaxis().SetTitle("1 / 2")
	ratio.GetYaxis().SetNdivisions(505)
	ratio.GetYaxis().SetTitleSize(0.10)
	ratio.GetYaxis().SetLabelSize(0.08)
	ratio.GetYaxis().SetTitleOffset(0.45)
	ratio.GetXaxis().SetTitle(h1_draw.GetXaxis().GetTitle())
	ratio.GetXaxis().SetTitleSize(0.12)
	ratio.GetXaxis().SetLabelSize(0.10)
	ratio.GetXaxis().SetTitleOffset(1.0)
	ratio.SetLineColor(ROOT.kBlack)
	ratio.SetMarkerColor(ROOT.kBlack)
	ratio.SetMarkerStyle(20)
	ratio.GetYaxis().SetRangeUser(0.5, 1.5)
	ratio.Draw("E1")

	x_min = ratio.GetXaxis().GetXmin()
	x_max = ratio.GetXaxis().GetXmax()
	line = ROOT.TLine(x_min, 1.0, x_max, 1.0)
	line.SetLineStyle(2)
	line.Draw()

	c.SaveAs(out_path)


def parse_args(argv: list[str]) -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Compare two TH1 histograms from ROOT files")
	parser.add_argument("--file1", required=True, help="ROOT file for input 1")
	parser.add_argument("--hist1", required=True, help="TH1 path in file1")
	parser.add_argument("--file2", required=True, help="ROOT file for input 2")
	parser.add_argument("--hist2", required=True, help="TH1 path in file2")

	group1 = parser.add_mutually_exclusive_group(required=True)
	group1.add_argument("--rate1", type=float, help="Explicit expected overall rate for input 1")
	group1.add_argument("--rate-hist1", help="TH1 path in file1 to derive expected overall rate")

	group2 = parser.add_mutually_exclusive_group(required=True)
	group2.add_argument("--rate2", type=float, help="Explicit expected overall rate for input 2")
	group2.add_argument("--rate-hist2", help="TH1 path in file2 to derive expected overall rate")

	parser.add_argument("--label1", default="Input 1", help="Legend label for input 1")
	parser.add_argument("--label2", default="Input 2", help="Legend label for input 2")
	parser.add_argument("--title", default="Histogram Comparison", help="Plot title")
	parser.add_argument("--output", default="comparison.png", help="Output plot path")
	parser.add_argument(
		"--warn-threshold",
		type=float,
		default=0.01,
		help="Fractional warning threshold (default: 0.01 for 1%%)",
	)
	return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
	args = parse_args(sys.argv[1:] if argv is None else argv)

	try:
		h1 = _open_hist(args.file1, args.hist1)
		h2 = _open_hist(args.file2, args.hist2)
		rate1 = _resolve_rate("input 1", args.file1, args.rate1, args.rate_hist1)
		rate2 = _resolve_rate("input 2", args.file2, args.rate2, args.rate_hist2)
	except RuntimeError as exc:
		print(f"ERROR: {exc}")
		return 2

	h1_cmp, h2_cmp, rebin_message = _try_match_binning(h1, h2)
	if 'already' not in rebin_message: print(rebin_message)

	int1 = h1_cmp.Integral()
	int2 = h2_cmp.Integral()
	_print_mismatch_warnings(int1, int2, rate1, rate2, args.warn_threshold)

	_make_plot(h1_cmp, h2_cmp, args.label1, args.label2, args.output, args.title)
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
