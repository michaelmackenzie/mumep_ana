#! /usr/bin/env python3
"""Compare the expected sensitivities of models: Asimov CLs upper limits per model.

    ./selection_comparisons.py mumem_75_r0106_hists mumep_40_r0106_env [--single-test]

A model is a name in models/, a model directory, or a card. For each model:
  1. the shape model with all uncertainties (the only limit with --single-test)
  2. the shape model with the constrained nuisances frozen (statistics only)
  3./4. the same for the cut-and-count model (cc/)
  5. a toy-based CLs limit (HybridNew) of the shape model, unless --skip-toys

The limits are Asimov (-t -1 --noFitAsimov): the observed limit on the Asimov data set is the median expected limit,
and the bands are the expected quantiles. The Combine logs and outputs go to logs/comparisons/ and the model
directories; a summary table is printed at the end.
"""

import argparse
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
QUANTILES = ["2.5", "16.0", "50.0", "84.0", "97.5"]


def model_card(name):
    """The shape card of a model name, directory, or card."""
    for path in (name, os.path.join(name, "combined.txt"), os.path.join(HERE, "models", name, "combined.txt")):
        if os.path.isfile(path):
            return os.path.abspath(path)
    sys.exit(f"selection_comparisons.py: no model '{name}' (a name in models/, a model directory, or a card)")


def reference_br(card):
    with open(card) as f:
        for line in f:
            if line.startswith("# R_mue used for signal:"):
                return float(line.split(":")[1])
    return 1.


def parse(text):
    out = {}
    m = re.findall(r"Observed Limit:\s*r\s*<\s*([-\d.eE+]+)", text)
    out["median"] = float(m[-1]) if m else None
    for q in QUANTILES:
        m = re.findall(rf"Expected\s*{re.escape(q)}%:\s*r\s*<\s*([-\d.eE+]+)", text)
        out[q] = float(m[-1]) if m else None
    m = re.findall(r"Limit:\s*r\s*<\s*([-\d.eE+]+)\s*\+/-", text)   # HybridNew
    if out["median"] is None and m:
        out["median"] = float(m[-1])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("models", nargs="+", help="model names in models/, model directories, or cards")
    parser.add_argument("--cl", default="0.9", help="confidence level of the CLs limits (default: %(default)s)")
    parser.add_argument("--single-test", "--single-fit", dest="single", action="store_true",
                        help="only the shape-model limit with all uncertainties, for each model")
    parser.add_argument("--skip-toys", action="store_true", help="skip the toy-based CLs limits")
    parser.add_argument("--fit-args", default="", help="additional Combine arguments for every fit, quoted")
    parser.add_argument("--dry-run", action="store_true", help="print the Combine commands without running them")
    args = parser.parse_args()

    cards = [model_card(m) for m in args.models]
    common = (f"-t -1 --rMin 0. --rMax {{r_max}} --cl {args.cl} --cminDefaultMinimizerStrategy=1"
              " --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rAbsAcc 0.001 --rRelAcc 0.001"
              " --cminRunAllDiscreteCombinations --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
              " --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --noFitAsimov")
    toys = (f"-M HybridNew --LHCmode LHC-limits -t -1 --rMin 0. --rMax 30. --cl {args.cl}"
            " --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1"
            " --clsAcc 0.005 --confidenceTolerance 0.05 --interpAcc 0.1 --rRelAcc 0.001 --rAbsAcc 0.01")
    extra = f" {args.fit_args}" if args.fit_args else ""
    freeze = " --freezeParameters allConstrainedNuisances"
    logdir = os.path.join(HERE, "logs", "comparisons")
    os.makedirs(logdir, exist_ok=True)
    print(f">>> Using {args.cl} CL CLs")

    rows = []
    for card in cards:
        name = os.path.basename(os.path.dirname(card)) if os.path.basename(card) == "combined.txt" else \
            os.path.splitext(os.path.basename(card))[0]
        unit = reference_br(card)
        r_max = "1000" if "mumep" in name else "100"
        cc = os.path.join(os.path.dirname(card), "cc", "combined.txt")
        tests = [("shape", card, common)]
        if not args.single:
            tests += [("shape, stat. only", card, common + freeze)]
            if os.path.isfile(cc):
                tests += [("cut-and-count", cc, common), ("cut-and-count, stat. only", cc, common + freeze)]
            if not args.skip_toys:
                tests += [("shape, toys", card, toys)]
        print("=" * 60 + f"\nModel {name} ({card})\n" + "=" * 60)
        for label, path, options in tests:
            tag = f"{name}_" + re.sub(r"[^a-z0-9]+", "_", label).strip("_")
            command = f"combine -d {os.path.basename(path)} {options.format(r_max=r_max)}{extra} -n .{tag}"
            print(f"(cd {os.path.dirname(path)} && {command})")
            if args.dry_run:
                continue
            log = os.path.join(logdir, f"{tag}.log")
            with open(log, "w") as f:
                status = subprocess.run(command, shell=True, cwd=os.path.dirname(path), stdout=f,
                                        stderr=subprocess.STDOUT).returncode
            res = parse(open(log).read())
            rows.append((name, label, unit, res, status, log))
            median = res["median"]
            print(f"  {label}: " + (f"r < {median:.4g} ({median * unit:.3e})" if median is not None
                                    else f"no limit found (status {status}), see {log}"))

    if rows:
        print("\nSummary: median expected (Asimov) limits and the expected bands, in units of r and of the branching "
              "fraction")
        print(f"{'model':28s} {'test':26s} {'median':>9s} {'-2s':>9s} {'-1s':>9s} {'+1s':>9s} {'+2s':>9s} {'B(median)':>11s}")
        for name, label, unit, res, status, log in rows:
            fmt = lambda v: f"{v:9.4g}" if v is not None else f"{'-':>9s}"
            br = f"{res['median'] * unit:11.3e}" if res["median"] is not None else f"{'-':>11s}"
            print(f"{name:28s} {label:26s} {fmt(res['median'])} {fmt(res['2.5'])} {fmt(res['16.0'])} "
                  f"{fmt(res['84.0'])} {fmt(res['97.5'])} {br}")
    print("Done.")


if __name__ == "__main__":
    main()
