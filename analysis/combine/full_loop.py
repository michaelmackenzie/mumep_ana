#! /usr/bin/env python3
"""Build the search model of a list of selection sets and optionally test it with Combine.

    ./full_loop.py -p mumep -s "40" [-c "1000 2000 4000"] [-t <tag>] [build_model.py options] [--test-combine]

1. build_model.py: the component shapes, the shape model of every set and control region (one channel each, in one
   card), and the cut-and-count model, in models/<process>_<sets>[_cr<offset>...][_<tag>]/ (see build_model.py -h for
   the model options: --envelope, --functions, --data, --2bb, --mc-stat, --exposure, ...).
2. --test-combine: text2workspace, the expected 90% CL limit (Asimov), and a fit with pre/post-fit figures, run in
   the model directory. --dry-run prints these commands instead.

Options not listed below are passed to build_model.py.
"""

import argparse
import json
import os
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def combine_commands(outdir, name, fit_args=""):
    """text2workspace, AsymptoticLimits, FitDiagnostics, and the post-fit figures, run from the model directory."""
    with open(os.path.join(outdir, "commands.sh")) as f:
        t2w = next(line.strip() for line in f if line.startswith("text2workspace.py"))
    with open(os.path.join(outdir, "model.json")) as f:
        info = json.load(f)
    lo, hi = next(iter(info["pois"].values()))
    envelope = any(p.get("kind") == "envelope" for ch in info["channels"].values() for p in ch["processes"].values())
    options = (f"--rMin {lo:g} --rMax {hi:g} --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1"
               " --cminPreScan --cminPreFit 1")
    if envelope:
        options += (" --cminRunAllDiscreteCombinations --X-rtd MINIMIZER_freezeDisassociatedParams"
                    " --X-rtd MINIMIZER_multiMin_hideConstants")
    if fit_args:
        options += f" {fit_args}"
    return [
        t2w,
        f"combine -M AsymptoticLimits model.root -t -1 --cl 0.9 {options} -n .{name}",
        f"combine -M FitDiagnostics model.root {options} --saveShapes --saveWithUncertainties -n .{name}",
        f"{shlex.quote(sys.executable)} -m mu2ecombine.plots fitDiagnostics.{name}.root",
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--skip-model", action="store_true", help="use the existing model directory")
    parser.add_argument("--test-combine", action="store_true", help="run text2workspace, limits, and a fit")
    parser.add_argument("--fit-args", default="", help="additional Combine arguments, quoted")
    parser.add_argument("--dry-run", action="store_true", help="print the Combine commands without running them")
    args, rest = parser.parse_known_args()
    if any(a in ("-h", "--help") for a in rest):
        parser.print_help()
        sys.exit(0)

    import build_model  # after the help: it loads ROOT and mu2ecombine
    model_args = build_model.parse_args(rest)
    outdir = model_args.output
    if args.skip_model:
        if not os.path.isfile(os.path.join(outdir, "combined.txt")):
            sys.exit(f"full_loop.py: --skip-model, but {outdir}/combined.txt does not exist")
        print(f"Using the existing model in {outdir}")
    else:
        build_model.run(model_args)

    if args.test_combine or args.dry_run:
        # mu2ecombine.plots runs in a subprocess: give it the module path this process uses
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join([os.path.dirname(os.path.dirname(build_model.mc.__file__)),
                                             env.get("PYTHONPATH", "")])
        for command in combine_commands(outdir, model_args.name, args.fit_args):
            print(f"(cd {outdir} && {command})")
            if args.dry_run:
                continue
            if subprocess.run(command, shell=True, cwd=outdir, env=env).returncode != 0:
                sys.exit(f"full_loop.py: command failed: {command}")


if __name__ == "__main__":
    main()
