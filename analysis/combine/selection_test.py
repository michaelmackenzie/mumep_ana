#! /usr/bin/env python3
"""Estimate the sensitivity of a test selection: re-fill the dataset histograms with the selection, build a
histogram-template model and an envelope model, and fit the envelope model to the Asimov data of the template model.

    ./selection_test.py "trk_t0 > 540" t0_540 [-p mumep] [-s "45"] [-c "1000"]

1. make_selection_hists.C writes the histograms to test_hists/<tag>/ (the selection uses the selection tree branch
   names; "1" for no extra cut).
2. build_model.py builds the template model <tag>_hists and the envelope model <tag>_env (fitted signal shape), whose
   data is the template model's Asimov data set (--data-model), in models/.
3. Combine: the limit of the envelope model on that data (its observed limit) is the sensitivity estimate; the
   template model's own Asimov limit is the reference.

Outputs: test_hists/<tag>/ holds the histograms, logs/, the Combine outputs in combine/, and summary.txt.
"""

import argparse
import glob
import os
import re
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run(command, log, cwd=HERE):
    """Run a command (a string, through the shell), logging it and its output; exit on failure."""
    print(command)
    with open(log, "w") as f:
        f.write(command + "\n\n")
        f.flush()
        status = subprocess.run(command, shell=True, cwd=cwd, stdout=f, stderr=subprocess.STDOUT).returncode
    if status != 0:
        sys.exit(f"selection_test.py: command failed, see {log}")
    return open(log).read()


def limit(text, label):
    values = re.findall(rf"{re.escape(label)}:\s*r\s*<\s*([-\d.eE+]+)", text)
    return values[-1] if values else "n/a"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("selection", help="cut on top of the set selection, e.g. \"trk_t0 > 540\" (\"1\": none)")
    parser.add_argument("tag", help="test name")
    parser.add_argument("-p", "--process", default="mumep", choices=("mumem", "mumep"))
    parser.add_argument("-s", "--selections", default="45", help="selection sets, quoted (default: %(default)s)")
    parser.add_argument("-c", "--control-regions", default="", help="control region offsets, quoted; their sets "
                                                                       "are re-filled with the same selection")
    parser.add_argument("--cl", default="0.9", help="confidence level (default: %(default)s)")
    parser.add_argument("--rmax", default=None, help="upper bound of r (default: 1000 for mumep, 100 otherwise)")
    parser.add_argument("--model-args", default="", help="additional build_model.py arguments for both models")
    parser.add_argument("--fit-args", default="", help="additional Combine arguments for the limits")
    parser.add_argument("--skip-hists", action="store_true", help="use the existing test_hists/<tag>/ histograms")
    parser.add_argument("--skip-models", action="store_true", help="use the existing models")
    args = parser.parse_args()
    rmax = args.rmax or ("1000" if args.process == "mumep" else "100")

    outdir = os.path.join(HERE, "test_hists", args.tag)
    logdir, combdir = os.path.join(outdir, "logs"), os.path.join(outdir, "combine")
    os.makedirs(logdir, exist_ok=True)
    os.makedirs(combdir, exist_ok=True)
    selections = args.selections.split()
    offsets = args.control_regions.split()
    sets = [str(int(s) + int(o)) for s in selections for o in ["0"] + offsets]
    print(f">>> Selection test {args.tag}: \"{args.selection}\" on sets {' '.join(sets)}")

    # 1) Dataset histograms with the test selection
    if not args.skip_hists:
        cut = args.selection.replace('"', '\\"')
        macro = f'make_selection_hists.C("{" ".join(sets)}", "{cut}", "{outdir}")'
        text = run(f"root.exe -q -b {shlex.quote(macro)}", os.path.join(logdir, "make_selection_hists.log"))
        print("\n".join(line for line in text.splitlines() if "selected" in line))

    # 2) Template and envelope models
    name = "_".join([args.process] + selections + [f"cr{o}" for o in offsets])
    hists_dir = os.path.join("models", f"{name}_{args.tag}_hists")
    env_dir = os.path.join("models", f"{name}_{args.tag}_env")
    common = (f"{shlex.quote(sys.executable)} build_model.py -p {args.process} -s {shlex.quote(args.selections)} "
              f"--hist-path {shlex.quote(outdir)} --no-cc --exposure run1a")
    if offsets:
        common += f" -c {shlex.quote(args.control_regions)}"
    if args.model_args:
        common += f" {args.model_args}"
    if not args.skip_models:
        # the environment would switch the model types: the options below set them
        for var in ("CONVANA_HIST_PDFS", "CONVANA_USE_ENV"):
            os.environ.pop(var, None)
        print(f">>> Building the models {args.tag}_hists and {args.tag}_env (logs in {logdir})")
        run(f"{common} -t {args.tag}_hists", os.path.join(logdir, "model_hists.log"))
        run(f"{common} -t {args.tag}_env --envelope --functions --data-model {hists_dir}",
            os.path.join(logdir, "model_env.log"))
    for d in (hists_dir, env_dir):
        if not os.path.isfile(os.path.join(HERE, d, "combined.txt")):
            sys.exit(f"selection_test.py: model {d} not found")

    # 3) Limits
    options = (f"--rMin 0. --rMax {rmax} --cl {args.cl} --cminDefaultMinimizerStrategy=1 --cminApproxPreFitTolerance 0.1"
               " --cminPreScan --cminPreFit 1 --rAbsAcc 0.001 --rRelAcc 0.001 --cminRunAllDiscreteCombinations"
               " --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams"
               " --X-rtd MINIMIZER_multiMin_hideConstants")
    if args.fit_args:
        options += f" {args.fit_args}"
    results = {}
    for kind, d, extra in (("env", env_dir, ""), ("hists", hists_dir, " -t -1 --noFitAsimov")):
        with open(os.path.join(HERE, d, "commands.sh")) as f:
            t2w = next(line.strip() for line in f if line.startswith("text2workspace.py"))
        run(t2w, os.path.join(logdir, f"t2w_{kind}.log"), cwd=os.path.join(HERE, d))
        label = f"{args.tag}_{kind}"
        what = ("the envelope model on the template-model Asimov data" if kind == "env"
                else "the template model on its own Asimov data")
        print(f">>> Limit of {what}")
        results[kind] = run(f"combine -M AsymptoticLimits model.root {options}{extra} -n .{label}",
                            os.path.join(logdir, f"limit_{kind}.log"), cwd=os.path.join(HERE, d))
        for path in glob.glob(os.path.join(HERE, d, f"higgsCombine.{label}.AsymptoticLimits.*.root")):
            shutil.move(path, os.path.join(combdir, os.path.basename(path)))

    # 4) Summary
    env = results["env"]
    lines = [f"tag            : {args.tag}",
             f"selection      : {args.selection}",
             f"sets           : {' '.join(sets)}",
             f"hists model    : {hists_dir}",
             f"env model      : {env_dir}",
             f"{args.cl} CL limits on r:",
             f"  envelope fit to the MC-template Asimov : {limit(env, 'Observed Limit')}",
             f"    expected 16% / 50% / 84%             : {limit(env, 'Expected 16.0%')} / "
             f"{limit(env, 'Expected 50.0%')} / {limit(env, 'Expected 84.0%')}",
             f"  histogram model on its own Asimov      : {limit(results['hists'], 'Observed Limit')}"]
    with open(os.path.join(outdir, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f">>> Summary written to {os.path.join(outdir, 'summary.txt')}")


if __name__ == "__main__":
    main()
