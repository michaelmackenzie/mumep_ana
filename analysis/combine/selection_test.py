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

Options: --crv-from-tree (sets without a CRV veto or CRV-tagged copy, e.g. 81: veto from trk_crv_raw_dt, and the
inverted veto fills S+1000, which only gives build_model.py the cosmic shape), --target NAME (mu- --> e- signal and DIO of
a stopping-target design, also re-filled), --no-limits, --work-dir DIR (test_hists/ and models/ go there), --datasets.
Example, a target design with set 81 (no CRV in the tree selection), models only:

    ./selection_test.py "trk_cluster > 0. && trk_opa == 0 && trk_ntsda == 0 && <cuts>" tgt85 -p mumem -s "81" \\
        --crv-from-tree --target nominal --no-limits --work-dir <scratch>

Outputs: test_hists/<tag>/ holds the histograms, logs/, the Combine outputs in combine/, and summary.txt.
"""

import argparse
import ast
import glob
import os
import re
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def read_samples():
    """TARGETS and TARGET_HIST_PATH of samples.py (read without importing it: that needs mu2ecombine and ROOT; the
    values are literals and arithmetic)."""
    with open(os.path.join(HERE, "samples.py")) as f:
        nodes = [n for n in ast.parse(f.read()).body if isinstance(n, ast.Assign) and isinstance(n.targets[0], ast.Name)
                 and n.targets[0].id in ("TARGETS", "TARGET_HIST_PATH")]
    values = {n.targets[0].id: eval(compile(ast.Expression(n.value), "samples.py", "eval"), {"__builtins__": {}})
              for n in nodes}
    return values["TARGETS"], values["TARGET_HIST_PATH"]


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
    TARGETS, TARGET_HIST_PATH = read_samples()
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
    parser.add_argument("--crv-from-tree", action="store_true", help="apply the CRV veto from trk_crv_raw_dt to each "
                        "set and fill S+1000 with the inverted veto (for sets without a CRV veto or tagged copy)")
    parser.add_argument("--target", default=None, choices=sorted(TARGETS), help="stopping-target design (mumem): re-fill "
                        "its signal and DIO histograms too (from CONVANA_TARGET_HIST_PATH or the default) and use them")
    parser.add_argument("--no-limits", action="store_true", help="build the histograms and models, skip Combine")
    parser.add_argument("--work-dir", default=HERE, help="base directory of test_hists/ and models/ (default: here)")
    parser.add_argument("--datasets", default="", help="dataset keys to re-fill, quoted (default: all)")
    parser.add_argument("--skip-hists", action="store_true", help="use the existing test_hists/<tag>/ histograms")
    parser.add_argument("--skip-models", action="store_true", help="use the existing models")
    args = parser.parse_args()
    if args.target and args.process != "mumem":
        parser.error("--target: the stopping-target samples are mu- --> e- (mumem) samples")
    if args.target:
        args.tag = f"{args.tag}_{args.target}"
    work = os.path.abspath(args.work_dir)
    rmax = args.rmax or ("1000" if args.process == "mumep" else "100")

    outdir = os.path.join(work, "test_hists", args.tag)
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
        tokens, set_cuts = list(sets), []
        if args.crv_from_tree:
            # the CRV veto from the closest cluster time difference: S gets the veto, S+1000 the inverted veto
            tagged = "(trk_crv_raw_dt > 0. && trk_crv_raw_dt < 150.)"
            for s in sets:
                tokens.append(f"{s}>{int(s) + 1000}")
                set_cuts += [f"{s}:!{tagged}", f"{int(s) + 1000}:{tagged}"]
        set_cuts = ";".join(set_cuts).replace('"', '\\"')
        tail = f', "{set_cuts}"' if set_cuts else ""
        macro = f'make_selection_hists.C("{" ".join(tokens)}", "{cut}", "{outdir}"'
        macro += f', "{args.datasets}", true{tail})' if (args.datasets or tail) else ")"
        text = run(f"root.exe -q -b {shlex.quote(macro)}", os.path.join(logdir, "make_selection_hists.log"))
        print("\n".join(line for line in text.splitlines() if "selected" in line))
        if args.target:
            # the target's signal and DIO come from the target histogram directory, written to the same directory
            tdir = os.path.join(os.environ.get("CONVANA_TARGET_HIST_PATH") or TARGET_HIST_PATH, "")
            names = " ".join([TARGETS[args.target][0][0], TARGETS[args.target][1][0]])
            # copy_all: the other sets of these files (e.g. set 10 for the DIO tail) are needed by the model
            macro = f'make_selection_hists.C("{" ".join(tokens)}", "{cut}", "{outdir}", "{names}", true, "{set_cuts}")'
            text = run(f"CONVANA_HIST_PATH={shlex.quote(tdir)} root.exe -q -b {shlex.quote(macro)}",
                       os.path.join(logdir, "make_selection_hists_target.log"))
            print("\n".join(line for line in text.splitlines() if "selected" in line))

    # 2) Template and envelope models
    name = "_".join([args.process] + selections + [f"cr{o}" for o in offsets])
    hists_dir = os.path.join(work, "models", f"{name}_{args.tag}_hists")
    env_dir = os.path.join(work, "models", f"{name}_{args.tag}_env")
    common = (f"{shlex.quote(sys.executable)} build_model.py -p {args.process} -s {shlex.quote(args.selections)} "
              f"--hist-path {shlex.quote(outdir)} --no-cc --exposure run1a")
    if args.target:
        common = f"CONVANA_TARGET_HIST_PATH={shlex.quote(outdir)} {common} --target {args.target}"
    if offsets:
        common += f" -c {shlex.quote(args.control_regions)}"
    if args.model_args:
        common += f" {args.model_args}"
    if not args.skip_models:
        # the environment would switch the model types: the options below set them
        for var in ("CONVANA_HIST_PDFS", "CONVANA_USE_ENV"):
            os.environ.pop(var, None)
        print(f">>> Building the models {args.tag}_hists and {args.tag}_env (logs in {logdir})")
        run(f"{common} -t {args.tag}_hists -o {shlex.quote(hists_dir)}", os.path.join(logdir, "model_hists.log"))
        run(f"{common} -t {args.tag}_env -o {shlex.quote(env_dir)} --envelope --functions "
            f"--data-model {shlex.quote(hists_dir)}",
            os.path.join(logdir, "model_env.log"))
    for d in (hists_dir, env_dir):
        if not os.path.isfile(os.path.join(d, "combined.txt")):
            sys.exit(f"selection_test.py: model {d} not found")

    # 3) Limits
    if args.no_limits:
        lines = [f"tag            : {args.tag}", f"selection      : {args.selection}", f"sets           : {' '.join(sets)}",
                 f"hists model    : {hists_dir}", f"env model      : {env_dir}", "limits         : skipped (--no-limits)"]
        with open(os.path.join(outdir, "summary.txt"), "w") as f:
            f.write("\n".join(lines) + "\n")
        print("\n".join(lines))
        return
    options = (f"--rMin 0. --rMax {rmax} --cl {args.cl} --cminDefaultMinimizerStrategy=1 --cminApproxPreFitTolerance 0.1"
               " --cminPreScan --cminPreFit 1 --rAbsAcc 0.001 --rRelAcc 0.001 --cminRunAllDiscreteCombinations"
               " --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams"
               " --X-rtd MINIMIZER_multiMin_hideConstants")
    if args.fit_args:
        options += f" {args.fit_args}"
    results = {}
    for kind, d, extra in (("env", env_dir, ""), ("hists", hists_dir, " -t -1 --noFitAsimov")):
        with open(os.path.join(d, "commands.sh")) as f:
            t2w = next(line.strip() for line in f if line.startswith("text2workspace.py"))
        run(t2w, os.path.join(logdir, f"t2w_{kind}.log"), cwd=d)
        label = f"{args.tag}_{kind}"
        what = ("the envelope model on the template-model Asimov data" if kind == "env"
                else "the template model on its own Asimov data")
        print(f">>> Limit of {what}")
        results[kind] = run(f"combine -M AsymptoticLimits model.root {options}{extra} -n .{label}",
                            os.path.join(logdir, f"limit_{kind}.log"), cwd=d)
        for path in glob.glob(os.path.join(d, f"higgsCombine.{label}.AsymptoticLimits.*.root")):
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
