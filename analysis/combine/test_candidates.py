#!/usr/bin/env python3
"""Combine model test of the selection candidates of optimize_selection_fast.C: the expected limits of the top
candidates, with models built from the selection trees.

    ./test_candidates.py <optimize dir> [-n 3]       # e.g. ../optimize/mumem_tgt_nominal

optimize_selection_fast.C writes <optimize dir>/candidates.json (the candidates ranked by the test-half figure of
merit) and runs this with its n_combine argument; run it on its own to (re)do the test, in the Combine environment
(build_model.py needs the full mu2ecombine package, the limits need combine). Each of the top N candidates, and a
baseline with no extra cut, goes through selection_test.py: the dataset histograms are re-filled from the selection
trees with the candidate's selection, the template (hists) and envelope (env) models are built, and the Combine
limits are evaluated: the template model on its own Asimov data (the expected limit), and the envelope model on the
template model's Asimov data. The models are built with build_model.py --plots (figures in <model>/figures/) for
inspection, unless --no-plots.

Outputs, in <optimize dir>:
  combine_test.txt   the candidates ranked by the expected limit, relative to the baseline, one line each by tag (c1,
                     c2, ..., baseline), then their selections and where each tag's details are
  combine_test/      the selection_test.py work directory: test_hists/<tag>/ (histograms, logs, summary.txt), models/,
                     and logs/<tag>.log
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.realpath(__file__))
SELECTION_TEST = os.path.join(HERE, "selection_test.py")


def read_summary(path):
    """Limits and model directories of a selection_test.py summary.txt: {hists, env, env_exp16, env_exp50, env_exp84,
    hists_model, env_model} (None if missing)."""
    out = dict.fromkeys(("hists", "env", "env_exp16", "env_exp50", "env_exp84", "hists_model", "env_model"))
    if not os.path.isfile(path):
        return out
    number = lambda v: None if v.strip() in ("", "n/a") else float(v)
    for line in open(path):
        if line.strip().startswith("envelope fit to the MC-template Asimov"):
            out["env"] = number(line.split(":")[-1])
        elif line.strip().startswith("expected 16% / 50% / 84%"):
            values = line.split(":")[-1].split("/")
            if len(values) == 3:
                out["env_exp16"], out["env_exp50"], out["env_exp84"] = (number(v) for v in values)
        elif line.strip().startswith("histogram model on its own Asimov"):
            out["hists"] = number(line.split(":")[-1])
        elif line.startswith("hists model"):
            out["hists_model"] = line.split(":", 1)[1].strip()
        elif line.startswith("env model"):
            out["env_model"] = line.split(":", 1)[1].strip()
    return out


def check_environment(no_limits):
    """Fail early, with a pointer to the Combine environment, if the model building or the limits can't run."""
    problems = []
    if not no_limits and shutil.which("combine") is None:
        problems.append("combine is not on PATH")
    try:
        sys.path.insert(0, HERE)
        import samples
        samples._import_mu2ecombine()
    except ImportError as e:
        problems.append(f"the mu2ecombine package can't be imported ({e})")
    if problems:
        sys.exit("test_candidates.py: " + "; ".join(problems) + ".\n  Run it in the Combine environment (pyenv rootana "
                 "+ CombinedLimit env_standalone_mu2e.sh, see the top-level README).")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", help="optimize_selection_fast.C output directory (with candidates.json)")
    parser.add_argument("-n", "--ncandidates", type=int, default=3, help="number of top candidates (default: %(default)s)")
    parser.add_argument("--no-baseline", action="store_true", help="skip the baseline (no extra cut) model")
    parser.add_argument("--no-limits", action="store_true", help="build the histograms and models, skip Combine")
    parser.add_argument("--work-dir", default=None, help="selection_test.py work directory (default: <directory>/combine_test)")
    parser.add_argument("--model-args", default="", help="additional build_model.py arguments (selection_test.py)")
    parser.add_argument("--no-plots", action="store_true", help="build the models without build_model.py --plots")
    parser.add_argument("--fit-args", default="", help="additional Combine arguments (selection_test.py)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="reuse the candidates whose selection_test.py summary exists")
    args = parser.parse_args()

    directory = os.path.abspath(args.directory)
    with open(os.path.join(directory, "candidates.json")) as f:
        config = json.load(f)
    if not config.get("supported", True):
        sys.exit("test_candidates.py: these candidates have no model equivalent (optimize_crv: the CRV-tagged set is "
                 "merged into the selection)")
    check_environment(args.no_limits)
    work = os.path.abspath(args.work_dir or os.path.join(directory, "combine_test"))
    logdir = os.path.join(work, "logs")
    os.makedirs(logdir, exist_ok=True)

    preselection, target = config["preselection"], config["target"]
    join = lambda *cuts: " && ".join(c for c in cuts if c) or "1"
    tests = [] if args.no_baseline else [{"name": "baseline", "rank": 0, "cut": "", "eff_test": 1., "sig_test": 1., "inv_limit_test": 1.,
                                          "selection": join(preselection)}]
    for rank, c in enumerate(config["candidates"][:args.ncandidates], start=1):
        tests.append(dict(c, rank=rank, selection=join(preselection, c["cut"])))

    # one thread per process (OpenBLAS in the fits)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    failed = []
    for test in tests:
        suffix = "baseline" if test["rank"] == 0 else f"c{test['rank']}"
        tag = test["tag"] = f"{config['tag']}_{suffix}"
        summary = os.path.join(work, "test_hists", f"{tag}_{target}" if target else tag, "summary.txt")
        test["summary"] = summary
        if args.skip_existing and os.path.isfile(summary):
            print(f">>> {tag}: using the existing {summary}")
        else:
            command = [sys.executable, SELECTION_TEST, test["selection"], tag, "-p", config["process"],
                       "-s", config["sets"], "--work-dir", work]
            if target:
                command += ["--target", target]
            if config.get("crv_from_tree"):
                command.append("--crv-from-tree")
            if args.no_limits:
                command.append("--no-limits")
            model_args = " ".join(a for a in (args.model_args, "" if args.no_plots else "--plots") if a)
            if model_args:
                command.append(f"--model-args={model_args}")
            if args.fit_args:
                command.append(f"--fit-args={args.fit_args}")
            if os.path.isfile(summary):
                os.remove(summary)   # a failed test must not report the limits of an earlier run
            log = os.path.join(logdir, f"{tag}.log")
            print(f">>> {tag}: {test['selection']}\n    (log: {log})", flush=True)
            with open(log, "w") as f:
                f.write(" ".join(command) + "\n\n")
                f.flush()
                status = subprocess.run(command, cwd=HERE, env=env, stdout=f, stderr=subprocess.STDOUT).returncode
            if status != 0:
                print(f"    failed (status {status}), see {log}")
                failed.append(tag)
        test.update(read_summary(summary))

    # Ranking by the expected limit (the template model on its own Asimov data)
    base = next((t["hists"] for t in tests if t["rank"] == 0), None)
    ranked = sorted(tests, key=lambda t: (t["hists"] is None, t["hists"] or 0.))
    fmt = lambda v, spec=".4g": "n/a" if v is None else format(v, spec)
    fom = "1/l(test)" if config.get("figure_of_merit") == "inv_limit" else "sig(test)"
    fom_key = "inv_limit_test" if config.get("figure_of_merit") == "inv_limit" else "sig_test"
    # One line per test, by its short tag; the selections follow, and the details are in the tag's directories
    short = lambda t: "baseline" if t["rank"] == 0 else f"c{t['rank']}"
    full = lambda tag: f"{config['tag']}_{tag}" + (f"_{target}" if target else "")
    work_rel = os.path.relpath(work, directory)
    lines = [f"# Combine model test of {config['tag']} (process {config['process']}, set {config['sets']}"
             f"{', target ' + target if target else ''}): 90% CL limits on r",
             "# expected: the template (hists) model on its own Asimov data; env: the envelope model on it (observed,"
             " expected median)",
             "# opt: the optimizer's rank (by the test-half figure of merit); /base: the expected limit relative to the"
             " baseline (no extra cut)",
             f"# details of tag T: {work_rel}/test_hists/{full('T')}/ (summary.txt, logs/), figures in "
             f"{work_rel}/models/*_{full('T')}_{{hists,env}}/figures/", "",
             f"{'rank':4s} {'tag':9s} {'opt':>4s} {'start':10s} {'eff(test)':>9s} {fom:>9s} {'expected':>9s} {'/base':>7s}"
             f" {'env':>9s} {'env exp':>9s}"]
    for i, t in enumerate(ranked, start=1):
        ratio = t["hists"] / base if (base and t["hists"] is not None) else None
        lines.append(f"{i:<4d} {short(t):9s} {('-' if t['rank'] == 0 else str(t['rank'])):>4s} "
                     f"{t.get('name', 'baseline'):10s} {fmt(t.get('eff_test'), '.4f'):>9s} "
                     f"{fmt(t.get(fom_key), '.4f'):>9s} {fmt(t['hists']):>9s} {fmt(ratio, '.3f'):>7s} "
                     f"{fmt(t['env']):>9s} {fmt(t['env_exp50']):>9s}")
    lines += ["", "# selections"]
    lines += [f"{short(t):9s} {t['selection']}" for t in sorted(tests, key=lambda t: t["rank"])]
    if failed:
        lines += ["", f"# failed: {' '.join(failed)} (see {os.path.relpath(logdir, directory)}/<tag>.log)"]
    out = os.path.join(directory, "combine_test.txt")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f">>> Written to {out}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
