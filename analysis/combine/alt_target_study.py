#! /usr/bin/env python3
"""Stopping-target design comparison: the mu- --> e- sensitivity of each design in samples.py TARGETS, for a set of
selections, with the selection-tree route of selection_test.py, and the change relative to the nominal target.

    ./alt_target_study.py                         # all tests, all targets, outputs in alt_target_tests/
    ./alt_target_study.py --tests "s89 cat_cos"   # a subset; --targets, --jobs, --work-dir, --force, --table-only
    ./alt_target_study.py --cos-splits "0.60 0.64 0.72 0.76" --d0-splits "50 75 125"   # scan the split points

Each (test, target) is one selection_test.py run (--target): the target's signal and DIO, the nominal other
backgrounds, a template model (its own Asimov limit) and an envelope model fit to that Asimov data. Runs whose
summary.txt exists are reused, so rerunning the script only fills in what is missing (--force redoes them).

Tests (see TESTS):
  ref75      set 75 as made by ConvAna (with its ST intersection/boundary cuts), no extra cut
  s89        set 81 with the common alternate-target selection (ConvAna set 89, median of the per-target optima) and the
             CRV veto from the tree (closest cluster); one category
  s89_d0open s89 without its upper d0 cut (d0 > -75 only): more signal of the alternate targets (forward, large +d0
             tracks lose the least energy in the foils) and more cosmics
  cat_cos    s89 in two categories, cos(theta) > 0.68 and <= 0.68, fit together
  cat_d0     s89_d0open in two categories, d0 > 100 mm and <= 100 mm, fit together

The categories ask whether separating the cosmic-rich region keeps the signal that a looser selection adds without the
cosmics it brings diluting the clean region.

Outputs in the work directory: test_hists/<test>_<target>/ and models/ (from selection_test.py), logs/<test>_<target>.log,
and the tables results.csv and results.txt (limits on r, at the 1e-15 reference, and the change in %).
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from selection_test import read_samples  # noqa: E402  (reads samples.py without ROOT)

# Set 81 preselection used by the per-target optimization (optimize_selection_fast.C) and the common set 89 cuts
PRE = "trk_cluster > 0. && trk_opa == 0 && trk_ntsda == 0 && trk_p > 100. && trk_p < 110."
S89 = ("trk_cos > 0.51 && trk_cos < 0.9 && trk_t0 > 530. && trk_t0 < 1585. && trk_rmax > 400. && trk_rmax < 637. && "
       "trk_ep > 0. && trk_ep < 1.2 && trk_dt > -4. && trk_dt < 4. && trk_qual > 0.065 && trk_pid > 0.49")
D0 = "trk_d0 > -75. && trk_d0 < 93."
D0_OPEN = "trk_d0 > -75."
COS_SPLIT, D0_SPLIT = 0.68, 100.

# name: (description, set, extra selection, categories, CRV veto from the tree)
TESTS = {
    "ref75":      ("set 75 (ConvAna)", "75", "1", "", False),
    "s89":        ("set 89 cuts", "81", f"{PRE} && {S89} && {D0}", "", True),
    "s89_d0open": ("set 89, d0 > -75 only", "81", f"{PRE} && {S89} && {D0_OPEN}", "", True),
    "cat_cos":    (f"set 89, cos > / <= {COS_SPLIT}", "81", f"{PRE} && {S89} && {D0}",
                   f"trk_cos > {COS_SPLIT}; trk_cos <= {COS_SPLIT}", True),
    "cat_d0":     (f"set 89 d0 open, d0 > / <= {D0_SPLIT:g}", "81", f"{PRE} && {S89} && {D0_OPEN}",
                   f"trk_d0 > {D0_SPLIT}; trk_d0 <= {D0_SPLIT}", True),
}
# The single-category test each category test is compared with
CATEGORY_BASE = {"cat_cos": "s89", "cat_d0": "s89_d0open"}


def split_name(kind, value):
    """Test name of a split scan point, e.g. cat_cos_0p72, cat_d0_125."""
    return f"cat_{kind}_{value:g}".replace(".", "p").replace("-", "m")


def add_split_tests(cos_splits, d0_splits):
    """Add the scan points: cat_cos at other cos(theta) splits (on s89), cat_d0 at other d0 splits (on s89_d0open).
    Returns {kind: [(split, test name), ...]} including the default split points."""
    scans = {"cos": [(COS_SPLIT, "cat_cos")], "d0": [(D0_SPLIT, "cat_d0")]}
    for kind, values, base, var in (("cos", cos_splits, "s89", "trk_cos"), ("d0", d0_splits, "s89_d0open", "trk_d0")):
        for v in values:
            if v == (COS_SPLIT if kind == "cos" else D0_SPLIT):
                continue
            name = split_name(kind, v)
            TESTS[name] = (f"{TESTS[base][0]}, {var} > / <= {v:g}", TESTS[base][1], TESTS[base][2],
                           f"{var} > {v}; {var} <= {v}", True)
            CATEGORY_BASE[name] = base
            scans[kind].append((v, name))
        scans[kind].sort()
    return scans

# One core per job: OpenBLAS (and ROOT's implicit MT) would otherwise use every core of the node
THREAD_ENV = {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}


def summary_path(work, test, target):
    return os.path.join(work, "test_hists", f"{test}_{target}", "summary.txt")


def read_summary(path):
    """(envelope limit, template Asimov limit) from a selection_test.py summary.txt, None if missing."""
    if not os.path.isfile(path):
        return None, None
    text = open(path).read()
    get = lambda label: (lambda m: float(m.group(1)) if m else None)(
        re.search(rf"{re.escape(label)}\s*:\s*([-\d.eE+]+)", text))
    return get("envelope fit to the MC-template Asimov"), get("histogram model on its own Asimov")


def run_test(work, test, target, force):
    _, sets, selection, categories, crv = TESTS[test]
    log = os.path.join(work, "logs", f"{test}_{target}.log")
    if not force and all(x is not None for x in read_summary(summary_path(work, test, target))):
        print(f"{test}_{target}: done earlier, reusing", flush=True)
        return 0
    cmd = [sys.executable, os.path.join(HERE, "selection_test.py"), selection, test, "-p", "mumem", "-s", sets,
           "--target", target, "--work-dir", work]
    if categories:
        cmd += ["--categories", categories]
    if crv:
        cmd.append("--crv-from-tree")
    print(f"{test}_{target}: running (log {log})", flush=True)
    with open(log, "w") as f:
        rc = subprocess.run(cmd, cwd=HERE, stdout=f, stderr=subprocess.STDOUT,
                            env={**os.environ, **THREAD_ENV}).returncode
    print(f"{test}_{target}: {'done' if rc == 0 else f'FAILED (rc {rc}), see {log}'}", flush=True)
    return rc


def pct(x, ref):
    return f"{100. * (x / ref - 1.):+.1f}" if x is not None and ref else "n/a"


def fmt(x):
    return f"{x:.3f}" if x is not None else "n/a"


def scan_tables(work, scans, tests, targets, results):
    """Limit vs split value for each target (template Asimov / envelope), with the single-category limit."""
    lines = []
    for kind, points in scans.items():
        points = [(v, t) for v, t in points if t in tests]
        if len(points) < 2:
            continue
        base = CATEGORY_BASE[points[0][1]]
        lines += ["", f"Split scan in {'cos(theta)' if kind == 'cos' else 'd0 (mm)'} on {base} "
                      f"(template Asimov / envelope limit; change vs {base} in %)"]
        head = f"{'split':>7s} " + " ".join(f"{g:>26s}" for g in targets)
        lines += [head, "-" * len(head)]
        rows = [(f"{base}", (lambda g: results.get((base, g), (None, None))))]
        rows += [(f"{v:g}", (lambda g, t=t: results.get((t, g), (None, None)))) for v, t in points]
        for label, get in rows:
            cells = []
            for g in targets:
                env, hists = get(g)
                b_env, b_hists = results.get((base, g), (None, None))
                cells.append(f"{fmt(hists)}/{fmt(env)} ({pct(hists, b_hists)})")
            lines.append(f"{label:>7s} " + " ".join(f"{c:>26s}" for c in cells))
        for g in targets:
            best = min(((results.get((t, g), (None, None))[1], v) for v, t in points
                        if results.get((t, g), (None, None))[1] is not None), default=None)
            if best:
                lines.append(f"  best {kind} split for {g}: {best[1]:g} (template Asimov limit {best[0]:.3f})")
    return lines


def tables(work, tests, targets, scans=None):
    results = {(t, g): read_summary(summary_path(work, t, g)) for t in tests for g in targets}
    rows = []
    for t in tests:
        for g in targets:
            env, hists = results[(t, g)]
            nom_env, nom_hists = results.get((t, "nominal"), (None, None))
            base = CATEGORY_BASE.get(t)
            b_env, b_hists = results.get((base, g), (None, None)) if base else (None, None)
            rows.append({"test": t, "target": g, "limit_hists": hists, "limit_env": env,
                         "vs_nominal_hists_pct": pct(hists, nom_hists), "vs_nominal_env_pct": pct(env, nom_env),
                         "vs_1cat_hists_pct": pct(hists, b_hists) if base else "",
                         "vs_1cat_env_pct": pct(env, b_env) if base else ""})
    with open(os.path.join(work, "results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    lines = ["Stopping-target sensitivity comparison (90% CL upper limit on r, reference 1e-15, Run 1A, unmixed)",
             "hists: template model on its own Asimov data; env: envelope model on that Asimov data",
             "vs nom: change relative to the nominal target with the same test; vs 1cat: category test relative "
             "to its single-category selection (cat_cos vs s89, cat_d0 vs s89_d0open), same target",
             ""]
    head = f"{'test':11s} {'target':15s} {'hists':>8s} {'env':>8s} {'vs nom h':>9s} {'vs nom e':>9s} " \
           f"{'vs 1cat h':>10s} {'vs 1cat e':>10s}"
    lines += [head, "-" * len(head)]
    for r in rows:
        lines.append(f"{r['test']:11s} {r['target']:15s} {fmt(r['limit_hists']):>8s} {fmt(r['limit_env']):>8s} "
                     f"{r['vs_nominal_hists_pct']:>9s} {r['vs_nominal_env_pct']:>9s} "
                     f"{r['vs_1cat_hists_pct']:>10s} {r['vs_1cat_env_pct']:>10s}")
        if r["target"] == targets[-1]:
            lines.append("")
    if scans:
        lines += scan_tables(work, scans, tests, targets, results) + [""]
    lines.append("Tests:")
    lines += [f"  {t:11s} {TESTS[t][0]}" + (f"  [categories: {TESTS[t][3]}]" if TESTS[t][3] else "") for t in tests]
    text = "\n".join(lines) + "\n"
    with open(os.path.join(work, "results.txt"), "w") as f:
        f.write(text)
    print(text)
    print(f"Tables written to {os.path.join(work, 'results.txt')} and results.csv")


def main():
    targets_all, _ = read_samples()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tests", default="", help="tests to run, quoted (default: all, with the scan points)")
    parser.add_argument("--targets", default=" ".join(targets_all), help="targets, quoted (default: all)")
    parser.add_argument("--work-dir", default=os.path.join(HERE, "alt_target_tests"),
                        help="output directory (default: %(default)s)")
    parser.add_argument("--jobs", type=int, default=5, help="parallel selection_test.py runs, one core each "
                                                            "(default: %(default)s, at most 20)")
    parser.add_argument("--cos-splits", default="", help="scan: further cos(theta) split points of cat_cos, quoted, "
                        "e.g. \"0.60 0.64 0.72 0.76\" (adds the tests cat_cos_<split>)")
    parser.add_argument("--d0-splits", default="", help="scan: further d0 split points (mm) of cat_d0, quoted, "
                        "e.g. \"50 75 125\" (adds the tests cat_d0_<split>)")
    parser.add_argument("--force", action="store_true", help="redo runs that already have a summary")
    parser.add_argument("--table-only", action="store_true", help="only remake the tables from existing runs")
    args = parser.parse_args()
    scans = add_split_tests([float(x) for x in args.cos_splits.split()], [float(x) for x in args.d0_splits.split()])
    tests, targets = (args.tests.split() or list(TESTS)), args.targets.split()
    for t in tests:
        if t not in TESTS:
            parser.error(f"unknown test {t} (available: {', '.join(TESTS)})")
    for g in targets:
        if g not in targets_all:
            parser.error(f"unknown target {g} (available: {', '.join(targets_all)})")
    if "nominal" not in targets:
        print("Warning: without the nominal target there is no reference for the changes")
    # a category test needs its single-category selection for the comparison
    tests += [CATEGORY_BASE[t] for t in tests if t in CATEGORY_BASE and CATEGORY_BASE[t] not in tests]
    tests = [t for t in TESTS if t in tests]
    work = os.path.abspath(args.work_dir)
    os.makedirs(os.path.join(work, "logs"), exist_ok=True)

    if not args.table_only:
        jobs = max(1, min(args.jobs, 20))
        runs = [(t, g) for g in targets for t in tests]
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            rcs = list(pool.map(lambda r: run_test(work, r[0], r[1], args.force), runs))
        failed = [f"{t}_{g}" for (t, g), rc in zip(runs, rcs) if rc != 0]
        if failed:
            print(f"Failed runs: {' '.join(failed)}")
    tables(work, tests, targets, scans)


if __name__ == "__main__":
    main()
