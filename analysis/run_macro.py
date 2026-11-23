#!/usr/bin/env python3
"""Run an analysis ROOT macro with the dataset list of combine/samples.py.

The macro is loaded into ROOT, samples.export_to_root sets its exposure (physics.C), histogram directories
(defaults.C), and dataset list (datasets.C), and the macro function is called. Use it in place of
root.exe -q -b '<macro>.C(<args>)':

    ./run_macro.py 'make_plots.C(true, {75})'
    ./run_macro.py --exposure run2 'make_plots.C(false, {40}, "", "run2")'
    ./run_macro.py --exposure mds3c 'make_plots.C(true, {75}, "mds3c")'
    ./run_macro.py --target st49hole1to5 'optimize_selection_fast.C(true, 81, false, "st49hole1to5", 1e3)'
    cd combine; ../run_macro.py 'make_selection_hists.C("40", "trk_qual > 0.9", "selection_hists/qual_0p9")'

The macros check that the exposure, 2BB setting, and target their tag/arguments name match the ones set here
(init_physics, check_dataset_exposure), so a mismatched --exposure, --2bb, or --target fails instead of
mis-normalizing. The macro path is relative to the working directory, as with root.exe. The exit status is 0 if the
macro returns 0.
"""

import argparse
import ctypes
import os
import re
import sys

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(HERE, "combine"))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("macro", help="macro call, '<file>.C(<args>)', or '<file>.C' to call it with no arguments")
    parser.add_argument("-e", "--exposure", default="run1a", help="samples.Exposure name (default: %(default)s)")
    parser.add_argument("--2bb", dest="two_bb", action="store_true",
                        default=os.environ.get("CONVANA_USE_2BB", "").lower() in ("1", "true", "yes"),
                        help="the 2BB exposure with the inputs re-weighted to 2BB (default: CONVANA_USE_2BB)")
    parser.add_argument("--target", default=None, help="stopping-target design of samples.TARGETS")
    parser.add_argument("--nominal-tracker", action="store_true",
                        help="use the nominal tracker samples instead of the reduced tracker ones")
    args = parser.parse_args()

    match = re.fullmatch(r"\s*([^()\s]+\.C)\s*(?:\((.*)\))?\s*", args.macro, re.S)
    if not match:
        parser.error(f"can't parse the macro call '{args.macro}': expected '<file>.C(<args>)'")
    path, call_args = match.group(1), match.group(2) or ""
    if not os.path.isfile(path):
        parser.error(f"macro {path} not found")

    import samples
    import ROOT
    ROOT.gROOT.SetBatch(True)
    exposure = samples.Exposure(args.exposure, two_bb=args.two_bb, target=args.target)

    if ROOT.gROOT.LoadMacro(path) != 0:
        sys.exit(f"run_macro.py: failed to load {path}")
    samples.export_to_root(exposure, reduced_tracker=not args.nominal_tracker)
    tracker = "nominal" if args.nominal_tracker else "reduced"
    print(f"run_macro.py: datasets from samples.py for {exposure} ({tracker} tracker)", flush=True)

    function = os.path.splitext(os.path.basename(path))[0]
    error = ctypes.c_int(0)
    status = ROOT.gROOT.ProcessLine(f"{function}({call_args});", error)
    sys.stdout.flush()
    if error.value != 0:
        sys.exit(f"run_macro.py: ROOT failed to run {function}({call_args}) (error {error.value})")
    if status != 0:
        print(f"run_macro.py: {function} returned {status}")
    sys.exit(0 if status == 0 else 1)


if __name__ == "__main__":
    main()
