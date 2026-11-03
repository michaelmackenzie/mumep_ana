#! /usr/bin/env python3
"""Signal strength bias test: generate toys from one model and fit them with another.

Toys are generated from the generation card with Combine (GenerateOnly, --saveToys), fit with the
fit card (FitDiagnostics), and the fitted r and its pull are plotted with ../tools/plot_combine_fits.C.
With a single card, the toys are fit with the model they were generated from (a closure test).

build_model.py writes every model as workspaces on the channel observables (obs_<set>), so a model can
generate the toys another model of the same channels fits, e.g. the template model for an envelope:
  ./bias_test.py models/mumep_40_r0106_hists/combined.txt models/mumep_40_r0106_env/combined.txt
  ./bias_test.py models/mumem_75_r0106_hists/combined.txt -t 500 --expect-signal 10
"""
import argparse
import glob
import os
import shlex
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PLOT_MACRO = os.path.normpath(os.path.join(SCRIPT_DIR, '..', 'tools', 'plot_combine_fits.C'))

# Fit settings of the toy fits, as in legacy/combine_bias_test.sh: minimizer strategy 1 with a pre-fit
# scan, and the discrete profiling settings, which have no effect on a card without an envelope
BASE_FIT_ARGS = ['--cminDefaultMinimizerStrategy=1', '--cminApproxPreFitTolerance', '0.1', '--cminPreScan',
                 '--cminPreFit', '1',
                 '--cminRunAllDiscreteCombinations', '--X-rtd', 'REMOVE_CONSTANT_ZERO_POINT=1',
                 '--X-rtd', 'MINIMIZER_freezeDisassociatedParams',
                 '--X-rtd', 'MINIMIZER_multiMin_hideConstants']


def card_name(card):
    base = os.path.splitext(os.path.basename(card))[0]
    if base == 'combined':  # models/<name>/combined.txt, as written by build_model.py
        return os.path.basename(os.path.dirname(os.path.abspath(card)))
    for prefix in ['combine_total_', 'combine_']:
        if base.startswith(prefix):
            return base[len(prefix):]
    return base


def card_info(card):
    """Reference signal branching fraction (the unit of r) and the discrete parameters of a card"""
    unit, discrete = None, []
    with open(card) as f:
        for line in f:
            if line.startswith('# R_mue used for signal:'):
                unit = float(line.split(':')[1])
            tokens = line.split()
            if len(tokens) == 2 and tokens[1] == 'discrete':
                discrete.append(tokens[0])
    return unit, discrete


def run(command, log, cwd, dry_run):
    """Run a command list in cwd with its output logged; returns the exit status"""
    text = ' '.join(shlex.quote(c) for c in command)
    if dry_run:
        print(f'(cd {cwd} && {text})')
        return 0
    print(f'Running: {text}\n  log: {log}')
    with open(log, 'w') as f:
        f.write(text + '\n\n')
        f.flush()
        result = subprocess.run(command, cwd=cwd, stdout=f, stderr=subprocess.STDOUT, text=True)
    return result.returncode


def tail(log, n=30):
    with open(log) as f:
        return ''.join(f.readlines()[-n:])


parser = argparse.ArgumentParser(description='Generate toys from one card and fit them with another',
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
parser.add_argument('gen_card', help='Combine data card the toys are generated from (e.g. models/<name>/combined.txt)')
parser.add_argument('fit_card', nargs='?', default=None, help='Combine data card the toys are fit with (default: the generation card)')
parser.add_argument('-t', '--toys', type=int, default=1000, help='Number of toys (default: %(default)s)')
parser.add_argument('-s', '--seed', type=int, default=123456, help='Random seed for the toys (default: %(default)s)')
parser.add_argument('-T', '--tag', default=None,
                    help='Tag for the outputs (default: closure_<card name>, or bias_<gen name>_vs_<fit name>)')
parser.add_argument('--expect-signal', type=float, default=0.,
                    help='Signal strength r the toys are generated with, in units of the generation card (default: %(default)s)')
parser.add_argument('--r-range', type=float, default=800., help='Fit range of r, [-r, r] (default: %(default)s)')
parser.add_argument('--gen-params', default='', help='--setParameters for the generation, e.g. "mumep_40_env_cat=1"')
parser.add_argument('--fit-params', default='', help='--setParameters for the fit')
parser.add_argument('-G', '--gen-args', default='', help='Additional Combine arguments for the generation, e.g. "--toysFrequentist"')
parser.add_argument('-F', '--fit-args', default='', help='Additional Combine arguments for the fit')
parser.add_argument('--skip-bad-fits', type=int, default=0, choices=[0, 1, 2],
                    help='plot_combine_fits.C: 1 = drop fits with asymmetric errors, 2 = also drop failed fits (default: %(default)s)')
parser.add_argument('--err-mode', type=int, default=0, choices=[0, 1, 2],
                    help='plot_combine_fits.C pull error: 0 = asymmetric, 1 = rHiErr, 2 = symmetric rErr (default: %(default)s)')
parser.add_argument('-o', '--out-dir', default=None, help='Output directory (default: bias_tests/<tag>)')
parser.add_argument('--skip-generate', action='store_true', help='Reuse the toys of an earlier run in the output directory')
parser.add_argument('--plot-only', action='store_true', help='Only remake the figure from an earlier fit')
parser.add_argument('--dry-run', action='store_true', help='Print the commands without running them')
args = parser.parse_args()

if args.toys <= 0: parser.error('--toys must be positive')
if args.r_range <= 0.: parser.error('--r-range must be positive')

gen_card = os.path.abspath(args.gen_card)
fit_card = os.path.abspath(args.fit_card) if args.fit_card else gen_card
for c in [gen_card, fit_card]:
    if not os.path.isfile(c): parser.error(f'card {c} does not exist')

gen_name, fit_name = card_name(gen_card), card_name(fit_card)
tag = args.tag if args.tag else (f'closure_{fit_name}' if gen_card == fit_card else f'bias_{gen_name}_vs_{fit_name}')
out_dir = os.path.abspath(args.out_dir if args.out_dir else os.path.join('bias_tests', tag))

# r is in units of each card's signal reference branching fraction: express the generated signal in the
# units of the fit card, which is what the fitted r is compared with
gen_unit, gen_discrete = card_info(gen_card)
fit_unit, fit_discrete = card_info(fit_card)
r_true = args.expect_signal
if gen_unit and fit_unit and gen_unit != fit_unit:
    r_true = args.expect_signal*gen_unit/fit_unit
    print(f'WARNING: the cards use different signal reference branching fractions ({gen_unit:.3e} vs. {fit_unit:.3e}):'
          f' r = {args.expect_signal:g} generated is r = {r_true:g} in the fit card')
if gen_discrete and not any(gen_discrete[0] in p for p in args.gen_params.split(',')):
    print(f'NOTE: the generation card has a discrete envelope index ({", ".join(gen_discrete)}); the toys use its'
          ' default function unless --gen-params sets it')

toy_file = os.path.join(out_dir, f'higgsCombine.{tag}.GenerateOnly.mH120.{args.seed}.root')
fit_file = os.path.join(out_dir, f'fitDiagnostics.{tag}.root')
figure   = os.path.join(out_dir, f'{tag}.png')
log_dir  = os.path.join(out_dir, 'logs')

print(f'{args.toys} toys from {gen_card} (r = {args.expect_signal:g}), fit with {fit_card}'
      + (f' (envelope index {", ".join(fit_discrete)})' if fit_discrete else ''))
print(f'Outputs in {out_dir}')
if not args.dry_run:
    os.makedirs(log_dir, exist_ok=True)

# Single-threaded Combine jobs (OpenBLAS otherwise starts a thread per core)
os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')

if not args.plot_only:
    # Generate the toys
    if args.skip_generate:
        if not args.dry_run and not os.path.isfile(toy_file):
            exit(f'Error: toy file {toy_file} not found for --skip-generate')
    else:
        if not args.dry_run:
            for old in glob.glob(os.path.join(out_dir, f'higgsCombine.{tag}.GenerateOnly.*.root')):
                os.remove(old)
        command = ['combine', '-M', 'GenerateOnly', '-d', gen_card, '-t', str(args.toys), '--saveToys',
                   '--expectSignal', str(args.expect_signal), '-n', f'.{tag}', '-s', str(args.seed)]
        if args.gen_params: command += ['--setParameters', args.gen_params]
        command += shlex.split(args.gen_args)
        status = run(command, os.path.join(log_dir, 'generate.log'), out_dir, args.dry_run)
        if not args.dry_run and (status != 0 or not os.path.isfile(toy_file)):
            print(tail(os.path.join(log_dir, 'generate.log')))
            exit(f'Error: toy file {toy_file} was not created')

    # Fit them. Only the signal+background fit is used (tree_fit_sb), so the background-only fit is skipped.
    if not args.dry_run and os.path.isfile(fit_file):
        os.remove(fit_file)
    command = ['combine', '-M', 'FitDiagnostics', '-d', fit_card, f'--toysFile={toy_file}', '-t', str(args.toys),
               '-n', f'.{tag}', '-s', str(args.seed), '--rMin', str(-args.r_range), '--rMax', str(args.r_range),
               '--skipBOnlyFit'] + BASE_FIT_ARGS
    if args.fit_params: command += ['--setParameters', args.fit_params]
    command += shlex.split(args.fit_args)
    status = run(command, os.path.join(log_dir, 'fit.log'), out_dir, args.dry_run)
    if not args.dry_run and (status != 0 or not os.path.isfile(fit_file)):
        print(tail(os.path.join(log_dir, 'fit.log')))
        exit(f'Error: fit result file {fit_file} was not created')

# Plot the fitted r and its pull
if not args.dry_run and not os.path.isfile(fit_file):
    exit(f'Error: fit result file {fit_file} not found')
macro = f'{PLOT_MACRO}("{fit_file}", {r_true}, "{figure}", {args.skip_bad_fits}, {args.err_mode})'
status = run(['root.exe', '-q', '-b', '-l', macro], os.path.join(log_dir, 'plot.log'), out_dir, args.dry_run)
if not args.dry_run:
    print(tail(os.path.join(log_dir, 'plot.log'), 6))
    if status != 0 or not os.path.isfile(figure):
        exit(f'Error: the figure {figure} was not created')
    print(f'Figure: {figure}')
