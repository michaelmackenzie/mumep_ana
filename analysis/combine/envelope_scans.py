"""Likelihood scans of the signal strength for a background envelope, and their overlay.

For a card whose background is a discretely profiled envelope (RooMultiPdf), scan r with each
envelope function held fixed, then with the function index profiled, and plot the curves together:
the profiled curve has to follow the lowest single-function curve at every point.

Example:
  python3 envelope_scans.py datacards/combine_mumep_40_evt_r0104_env.txt -r 200
  python3 envelope_scans.py datacards/combine_mumep_40_evt_r0104_env.txt --plot-only
"""
import argparse
import os
import re
import shlex
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Fixed categorical order for the envelope functions (validated palette), never cycled
COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
INK, MUTED_INK, GRID = '#1F1F1F', '#5C5C5C', '#D9D9D9'

# Discrete profiling settings. Not MINIMIZER_multiMin_maskConstraints: masking the constraints for
# the loop over the envelope functions leaves an NLL offset behind whenever the parameters change
# during the loop, which accumulates over the points of the profiled-index scan.
BASE_FIT_ARGS = ['--cminDefaultMinimizerStrategy=1', '--cminApproxPreFitTolerance', '0.1', '--cminPreScan',
                 '--cminPreFit', '1',
                 '--X-rtd', 'REMOVE_CONSTANT_ZERO_POINT=1',
                 '--X-rtd', 'MINIMIZER_freezeDisassociatedParams',
                 '--X-rtd', 'MINIMIZER_multiMin_hideConstants',
                 '--saveNLL']


def card_base_name(card):
    base = os.path.splitext(os.path.basename(card))[0]
    return base[len('combine_'):] if base.startswith('combine_') else base


def default_index(card):
    """Envelope index from the card name: build_model names it <process>_<selection>_env_cat"""
    match = re.search(r'(mumem|mumep)_(\d+)', os.path.basename(card))
    return f'{match.group(1)}_{match.group(2)}_env_cat' if match else None


def card_contents(card):
    """Discrete parameters, the workspace file, and the observable of the card"""
    discrete, ws_file, obs = [], None, None
    with open(card) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) == 2 and tokens[1] == 'discrete':
                discrete.append(tokens[0])
            if len(tokens) >= 4 and tokens[0] == 'shapes' and tokens[1] == '*' and ws_file is None:
                obs, ws_file = tokens[2], tokens[3]
    # the workspace path is relative to the working directory, the card, or the directory above it
    if ws_file and not os.path.isabs(ws_file) and not os.path.exists(ws_file):
        card_dir = os.path.dirname(os.path.abspath(card))
        for alt in [os.path.join(card_dir, ws_file), os.path.join(os.path.dirname(card_dir), ws_file)]:
            if os.path.exists(alt):
                ws_file = alt
                break
    return discrete, ws_file, obs


def envelope_functions(ws_file, index, obs_name):
    """Names of the envelope functions selected by the index, in index order, and the floating
    parameters of each"""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
    f = ROOT.TFile.Open(ws_file)
    ws = f.Get('workspace') if f else None
    cat = ws.cat(index) if ws else None
    if not cat:
        return None
    names, params = None, None
    for pdf in ws.allPdfs():
        if pdf.InheritsFrom('RooMultiPdf') and pdf.dependsOn(cat):
            names = [pdf.getPdf(i).GetName() for i in range(pdf.getNumPdfs())]
            observables = ROOT.RooArgSet(ws.var(obs_name)) if obs_name and ws.var(obs_name) else ROOT.RooArgSet()
            params = []
            for i in range(pdf.getNumPdfs()):
                pars = pdf.getPdf(i).getParameters(observables)
                params.append([p.GetName() for p in pars if p.InheritsFrom('RooRealVar') and not p.isConstant()])
            break
    f.Close()
    return names, params


def scan_file(tag, label):
    return f'higgsCombine.{tag}_{label}.MultiDimFit.mH120.root'


def read_scan(path, index=None):
    """Best fit and scan points of a MultiDimFit grid scan: r, absolute NLL, and the chosen index.
    Combine stores plain NLL values (not 2 NLL): nll0 before the initial fit, nll = NLL(best fit) - nll0,
    and deltaNLL = NLL(point) - NLL(best fit)."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile.Open(path)
    tree = f.Get('limit') if f else None
    if not tree:
        raise SystemExit(f'No scan tree in {path}')
    has_index = bool(index) and bool(tree.GetBranch(index))
    rows = []
    for entry in tree:
        # with --saveNLL the absolute NLL is nll0 + nll + deltaNLL, comparable between the scans
        rows.append((entry.r, entry.nll0 + entry.nll + entry.deltaNLL, getattr(entry, index) if has_index else -1))
    f.Close()
    best, points = rows[0], sorted(rows[1:])
    return {'best_r': best[0], 'best_nll': best[1],
            'r': np.array([p[0] for p in points]), 'nll': np.array([p[1] for p in points]),
            'index': np.array([p[2] for p in points], dtype=int)}


def crossings(r, y, level):
    """r values where the curve crosses the level, by linear interpolation"""
    res = []
    for i in range(1, len(r)):
        if (y[i-1] - level)*(y[i] - level) < 0.:
            res.append(r[i-1] + (level - y[i-1])*(r[i] - r[i-1])/(y[i] - y[i-1]))
    return res


def run_combine(command, log, dry_run):
    print(' '.join(shlex.quote(c) for c in command))
    if dry_run:
        return
    with open(log, 'w') as f:
        status = subprocess.run(command, stdout=f, stderr=subprocess.STDOUT).returncode
    if status != 0:
        with open(log) as f:
            print(''.join(f.readlines()[-20:]))
        raise SystemExit(f'Combine failed (status {status}), see {log}')


def make_plot(card_base, labels, scans, total, figure, ymax):
    # common reference: the lowest NLL of all scans
    nll_min = min([s['best_nll'] for s in scans + [total]] + [s['nll'].min() for s in scans + [total]])
    fig, ax = plt.subplots(figsize=(8.5, 6))
    for i, (label, scan) in enumerate(zip(labels, scans)):
        color = COLORS[i % len(COLORS)]
        ax.plot(scan['r'], 2.*(scan['nll'] - nll_min), color=color, linewidth=2.0, label=label, zorder=2)
        ax.plot([scan['best_r']], [2.*(scan['best_nll'] - nll_min)], marker='*', markersize=12, color=color,
                markeredgecolor='white', markeredgewidth=1.0, linestyle='none', zorder=4)
    y_total = 2.*(total['nll'] - nll_min)
    ax.plot(total['r'], y_total, color=INK, linewidth=2.0, linestyle='--', label='Envelope (index profiled)', zorder=3)
    # the function the envelope chose at each point
    if (total['index'] >= 0).any():
        for i in sorted(set(total['index'])):
            sel = total['index'] == i
            ax.plot(total['r'][sel], y_total[sel], marker='o', markersize=5, linestyle='none', zorder=5,
                    color=COLORS[i % len(COLORS)], markeredgecolor='white', markeredgewidth=0.8)
    ax.plot([total['best_r']], [2.*(total['best_nll'] - nll_min)], marker='*', markersize=14, color=INK,
            markeredgecolor='white', markeredgewidth=1.0, linestyle='none', label='Best fits', zorder=6)
    for level, text in [(1., '68% CL'), (4., '95% CL')]:
        ax.axhline(level, color=MUTED_INK, linewidth=1.0, linestyle=':', zorder=1)
        ax.text(total['r'].min(), level, text, va='bottom', ha='left', fontsize=9, color=MUTED_INK)
    if ymax is None:
        # show the 95% CL band and every function's own minimum
        ymax = max(10., max(2.*(s['best_nll'] - nll_min) for s in scans) + 5.)
    ax.set_ylim(-0.05*ymax, ymax)
    ax.set_xlabel('Signal strength r', fontsize=14, color=INK)
    ax.set_ylabel(r'$-2\Delta\ln L$', fontsize=14, color=INK)
    ax.set_title(f'Envelope likelihood scans ({card_base})', fontsize=12, color=INK)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.tick_params(colors=MUTED_INK)
    ax.legend(fontsize=9, frameon=False, loc='upper center', ncol=2)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(figure)), exist_ok=True)
    fig.savefig(figure)
    plt.close(fig)
    print(f'Figure written to {figure}')
    return nll_min


def summarize(labels, scans, total, nll_min):
    print(f"Envelope best fit: r = {total['best_r']:.3f}")
    y_total = 2.*(total['nll'] - nll_min)
    for level, text in [(1., '68%'), (4., '95%')]:
        edges = ', '.join(f'{x:.2f}' for x in crossings(total['r'], y_total, level))
        print(f'  {text} CL crossings of the envelope (-2 Delta ln L = {level:g}): r = {edges or "none in the scan range"}')
    for label, scan in zip(labels, scans):
        print(f"  {label}: best fit r = {scan['best_r']:.3f}, -2 Delta ln L at its minimum = {2.*(scan['best_nll'] - nll_min):.3f}")
    # the envelope has to follow the lowest function at every point; compared at the r points all
    # scans have (a fit that fails at some point drops it from that scan)
    common = [r for r in total['r'] if all(np.isclose(s['r'], r, atol=1.e-3).any() for s in scans)]
    if common:
        def at(scan, r):
            return scan['nll'][np.argmin(np.abs(scan['r'] - r))]
        diff = np.array([at(total, r) - min(at(s, r) for s in scans) for r in common])
        worst = np.argmax(np.abs(diff))
        status = 'OK' if abs(diff[worst]) < 0.01 else 'MISMATCH'
        print(f"  Envelope vs. lowest function at {len(common)} / {len(total['r'])} points: largest NLL difference "
              f"{diff[worst]:+.4f} at r = {common[worst]:.2f} [{status}]")
    else:
        print('  The scans share no r points, envelope vs. lowest function not compared')
    missing = [f'{label} ({len(total["r"]) - len(s["r"])})' for label, s in zip(labels, scans) if len(s['r']) < len(total['r'])]
    if missing:
        print(f"  Scan points lost to failed fits: {', '.join(missing)}")


parser = argparse.ArgumentParser(description='Likelihood scans of r for each envelope function and with the index profiled')
parser.add_argument('card', help='Combine card with a discretely profiled envelope (e.g. datacards/combine_mumep_40_evt_r0104_env.txt)')
parser.add_argument('-i', '--index', default=None,
                    help='Envelope index parameter (default: <process>_<selection>_env_cat from the card name)')
parser.add_argument('-n', '--npdfs', type=int, default=None,
                    help='Number of envelope functions (default: read from the workspace)')
parser.add_argument('-r', '--range', type=float, default=100., help='Scan r in [-range, range] (default: %(default)s)')
parser.add_argument('--points', type=int, default=40, help='Scan points (default: %(default)s)')
parser.add_argument('-T', '--tag', default=None, help='Tag of the Combine outputs (default: envelope_<card name>)')
parser.add_argument('-o', '--figure', default=None, help='Output figure (default: envelope_<card name>.png)')
parser.add_argument('--fit-args', default='', help='Additional Combine arguments for every fit, quoted')
parser.add_argument('--ymax', type=float, default=None,
                    help="Upper edge of the -2 Delta ln L axis (default: at least 10, and above every function's minimum)")
parser.add_argument('--plot-only', action='store_true', help='Only remake the figure from the existing scans')
parser.add_argument('--dry-run', action='store_true', help='Only print the Combine commands')
args = parser.parse_args()

if not os.path.isfile(args.card):
    parser.error(f'card {args.card} does not exist')
card_base = card_base_name(args.card)
tag = args.tag or f'envelope_{card_base}'
figure = args.figure or f'envelope_{card_base}.png'

discrete, ws_file, obs_name = card_contents(args.card)
index = args.index or default_index(args.card)
if not index:
    parser.error('no envelope index in the card name, pass it with --index')
if discrete and index not in discrete:
    parser.error(f'index {index} is not a discrete parameter of the card ({", ".join(discrete)}), pass it with --index')

labels, function_params = None, None
if ws_file and os.path.exists(ws_file):
    labels, function_params = envelope_functions(ws_file, index, obs_name)
npdfs = args.npdfs or (len(labels) if labels else None)
if not npdfs:
    parser.error('could not read the envelope functions from the workspace, pass their number with --npdfs')
if not labels or len(labels) != npdfs:
    labels = [f'PDF {i}' for i in range(npdfs)]
print(f'Card {args.card}: envelope index {index} with {npdfs} functions: {", ".join(labels)}')

if not args.plot_only:
    grid = ['combine', '-M', 'MultiDimFit', args.card, '--algo', 'grid', '--points', str(args.points),
            '--rMin', str(-args.range), '--rMax', str(args.range)]
    extra = shlex.split(args.fit_args)
    for i in range(npdfs):
        print(f'Scan with envelope function {i} ({labels[i]})')
        # With the index frozen, MINIMIZER_freezeDisassociatedParams leaves the parameters of the
        # other functions floating; the likelihood does not depend on them, so their zero curvature
        # makes the strategy 1 fits fail. Freeze them explicitly.
        freeze = [index]
        if function_params and len(function_params) == npdfs:
            freeze += [p for j, pars in enumerate(function_params) if j != i for p in pars if p not in function_params[i]]
        run_combine(grid + ['-n', f'.{tag}_pdf{i}', '--setParameters', f'{index}={i}', '--freezeParameters', ','.join(freeze)]
                    + BASE_FIT_ARGS + extra, f'{tag}_pdf{i}.log', args.dry_run)
    # the index floats, so the envelope picks the best function at each point; the choice is saved
    print('Scan with the envelope index profiled')
    run_combine(grid + ['-n', f'.{tag}_total', '--cminRunAllDiscreteCombinations', '--saveSpecifiedIndex', index]
                + BASE_FIT_ARGS + extra, f'{tag}_total.log', args.dry_run)
    if args.dry_run:
        raise SystemExit(0)

missing = [scan_file(tag, l) for l in [f'pdf{i}' for i in range(npdfs)] + ['total'] if not os.path.exists(scan_file(tag, l))]
if missing:
    raise SystemExit(f'Missing scan outputs: {", ".join(missing)}')
scans = [read_scan(scan_file(tag, f'pdf{i}')) for i in range(npdfs)]
total = read_scan(scan_file(tag, 'total'), index)
nll_min = make_plot(card_base, labels, scans, total, figure, args.ymax)
summarize(labels, scans, total, nll_min)
