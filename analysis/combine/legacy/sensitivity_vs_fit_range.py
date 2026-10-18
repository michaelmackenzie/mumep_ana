"""Expected sensitivity vs. the momentum fit range of one set.

Starts from the shape card and workspace build_model writes for the full fit range, and makes a
model for each requested range with narrow_model.C instead of refitting the components: the PDFs
keep their full-range shapes, the rates are scaled to the fraction inside the range, and Combine
renormalizes the PDFs over the narrowed observable. Each range then gets an Asimov expected limit.

Example:
  python3 sensitivity_vs_fit_range.py datacards/combine_mumem_75_evt_r0105_hists.txt \\
      --low-edges "100 101 102 103" --high-edges "105 110"
  # remake the figures from the results table of an earlier scan
  python3 sensitivity_vs_fit_range.py datacards/combine_mumem_75_evt_r0105_hists.txt --plot-only
"""
import argparse
import csv
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

QUANTILES = {'2.5': r' 2\.5', '16.0': r'16\.0', '50.0': r'50\.0', '84.0': r'84\.0', '97.5': r'97\.5'}


def parse_observed(text):
    """Observed limit; with -t -1 this is the limit on the Asimov data set, used as the median"""
    match = re.search(r'Observed Limit:\s*r\s*<\s*([\d\.eE+-]+)', text)
    return float(match.group(1)) if match else -1.


def range_name(lo, hi):
    return f'range_{lo:.2f}_{hi:.2f}'.replace('.', 'p')


def parse_limit(text, level):
    match = re.search(r'Expected ' + level + r'%:\s*r\s*<\s*([\d\.eE+-]+)', text)
    return float(match.group(1)) if match else -1.


def parse_edges(text):
    return sorted({float(x) for x in text.replace(',', ' ').split()})


def card_info(card):
    """Reference signal branching fraction and the shape card's observable and workspace"""
    ref_br, shapes = None, None
    with open(card) as f:
        for line in f:
            if line.startswith('# R_mue used for signal:'):
                ref_br = float(line.split(':')[1])
            tokens = line.split()
            if len(tokens) >= 5 and tokens[0] == 'shapes' and tokens[1] == '*':
                shapes = tokens
    if shapes is None:
        raise SystemExit(f'{card} is not a single-channel shape card written by build_model')
    if ref_br is None:
        raise SystemExit(f'{card} has no "# R_mue used for signal" line')
    # relative to the working directory, the card, or the directory above the card, as in narrow_model.C
    ws_file = shapes[3]
    if not os.path.isabs(ws_file) and not os.path.exists(ws_file):
        card_dir = os.path.dirname(card)
        for alt in [os.path.join(card_dir, ws_file), os.path.join(os.path.dirname(card_dir), ws_file)]:
            if os.path.exists(alt):
                ws_file = alt
                break
    return ref_br, shapes[2], ws_file


def obs_range(ws_file, obs_name):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile.Open(ws_file)
    obs = f.Get('workspace').var(obs_name)
    res = (obs.getMin(), obs.getMax(), (obs.getMax() - obs.getMin()) / obs.getBins())
    f.Close()
    return res


def run(command, cwd, log):
    with open(log, 'w') as f:
        result = subprocess.run(command, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        f.write(result.stdout)
    return result.returncode, result.stdout


TABLE_NAME = 'sensitivity_vs_fit_range.csv'
INK, MUTED_INK = '#1F1F1F', '#5C5C5C'


def write_table(results, table):
    with open(table, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['low', 'high', 'status', 'median'] + list(QUANTILES))
        writer.writeheader()
        for res in results:
            writer.writerow({key: res.get(key, '') for key in writer.fieldnames})


def read_table(table):
    """Results of an earlier scan. Tables from before the Asimov median was recorded get it from
    the scan's Combine logs, converted with the same factor as the table's expected quantiles."""
    results = []
    recovered, missing = 0, 0
    with open(table, newline='') as f:
        for row in csv.DictReader(f):
            res = {'low': float(row['low']), 'high': float(row['high']), 'status': row['status']}
            for quantile in QUANTILES:
                res[quantile] = float(row[quantile]) if row.get(quantile) not in (None, '') else -1.
            res['median'] = float(row['median']) if row.get('median') not in (None, '') else -1.
            if 'median' not in row and res['50.0'] > 0.:
                log = os.path.join(os.path.dirname(table), range_name(res['low'], res['high']) + '_combine.log')
                text = open(log).read() if os.path.exists(log) else ''
                observed, expected = parse_observed(text), parse_limit(text, QUANTILES['50.0'])
                if observed > 0. and expected > 0.:
                    res['median'] = observed * res['50.0'] / expected
                    recovered += 1
                else:
                    missing += 1
            results.append(res)
    if recovered or missing:
        print(f'Recovered the Asimov median of {recovered} ranges from the Combine logs'
              + (f', {missing} without a usable log are left out' if missing else ''))
    return results


# Largest relative difference between the Asimov (observed) limit and the expected median for a range
# to count as a valid result; a larger one means a profiled fit missed the global minimum
MAX_ASIMOV_MISMATCH = 0.05


def asimov_mismatch(res):
    return res.get('median', -1.) > 0. and res.get('50.0', -1.) > 0. and \
        abs(res['median']/res['50.0'] - 1.) > MAX_ASIMOV_MISMATCH


def valid(results):
    return [r for r in results if r.get('median', -1.) > 0. and not asimov_mismatch(r)]


def plot_vs_low(results, out_dir, card_base, max_curves=6):
    """Median expected limit (the Asimov limit) vs. the lower range edge, one curve per upper edge,
    with the expected +-1 and 2 sigma bands when there is a single curve"""
    points_all = valid(results)
    highs = sorted({r['high'] for r in points_all})
    if not highs:
        return None
    if len(highs) > max_curves:
        print(f'{len(highs)} upper edges: skipping the curves vs. the lower edge, see the 2D figure')
        return None
    fig, ax = plt.subplots(figsize=(8, 6))
    for hi in highs:
        points = sorted((r for r in points_all if r['high'] == hi), key=lambda r: r['low'])
        x = [r['low'] for r in points]
        if len(highs) == 1:
            ax.fill_between(x, [r['2.5'] for r in points], [r['97.5'] for r in points], color='#FFDF7F',
                            label=r'$\pm 2\sigma$ expected limit')
            ax.fill_between(x, [r['16.0'] for r in points], [r['84.0'] for r in points], color='#85D1FB',
                            label=r'$\pm 1\sigma$ expected limit')
        ax.plot(x, [r['median'] for r in points], marker='o', markersize=5.0, linewidth=2.0,
                label=f'Median expected limit (Asimov), p < {hi:g} MeV/c')
    ax.set_xlabel('Lower edge of the momentum fit range (MeV/c)', fontsize=15)
    ax.set_ylabel('Expected 90% CL upper limit', fontsize=15)
    ax.set_title(f'Expected limit vs. fit range ({card_base})')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    figure = os.path.join(out_dir, 'sensitivity_vs_fit_range.png')
    fig.savefig(figure)
    plt.close(fig)
    return figure


def bin_edges(centers):
    """Cell edges around each scanned edge value, halfway to its neighbors"""
    c = np.asarray(centers, dtype=float)
    mid = 0.5*(c[1:] + c[:-1])
    return np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])


def plot_2d(results, out_dir, card_base, max_annotated=150):
    """Median expected limit (the Asimov limit) vs. (lower edge, upper edge), with the best range outlined"""
    points = valid(results)
    lows = sorted({r['low'] for r in points})
    highs = sorted({r['high'] for r in points})
    if len(lows) < 2 or len(highs) < 2:
        return None
    z = np.full((len(highs), len(lows)), np.nan)
    for r in points:
        z[highs.index(r['high']), lows.index(r['low'])] = r['median']
    zmin, zmax = np.nanmin(z), np.nanmax(z)
    norm = LogNorm(vmin=zmin, vmax=zmax) if zmax > zmin else LogNorm(vmin=0.9*zmin, vmax=1.1*zmax)
    cmap = plt.get_cmap('Blues').copy()  # sequential: one hue, light (small limit) to dark (large)
    cmap.set_bad('white')                # ranges not evaluated

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    x_edges, y_edges = bin_edges(lows), bin_edges(highs)
    mesh = ax.pcolormesh(x_edges, y_edges, np.ma.masked_invalid(z), cmap=cmap, norm=norm,
                         edgecolors='white', linewidth=2.0)
    colorbar = fig.colorbar(mesh, ax=ax)
    colorbar.set_label('Median expected 90% CL upper limit (Asimov)', fontsize=13, color=INK)

    # cell values in a common power of ten, with the ink picked for the cell's background
    exponent = int(np.floor(np.log10(zmin)))
    if z.size <= max_annotated:
        for iy, hi in enumerate(highs):
            for ix, lo in enumerate(lows):
                if np.isnan(z[iy, ix]):
                    continue
                dark_cell = norm(z[iy, ix]) > 0.55
                ax.text(lo, hi, f'{z[iy, ix] / 10**exponent:.2f}', ha='center', va='center', fontsize=8,
                        color='white' if dark_cell else INK)

    # outline the best range
    iy, ix = np.unravel_index(np.nanargmin(z), z.shape)
    ax.add_patch(Rectangle((x_edges[ix], y_edges[iy]), x_edges[ix+1] - x_edges[ix], y_edges[iy+1] - y_edges[iy],
                           fill=False, edgecolor=INK, linewidth=2.0))
    best = (f'best: {lows[ix]:g} < p < {highs[iy]:g} MeV/c, '
            f'limit {z[iy, ix] / 10**exponent:.2f}$\\times 10^{{{exponent}}}$')

    ax.set_xlabel('Lower edge of the momentum fit range (MeV/c)', fontsize=15)
    ax.set_ylabel('Upper edge of the momentum fit range (MeV/c)', fontsize=15)
    note = f'cell values $\\times 10^{{{exponent}}}$; ' if z.size <= max_annotated else ''
    ax.set_title(f'Expected limit vs. fit range ({card_base})\n{note}{best}', fontsize=11, color=INK)
    ax.tick_params(colors=MUTED_INK)
    fig.tight_layout()
    figure = os.path.join(out_dir, 'sensitivity_vs_fit_range_2d.png')
    fig.savefig(figure)
    plt.close(fig)
    return figure


def make_plots(results, out_dir, card_base):
    for figure in [plot_vs_low(results, out_dir, card_base), plot_2d(results, out_dir, card_base)]:
        if figure:
            print(f'Figure written to {figure}')


parser = argparse.ArgumentParser(description='Expected sensitivity vs. the momentum fit range of one set')
parser.add_argument('card', help='Shape card of one set, written by build_model for the full fit range '
                                 '(e.g. datacards/combine_mumem_75_evt_r0105_hists.txt)')
parser.add_argument('--low-edges', default=None,
                    help='Lower range edges to scan, quoted (default: the full range minimum to its midpoint, in --step steps)')
parser.add_argument('--high-edges', default=None,
                    help='Upper range edges to scan, quoted (default: the full range maximum); one curve per edge')
parser.add_argument('--step', type=float, default=0.5, help='Default lower-edge step in MeV/c (default: %(default)s)')
parser.add_argument('--out-dir', default=None, help='Output directory (default: range_scans/<card name>)')
parser.add_argument('--jobs', type=int, default=1, help='Ranges to evaluate in parallel (default: %(default)s)')
parser.add_argument('--no-combine', action='store_true',
                    help='Only build the narrowed cards and workspaces and print the Combine commands')
parser.add_argument('--fit-args', default='', help='Additional Combine arguments for every fit, quoted')
parser.add_argument('--plot-only', action='store_true',
                    help=f'Only remake the figures from the {TABLE_NAME} of an earlier scan in the output directory')
args = parser.parse_args()
if args.jobs < 1 or args.jobs > 20:
    parser.error('--jobs must be between 1 and 20')
if args.plot_only and args.no_combine:
    parser.error('--plot-only and --no-combine are exclusive')

macro_dir = os.path.dirname(os.path.abspath(__file__))
card = os.path.abspath(args.card)
card_base = os.path.splitext(os.path.basename(card))[0]
if card_base.startswith('combine_'):
    card_base = card_base[len('combine_'):]
out_dir = os.path.abspath(args.out_dir or os.path.join('range_scans', card_base))
table = os.path.join(out_dir, TABLE_NAME)

if args.plot_only:
    if not os.path.exists(table):
        parser.error(f'no results table {table} to plot; run the scan first')
    results = read_table(table)
    print(f'Read {len(results)} ranges ({len(valid(results))} with a limit) from {table}')
    make_plots(results, out_dir, card_base)
    raise SystemExit(0)

os.makedirs(out_dir, exist_ok=True)
ref_br, obs_name, ws_file = card_info(card)
full_min, full_max, bin_width = obs_range(ws_file, obs_name)
lows = parse_edges(args.low_edges) if args.low_edges else \
    list(np.round(np.arange(full_min, 0.5*(full_min + full_max) + 1.e-6, args.step), 6))
highs = parse_edges(args.high_edges) if args.high_edges else [full_max]
ranges = [(lo, hi) for hi in highs for lo in lows if hi - lo >= bin_width - 1.e-6]
if not ranges:
    parser.error('no range with a lower edge below its upper edge')
print(f'Card {card}: {obs_name} in [{full_min}, {full_max}] MeV/c, {bin_width:.3g} MeV/c bins, '
      f'{len(ranges)} ranges, output in {out_dir}')

r_max = 100. if 'mumem' in card_base else 500.
# Minimizer strategy 1: with strategy 0, the profiled fits of the Asimov limit search can get stuck
# in a worse local minimum of the background (envelope) shape parameters after the large-r points,
# so q_mu never returns to 0 at small r and the Asimov limit collapses towards 0.
# --noFitAsimov: the expected bands use the same pre-fit Asimov data set as the Asimov (observed)
# limit. Otherwise they use an Asimov data set regenerated from a background-only refit, and with a
# discretely profiled envelope the two can follow different functions (in narrow ranges a simpler
# function can describe the pre-fit Asimov almost as well and wins on the penalty), so the median
# and the bands describe different background hypotheses.
combine_options = (f'--rMin 0. --rMax {r_max} -t -1 --cl 0.9 --cminDefaultMinimizerStrategy=1'
                   ' --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rAbsAcc 0.01 --rRelAcc 0.01'
                   ' --cminRunAllDiscreteCombinations --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams'
                   ' --X-rtd MINIMIZER_multiMin_hideConstants --noFitAsimov'
                   + (f' {args.fit_args}' if args.fit_args else ''))


def evaluate(lo_hi):
    lo, hi = lo_hi
    name = range_name(lo, hi)
    narrow_card = os.path.join(out_dir, name + '.txt')
    status, _ = run(f'root.exe -q -b -l \'{macro_dir}/narrow_model.C("{card}", {lo}, {hi}, "{narrow_card}")\'',
                    os.getcwd(), os.path.join(out_dir, name + '_narrow.log'))
    if status != 0 or not os.path.exists(narrow_card):
        return {'low': lo, 'high': hi, 'status': 'narrow_model failed'}
    command = f'combine -d {narrow_card} -n .{name} {combine_options}'
    if args.no_combine:
        print(command)
        return {'low': lo, 'high': hi, 'status': 'not run'}
    status, text = run(command, out_dir, os.path.join(out_dir, name + '_combine.log'))
    res = {'low': lo, 'high': hi, 'status': 'ok' if status == 0 else 'combine failed'}
    for quantile, level in QUANTILES.items():
        value = parse_limit(text, level)
        res[quantile] = ref_br*value if value > 0. else -1.
    # with -t -1 the observed limit is the limit on the Asimov data set: the median expectation
    value = parse_observed(text)
    res['median'] = ref_br*value if value > 0. else -1.
    if res['median'] <= 0.:
        res['status'] = 'no limit found'
    elif asimov_mismatch(res):
        res['status'] = 'asimov mismatch'
        print(f"WARNING: [{lo:g}, {hi:g}] MeV/c: Asimov limit {res['median']:.3e} differs from the expected median "
              f"{res['50.0']:.3e} by more than {100*MAX_ASIMOV_MISMATCH:.0f}%, left out of the figures; see {name}_combine.log")
    return res


with ThreadPoolExecutor(max_workers=args.jobs) as pool:
    results = list(pool.map(evaluate, ranges))

write_table(results, table)
for res in results:
    median = res.get('median', -1.)
    print(f"[{res['low']:7.2f}, {res['high']:7.2f}] MeV/c: " +
          (f'median expected limit = {median:.3e}' if median > 0. else res['status'])
          + (f" ({res['status']})" if median > 0. and res['status'] != 'ok' else ''))
print(f'Results written to {table}')
if not args.no_combine:
    make_plots(results, out_dir, card_base)
