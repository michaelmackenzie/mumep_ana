"""Expected sensitivity vs. a livetime scale factor (the yieldScale rateParam).

For each scale factor, a background-only Asimov data set is generated from a generation card
(default: the fit card) with yieldScale set to the scale factor, and the limit is evaluated on it
with the fit card. The generation card can be a different model of the same channels, e.g. the
histogram-template model generating the background that a discretely profiled envelope fits:
  python3 sensitivity_vs_livetime.py models/mumep_40_env/combined.txt \\
      --gen-card models/mumep_40_hists/combined.txt
build_model.py writes every model as workspaces on the channel observables (obs_<set>), so any of
its models can generate the data another one fits.

Outputs: the Combine files, logs (logs/) and the results table (results.json) go in
livetime_scans/<tag>/ (-o), and the figures in figures/livetime/<tag>/ (--fig-dir). --plot-only remakes the
figures from the results table of an earlier scan.

The livetime plots show the running time in weeks (scale factor x --nominal-weeks, default the 28 days
of Run 1A) and mark the nominal running time. -j N runs N scale points in parallel.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import glob
import json
import os
import subprocess
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def parse_limit(text, level = r'50\.0', unit = 1.e-15):
    match = re.search(r"Expected "+level+r"%:\s*r\s*<\s*([\d\.eE+-]+)", text)
    if match:
        extracted_value = unit*float(match.group(1))
        return extracted_value
    return -1.

def parse_observed(text, unit = 1.e-15):
    """Observed limit; with the Asimov data set as the data, this is the median expected limit"""
    match = re.search(r"Observed Limit:\s*r\s*<\s*([\d\.eE+-]+)", text)
    return unit*float(match.group(1)) if match else -1.

def card_name(card):
    base = os.path.splitext(os.path.basename(card))[0]
    if base == 'combined':  # models/<name>/combined.txt, as written by build_model.py
        return os.path.basename(os.path.dirname(os.path.abspath(card)))
    for prefix in ['combine_total_', 'combine_']:
        if base.startswith(prefix):
            return base[len(prefix):]
    return base

def card_info(card):
    """Reference signal branching fraction (units of r) and the discrete parameters of a card"""
    unit, discrete = 1.e-15, []
    with open(card) as f:
        for line in f:
            if line.startswith('# R_mue used for signal:'):
                unit = float(line.split(':')[1])
            tokens = line.split()
            if len(tokens) == 2 and tokens[1] == 'discrete':
                discrete.append(tokens[0])
    return unit, discrete

# Combine jobs run one thread each, so -j counts cores
SINGLE_THREAD_ENV = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')

def run(command, log, cwd):
    result = subprocess.run(command, shell=True, cwd=cwd, env=SINGLE_THREAD_ENV, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open(log, 'w') as f:
        f.write(command + '\n\n' + result.stdout)
    return result.returncode, result.stdout


# Command line arguments
parser = argparse.ArgumentParser(description='Expected sensitivity vs. livetime scale factor')
parser.add_argument('card', help='Combine data card used in the fit (e.g. datacards/combine_total_mumep_40_evt_r0105_env.txt)')
parser.add_argument('--gen-card', default=None, help='Combine data card used to generate the Asimov data set (default: the fit card)')
parser.add_argument('--step', type=float, default=0.1, help='Livetime scale factor step size (default: %(default)s)')
parser.add_argument('--min-scale', type=float, default=None, help='Minimum livetime scale factor (default: step size)')
parser.add_argument('--max-scale', type=float, default=2.0, help='Maximum livetime scale factor, inclusive (default: %(default)s)')
parser.add_argument('--gen-params', default='', help='Extra --setParameters for the generation, e.g. "mumep_40_env_cat=1"')
parser.add_argument('--fit-params', default='', help='Extra --setParameters for the fit')
parser.add_argument('--gen-args', default='', help='Additional arguments for the generation')
parser.add_argument('--fit-args', default='', help='Additional arguments for the fit')
parser.add_argument('-T', '--tag', default=None, help='Tag for the outputs (default: the fit card name)')
parser.add_argument('--nominal-weeks', type=float, default=4.,
                    help='Running time in weeks of the nominal (scale factor 1) exposure (default: %(default)s, Run 1A)')
parser.add_argument('-j', '--threads', type=int, default=1, help='Number of scale points to run in parallel (default: %(default)s)')
parser.add_argument('-o', '--out-dir', default=None,
                    help='Directory for the Combine outputs, logs and results table (default: livetime_scans/<tag>)')
parser.add_argument('--fig-dir', default=None, help='Directory for the figures (default: figures/livetime/<tag>)')
parser.add_argument('--log-dir', default=None, help='Directory for the Combine logs (default: <out-dir>/logs)')
parser.add_argument('--plot-only', action='store_true', help='Only remake the figures from <out-dir>/results.json')
args = parser.parse_args()

step      = args.step
min_scale = args.min_scale if args.min_scale is not None else step
max_scale = args.max_scale
if step <= 0.: parser.error('--step must be positive')
if min_scale <= 0.: parser.error('--min-scale must be positive')
if max_scale < min_scale: parser.error('--max-scale must be >= --min-scale')
if args.nominal_weeks <= 0.: parser.error('--nominal-weeks must be positive')
if args.threads < 1: parser.error('--threads must be >= 1')
if args.threads > 20: print(f'WARNING: {args.threads} parallel Combine jobs, more than 20 cores on a shared node')

# Input data cards
card     = os.path.abspath(args.card)   # absolute: Combine runs in the output directory
gen_card = os.path.abspath(args.gen_card) if args.gen_card else card
for c in [card, gen_card]:
    if not os.path.isfile(c): parser.error(f'card {c} does not exist')
card_base = args.tag if args.tag else card_name(card)

# Output directories
out_dir      = os.path.abspath(args.out_dir if args.out_dir else os.path.join('livetime_scans', card_base))
fig_dir      = os.path.abspath(args.fig_dir if args.fig_dir else os.path.join('figures', 'livetime', card_base))
log_dir      = os.path.abspath(args.log_dir if args.log_dir else os.path.join(out_dir, 'logs'))
results_file = os.path.join(out_dir, 'results.json')
signal = 'mumem' if 'mumem' in card_name(card) else 'mumep'   # from the card, not the --tag
r_range = 100. if 'mumem' in signal else 500.

# Units for the input card: r is in units of the signal branching fraction used in the card
unit, discrete = card_info(card)
gen_unit, _    = card_info(gen_card)
if gen_unit != unit:
    print(f'WARNING: the cards use different reference branching fractions ({gen_unit:.3e} vs. {unit:.3e}),'
          ' the signal is not generated in this background-only Asimov data set')

# Generate a background-only Asimov data set from the pre-fit generation model
gen_command = f'combine -M GenerateOnly -d {gen_card} -t -1 --expectSignal 0 --saveToys -s 123456 {args.gen_args}'

# Fit the Asimov data set with the fit card. The limit on the Asimov data set (the observed limit)
# is the median expected limit. The bands use the fit model's Asimov data set from a
# background-only fit to the generated one (no --noFitAsimov), so they describe the generated
# background even when the fit model differs from the generation model.
fit_command = (f'combine -M AsymptoticLimits -d {card} --rMin 0. --rMax {r_range} -t -1 --cl 0.9 --cminApproxPreFitTolerance 0.1'
               ' --cminPreScan --cminPreFit 1 --rAbsAcc 0.01 --rRelAcc 0.01')
if discrete:
    # Discretely profiled envelope: minimizer strategy 0 can get stuck in a worse local minimum of
    # the envelope shape parameters, and the envelope needs every function tried in each fit
    fit_command += (' --cminDefaultMinimizerStrategy=1 --cminRunAllDiscreteCombinations'
                    ' --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_freezeDisassociatedParams'
                    ' --X-rtd MINIMIZER_multiMin_hideConstants')
else:
    fit_command += ' --cminDefaultMinimizerStrategy=0'
fit_command += f' {args.fit_args}'

# Allowed difference between the limit on the Asimov data set and the median of the expected bands
MAX_ASIMOV_MISMATCH = 0.05

# 90% CL from SINDRUM II
published_limit = 7.e-13 if 'mumem' in signal else 1.7e-12


def scan_point(scale):
    """Generate and fit the Asimov data set at one livetime scale factor: (limits, messages)"""
    name = f'livetime_{card_base}_{scale:.3f}'.replace('.', 'p')
    gen_params = ','.join(p for p in [f'yieldScale={scale}', args.gen_params] if p)
    fit_params = ','.join(p for p in [f'yieldScale={scale}', args.fit_params] if p)
    messages = []

    # Generate the Asimov data set
    toy_pattern = os.path.join(out_dir, f'higgsCombine.{name}.GenerateOnly.mH120*.root')
    for old in glob.glob(toy_pattern): os.remove(old)
    command = gen_command + f' -n .{name} --setParameters {gen_params}'
    status, output_text = run(command, os.path.join(log_dir, name + '_gen.log'), out_dir)
    toy_files = glob.glob(toy_pattern)
    if status != 0 or not toy_files:
        messages += [f"Scale {scale:.2f}: WARNING - Failed to generate the Asimov data set!", command, output_text]
        return None, messages

    # Fit it, with the signal and the background scaled the same way in the fit model
    command = fit_command + f' -n .{name} --toysFile={toy_files[0]} --setParameters {fit_params}'
    status, output_text = run(command, os.path.join(log_dir, name + '_fit.log'), out_dir)

    # Extract the median expected limit
    value_median = parse_observed(output_text, unit)
    if value_median <= 0.:
        messages += [f"Scale {scale:.2f}: WARNING - Failed to parse the median expected limit!", command, output_text]
        return None, messages

    value_bands = parse_limit(output_text, r'50\.0', unit)
    if value_bands > 0. and abs(value_median/value_bands - 1.) > MAX_ASIMOV_MISMATCH:
        messages.append(f"Scale {scale:.2f}: WARNING - Asimov limit {value_median:.3e} differs from the median of the bands"
                        f" {value_bands:.3e} by more than {100*MAX_ASIMOV_MISMATCH:.0f}%")
    messages.append(f"Scale {scale:.2f} ({scale*args.nominal_weeks:.2f} weeks): Limit = {value_median:.3e}"
                    f" Improvement = {published_limit/value_median:.1f}")
    limits = dict(median = value_median,
                  up_1   = parse_limit(output_text, r'84\.0', unit),
                  up_2   = parse_limit(output_text, r'97\.5', unit),
                  down_1 = parse_limit(output_text, r'16\.0', unit),
                  down_2 = parse_limit(output_text, r' 2\.5', unit))
    return limits, messages

def report(scale, result):
    limits, messages = result
    print('\n'.join(messages), flush=True)
    return scale, limits

if args.plot_only:
    # Remake the figures from the results table of an earlier scan
    if not os.path.isfile(results_file): parser.error(f'--plot-only: {results_file} does not exist')
    with open(results_file) as f:
        points = [(point['scale'], point['limits']) for point in json.load(f)['points']]
    print(f'Read {len(points)} scale points from {results_file}')
else:
    print(f'Generating with {gen_card}, fitting with {card}' + (f' (envelope index {", ".join(discrete)})' if discrete else ''))
    print(f'Combine outputs in {out_dir}')
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # Loop over livetime scale factors, in parallel with --threads; each point has its own output names
    scan = [float(scale) for scale in np.arange(min_scale, max_scale + step/2., step)]
    if args.threads > 1:
        print(f'Running {len(scan)} scale points with {args.threads} parallel jobs')
        with ThreadPoolExecutor(max_workers=args.threads) as pool:
            futures = [(scale, pool.submit(scan_point, scale)) for scale in scan]
            points = [report(scale, future.result()) for scale, future in futures]
    else:
        points = [report(scale, scan_point(scale)) for scale in scan]
    points = [(scale, limits) for scale, limits in points if limits is not None]

    # Results table, for --plot-only
    with open(results_file, 'w') as f:
        json.dump(dict(fit_card=card, gen_card=gen_card, unit=unit, nominal_weeks=args.nominal_weeks,
                       points=[dict(scale=scale, limits=limits) for scale, limits in points]), f, indent=1)
    print(f'Results in {results_file}')

# Sensitivity vs. time
scales       = [scale for scale, _ in points]                             # Scale factor to livetime
weeks        = [scale*args.nominal_weeks for scale in scales]             # Running time in weeks
values       = [limits['median'] for _, limits in points]                 # Expected limit vs. livetime
ups_1        = [limits['up_1']   for _, limits in points]                 # +1 sigma expectation
ups_2        = [limits['up_2']   for _, limits in points]                 # +2 sigma expectation
downs_1      = [limits['down_1'] for _, limits in points]                 # -1 sigma expectation
downs_2      = [limits['down_2'] for _, limits in points]                 # -2 sigma expectation
improvements = [published_limit/value for value in values]                # Improvement factor vs. livetime
bkg_only     = [value*scale for scale, value in zip(scales, values)]      # scale signal by 1/scale --> limit by scale

if len(values) == 0:
    print("Failed to retrieve data!")
    exit()

# Plot style: inward ticks on all sides, minor ticks, no title, Mu2e label above the frame
plt.rcParams.update({
    'font.size': 13, 'axes.labelsize': 15, 'xtick.labelsize': 13, 'ytick.labelsize': 13, 'legend.fontsize': 12,
    'axes.linewidth': 1.0, 'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'xtick.minor.visible': True, 'ytick.minor.visible': True, 'xtick.major.size': 6, 'ytick.major.size': 6,
    'xtick.minor.size': 3, 'ytick.minor.size': 3, 'legend.frameon': False, 'mathtext.default': 'regular',
})
process   = r'\mu^{-} \rightarrow e^{-}' if signal == 'mumem' else r'\mu^{-} \rightarrow e^{+}'
rate      = r'R_{\mu^{-}e^{-}}' if signal == 'mumem' else r'R_{\mu^{-}e^{+}}'
limit_lbl = rf'Expected 90% CL upper limit on ${rate}$'
weeks_lbl = 'Running time [weeks]'
nominal   = args.nominal_weeks
NOMINAL_COLOR = '#555555'

FIGSIZE = (8, 6)   # standard 4:3 canvas; fixed margins (no tight bbox) keep every output exactly this size

def new_axes():
    fig, ax = plt.subplots(figsize=FIGSIZE)
    fig.subplots_adjust(left=0.14, right=0.96, bottom=0.12, top=0.92)
    ax.text(0.0, 1.015, r'$\bf{Mu2e}$ $\it{Simulation}$', transform=ax.transAxes, ha='left', va='bottom', fontsize=15)
    ax.text(1.0, 1.015, f'${process}$', transform=ax.transAxes, ha='right', va='bottom', fontsize=13)
    return fig, ax

def mark_nominal(ax, x, label):
    """Vertical line at the nominal exposure, named in the legend"""
    ax.axvline(x, color=NOMINAL_COLOR, linestyle=':', linewidth=1.8, zorder=4, label=label)

def limit_scale(limits):
    """Power of ten for the linear limit axes, so the axis label carries it instead of an offset text"""
    return 10.**np.floor(np.log10(max(limits)))

def power_label(factor):
    return rf'[$10^{{{int(round(np.log10(factor)))}}}$]'

def save(fig, name):
    os.makedirs(fig_dir, exist_ok=True)
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(fig_dir, f'{name}.{ext}'), dpi=300)

def limit_at_nominal():
    """Median expected limit at the nominal running time, if it is a scan point"""
    for w, v in zip(weeks, values):
        if abs(w - nominal) < 1.e-6*nominal: return v
    return None

# Expected limit vs. running time, with the +-1 and +-2 sigma bands
lin = limit_scale(ups_2)
fig, ax = new_axes()
ax.fill_between(weeks, np.array(downs_2)/lin, np.array(ups_2)/lin, color='#FFDF00', edgecolor='none', zorder=1, label=r'Expected $\pm$ 2 s.d.')
ax.fill_between(weeks, np.array(downs_1)/lin, np.array(ups_1)/lin, color='#00A86B', edgecolor='none', zorder=2, label=r'Expected $\pm$ 1 s.d.')
ax.plot(weeks, np.array(values)/lin, color='black', linestyle='--', linewidth=2.0, zorder=3, label=r'Median expected (CL$_{s}$)')
nominal_value = limit_at_nominal()
if nominal_value is not None:
    ax.plot([nominal], [nominal_value/lin], marker='o', markersize=8, color='black', markeredgecolor='white',
            markeredgewidth=1.5, linestyle='none', zorder=5)
mark_nominal(ax, nominal, f'Nominal ({nominal:g} weeks)')
ax.set_xlim(0., max(weeks)*1.02)
ax.set_xlabel(weeks_lbl, loc='right')
ax.set_ylabel(f'{limit_lbl} {power_label(lin)}', loc='top')
ax.legend(loc='upper right', title='Background-only Asimov', title_fontsize=12)
ax.set_ylim(0., 1.3*max(ups_2)/lin)   # room for the legend
save(fig, 'sensitivity_vs_livetime')
plt.close(fig)

# Same, log scale (in absolute units, no rescaling needed)
fig, ax = new_axes()
ax.fill_between(weeks, downs_2, ups_2, color='#FFDF00', edgecolor='none', zorder=1, label=r'Expected $\pm$ 2 s.d.')
ax.fill_between(weeks, downs_1, ups_1, color='#00A86B', edgecolor='none', zorder=2, label=r'Expected $\pm$ 1 s.d.')
ax.plot(weeks, values, color='black', linestyle='--', linewidth=2.0, zorder=3, label=r'Median expected (CL$_{s}$)')
if nominal_value is not None:
    ax.plot([nominal], [nominal_value], marker='o', markersize=8, color='black', markeredgecolor='white',
            markeredgewidth=1.5, linestyle='none', zorder=5)
mark_nominal(ax, nominal, f'Nominal ({nominal:g} weeks)')
ax.set_xlim(0., max(weeks)*1.02)
ax.set_yscale('log')
ax.set_ylim(0.7*min(downs_2), 4.*max(ups_2))
ax.set_xlabel(weeks_lbl, loc='right')
ax.set_ylabel(limit_lbl, loc='top')
ax.legend(loc='upper right', title='Background-only Asimov', title_fontsize=12)
save(fig, 'sensitivity_vs_livetime_log')
plt.close(fig)

# Improvement over the SINDRUM II limit vs. running time
fig, ax = new_axes()
ax.plot(weeks, improvements, color='black', linewidth=2.0, marker='o', markersize=6, markeredgecolor='white',
        markeredgewidth=1.0, zorder=3, label='Median expected limit')
mark_nominal(ax, nominal, f'Nominal ({nominal:g} weeks)')
ax.set_xlim(0., max(weeks)*1.02)
ax.set_ylim(0., 1.15*max(improvements))
ax.set_xlabel(weeks_lbl, loc='right')
ax.set_ylabel('SINDRUM II limit / expected limit', loc='top')
ax.legend(loc='lower right')
save(fig, 'improvement_vs_livetime')
plt.close(fig)

# Expected limit vs. background rate, at fixed signal exposure (signal scaled by 1/scale)
fig, ax = new_axes()
lin_bkg = limit_scale(bkg_only)
ax.plot(scales, np.array(bkg_only)/lin_bkg, color='black', linewidth=2.0, marker='o', markersize=6, markeredgecolor='white',
        markeredgewidth=1.0, zorder=3, label='Median expected limit')
mark_nominal(ax, 1., 'Nominal background')
ax.set_xlim(0., max(scales)*1.02)
ax.set_ylim(0., 1.15*max(bkg_only)/lin_bkg)
ax.set_xlabel('Background rate / nominal', loc='right')
ax.set_ylabel(f'{limit_lbl} {power_label(lin_bkg)}', loc='top')
ax.legend(loc='lower right', title=f'Nominal signal exposure ({nominal:g} weeks)', title_fontsize=12)
save(fig, 'sensitivity_vs_background')
plt.close(fig)
print(f'Figures in {fig_dir}')
