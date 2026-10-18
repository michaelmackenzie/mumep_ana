"""Expected sensitivity vs. a livetime scale factor (the yieldScale rateParam).

For each scale factor, a background-only Asimov data set is generated from a generation card
(default: the fit card) with yieldScale set to the scale factor, and the limit is evaluated on it
with the fit card. The generation card can be a different model of the same channels, e.g. the
histogram-template model generating the background that a discretely profiled envelope fits:
  python3 sensitivity_vs_livetime.py models/mumep_40_env/combined.txt \\
      --gen-card models/mumep_40_hists/combined.txt
build_model.py writes every model as workspaces on the channel observables (obs_<set>), so any of
its models can generate the data another one fits.
"""
import argparse
import glob
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

def run(command, log):
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
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
parser.add_argument('--log-dir', default='logs', help='Directory for the Combine logs (default: %(default)s)')
args = parser.parse_args()

step      = args.step
min_scale = args.min_scale if args.min_scale is not None else step
max_scale = args.max_scale
if step <= 0.: parser.error('--step must be positive')
if min_scale <= 0.: parser.error('--min-scale must be positive')
if max_scale < min_scale: parser.error('--max-scale must be >= --min-scale')

# Input data cards
card     = args.card
gen_card = args.gen_card if args.gen_card else card
for c in [card, gen_card]:
    if not os.path.isfile(c): parser.error(f'card {c} does not exist')
card_base = args.tag if args.tag else card_name(card)
signal = 'mumem' if 'mumem' in card_base else 'mumep'
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

print(f'Generating with {gen_card}, fitting with {card}' + (f' (envelope index {", ".join(discrete)})' if discrete else ''))
os.makedirs(args.log_dir, exist_ok=True)

# Lists to store the sensitivity vs. time
scales       = [] # Scale factor to livetime
values       = [] # Expected limit vs. livetime
ups_1        = [] # +1 sigma expectation
ups_2        = [] # +2 sigma expectation
downs_1      = [] # -1 sigma expectation
downs_2      = [] # -2 sigma expectation
improvements = [] # Improvement factor vs. livetime
bkg_only     = [] # Expected limit vs. background-only scale factor

# Loop over livetime scale factors
for scale in np.arange(min_scale, max_scale + step/2., step):
    name = f'livetime_{card_base}_{scale:.3f}'.replace('.', 'p')
    gen_params = ','.join(p for p in [f'yieldScale={scale}', args.gen_params] if p)
    fit_params = ','.join(p for p in [f'yieldScale={scale}', args.fit_params] if p)

    # Generate the Asimov data set
    for old in glob.glob(f'higgsCombine.{name}.GenerateOnly.mH120*.root'): os.remove(old)
    command = gen_command + f' -n .{name} --setParameters {gen_params}'
    status, output_text = run(command, os.path.join(args.log_dir, name + '_gen.log'))
    toy_files = glob.glob(f'higgsCombine.{name}.GenerateOnly.mH120*.root')
    if status != 0 or not toy_files:
        print(f"Scale {scale:.2f}: WARNING - Failed to generate the Asimov data set!")
        print(command)
        print(output_text)
        continue

    # Fit it, with the signal and the background scaled the same way in the fit model
    command = fit_command + f' -n .{name} --toysFile={toy_files[0]} --setParameters {fit_params}'
    status, output_text = run(command, os.path.join(args.log_dir, name + '_fit.log'))

    # Extract the median expected limit
    value_median = parse_observed(output_text, unit)

    if value_median > 0.:
        value_bands = parse_limit(output_text, r'50\.0', unit)
        if value_bands > 0. and abs(value_median/value_bands - 1.) > MAX_ASIMOV_MISMATCH:
            print(f"Scale {scale:.2f}: WARNING - Asimov limit {value_median:.3e} differs from the median of the bands"
                  f" {value_bands:.3e} by more than {100*MAX_ASIMOV_MISMATCH:.0f}%")
        scales.append(scale)
        values.append(value_median)
        improvements.append(published_limit/value_median)
        bkg_only.append(value_median*scale) # scale signal by 1/scale --> limit by scale
        ups_1.  append(parse_limit(output_text, r'84\.0', unit))
        ups_2.  append(parse_limit(output_text, r'97\.5', unit))
        downs_1.append(parse_limit(output_text, r'16\.0', unit))
        downs_2.append(parse_limit(output_text, r' 2\.5', unit))
        print(f"Scale {scale:.2f}: Limit = {value_median:.3e} Improvement = {published_limit/value_median:.1f}")
    else:
        print(f"Scale {scale:.2f}: WARNING - Failed to parse the median expected limit!")
        print(command)
        print(output_text)

if len(values) == 0:
    print("Failed to retrieve data!")
    exit()

# Plot limit vs. time
fig, ax = plt.subplots(figsize=(8, 5))

# Force all grid lines and ticks to render BELOW the plot elements (zorder < 1)
ax.set_axisbelow(True)
ax.grid(True, linestyle="--", alpha=0.6)

# 1. Plot the 2-sigma uncertainty band 
ax.fill_between(
    scales, downs_2, ups_2,
    color="#FFDF00", alpha=1.0, edgecolor="none", zorder=1,
    label=r"$\pm 2\sigma$ expected limit"
)

# 2. Plot the 1-sigma uncertainty band
ax.fill_between(
    scales, downs_1, ups_1,
    color="#00A86B", alpha=1.0, edgecolor="none", zorder=2,
    label=r"$\pm 1\sigma$ expected limit"
)

# 3. Plot the central median expected limit line on top
ax.plot(
    scales, values, 
    color="black", linewidth=2.0, zorder=3,
    label=r"CL$_{s}$ median expected limit"
)

# Apply clear labels and aesthetics
ax.set_xlabel("Relative running time", fontsize=14)
ax.set_ylabel("Median expected upper limit", fontsize=15)
ax.set_title("Expected Upper Limits vs. Running Time", fontsize=16)
ax.legend(loc="upper right", fontsize=14)

# Save or show the plot
plt.savefig(f"sensitivity_vs_livetime_{card_base}.png", dpi=300)
ax.set_yscale("log")
plt.savefig(f"sensitivity_vs_livetime_{card_base}_log.png", dpi=300)

# Plot improvement vs time
plt.figure(figsize=(8, 5))
plt.plot(scales, improvements, color="#1F77B4", linewidth=2.0, marker="o", markersize=5.0, label=r"Limit reduction factor")

# Apply clear labels and aesthetics
plt.xlabel("Relative running time")
plt.ylabel("Median expected upper limit improvement")
plt.title("Expected limit improvement vs. Running Time")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()

# Save or show the plot
plt.savefig(f"improvement_vs_livetime_{card_base}.png", dpi=300)

# Plot sensitivity vs background scaling
plt.figure(figsize=(8, 5))
plt.plot(scales, bkg_only, color="#1F77B4", linewidth=2.0, marker="o", markersize=5.0, label=r"CL$_{s}$ median expected limit")

# Apply clear labels and aesthetics
plt.xlabel("Relative background rate")
plt.ylabel("Median expected upper limit")
plt.title("Expected limit vs. Background rate")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()

# Save or show the plot
plt.savefig(f"sensitivty_vs_background_{card_base}.png", dpi=300)
