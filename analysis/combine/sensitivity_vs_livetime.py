import subprocess
import re
import matplotlib.pyplot as plt
import numpy as np

def parse_limit(text, level = r'50\.0', unit = 1.e-15):
    match = re.search(r"Expected "+level+r"%:\s*r\s*<\s*([\d\.-]+)", output_text)
    if match:
        extracted_value = unit*float(match.group(1))
        return extracted_value
    return -1.


# Input data card
signal   = 'mumep'
card_tag = 'evt_r0104'
card_set = '40'
card_base = f'{signal}_{card_set}_{card_tag}'
card = f'datacards/combine_total_{card_base}.txt'
r_range = 100. if 'mumem' in signal else 500.

# Base Combine command
base_command = f'combine -d {card} --rMin 0. --rMax {r_range} -n .{card_base} -t -1 --cl 0.9 --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rAbsAcc 0.001 --rRelAcc 0.001'

# 90% CL from SINDRUM II
published_limit = 7.e-13 if 'mumem' in signal else 1.7e-12

# Units for the input card
unit = 1.e-15

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
step = 0.1
index_nominal = -1
for index, scale in enumerate(np.arange(step, 2.+step, step)):
    if abs(scale - 1.) < step/2.: index_nominal = index
    command = base_command + f' --setParameters yieldScale={scale}'
    
    # Run the combine fit
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    output_text = result.stdout
    
    # Extract the median expected limit
    value_median = parse_limit(output_text, r'50\.0', unit)
    
    if value_median > 0.:
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
ax.set_xlabel("Relative running time")
ax.set_ylabel("Median expected upper limit")
ax.set_title("Expected 50.0% Upper Limits vs. Running Time")
ax.legend(loc="upper right")

# Save or show the plot
plt.savefig(f"sensitivity_vs_livetime_{card_base}.png", dpi=300)

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
