import subprocess
import re
import matplotlib.pyplot as plt
import numpy as np

# Input data card
card_tag = 'evt_r0101_cc'
card_set = '75'
card_base = f'{card_set}_{card_tag}'
card = f'datacards/combine_total_mumem_{card_base}.txt'

# Base Combine command
base_command = f'combine -d {card} --rMin -100. --rMax 100. -n .{card_base} -t -1 --cl 0.9 --cminDefaultMinimizerStrategy=0 --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --rAbsAcc 0.0001 --rRelAcc 0.001'

# 90% CL from SINDRUM II
published_limit = 7.e-13

# Units for the input card
unit = 1.e-15

# Lists to store the sensitivity vs. time
scales       = [] # Scale factor to livetime
values       = [] # Expected limit vs. livetime
improvements = [] # Improvement factor vs. livetime
bkg_only     = [] # Expected limit vs. background-only scale factor

# Loop over livetime scale factors
step = 0.1
for scale in np.arange(step, 2.+step, step):
    command = base_command + f' --setParameters yieldScale={scale}'
    
    # Run the combine fit
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    output_text = result.stdout
    
    # Extract the median expected limit
    match = re.search(r"Expected 50\.0%:\s*r\s*<\s*([\d\.-]+)", output_text)
    
    if match:
        extracted_value = unit*float(match.group(1))
        scales.append(scale)
        values.append(extracted_value)
        improvements.append(published_limit/extracted_value)
        bkg_only.append(extracted_value*scale) # scale signal by 1/scale --> limit by scale
        print(f"Scale {scale:.2f}: Limit = {extracted_value:.3e} Improvement = {published_limit/extracted_value:.1f}")
    else:
        print(f"Scale {scale:.2f}: WARNING - Failed to parse the expected limit!")
        print(output_text)

# Plot limit vs. time
plt.figure(figsize=(8, 5))
plt.plot(scales, values, color="#1F77B4", linewidth=2.0, marker="o", markersize=5.0, label=r"CL$_{s}$ median expected limit")

# Apply clear labels and aesthetics
plt.xlabel("Relative running time")
plt.ylabel("Median expected upper limit")
plt.title("Expected 50.0% Upper Limits vs. Running Time")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()

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
