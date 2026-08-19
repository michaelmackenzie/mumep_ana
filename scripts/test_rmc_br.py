# Evaluate the Phase-space model rates above different energy thresholds

# Integral of the normalized model from K_1 to K_2
def plestid_integral(K_1, K_2, KMax, knockout):
    if KMax <= 0.0:  return 0.0
    if knockout < 0: return 0.0
    K_1 = max(0.0, min(KMax, K_1))
    K_2 = max(0.0, min(KMax, K_2))
    
    if K_1 >= K_2: return 0.0
    power = 2.0 + 1.5 * knockout
    x_1 = K_1 / KMax
    x_2 = K_2 / KMax
    
    val_1 = (x_1 - 1.0) * pow(1.0 - x_1, power) * (power * x_1 + x_1 + 1.0)
    val_2 = (x_2 - 1.0) * pow(1.0 - x_2, power) * (power * x_2 + x_2 + 1.0)
    
    integral = val_2 - val_1
    return integral

# Begin main section
if __name__ == "__main__":
    kmax_0n = 101.8667 # MeV, no knockout endpoint on Al-27
    kmax_1n =  95.4489 # MeV, 1 neutron knockout endpoint on Al-27

    frac_0_0n  = plestid_integral( 0., kmax_0n, kmax_0n, 0)
    frac_0_1n  = plestid_integral( 0., kmax_1n, kmax_1n, 1)
    frac_57_0n = plestid_integral(57., kmax_0n, kmax_0n, 0)
    frac_57_1n = plestid_integral(57., kmax_1n, kmax_1n, 1)
    frac_80_0n = plestid_integral(80., kmax_0n, kmax_0n, 0)
    frac_80_1n = plestid_integral(80., kmax_1n, kmax_1n, 1)

    print(f'R(0 knockout | E >  0) / R(0 knockout) = {frac_0_0n:.5g}')
    print(f'R(1 knockout | E >  0) / R(1 knockout) = {frac_0_1n:.5g}')
    print(f'R(0 knockout | E > 57) / R(0 knockout) = {frac_57_0n:.5g}')
    print(f'R(1 knockout | E > 57) / R(1 knockout) = {frac_57_1n:.5g}')
    print(f'R(0 knockout | E > 80) / R(0 knockout) = {frac_80_0n:.5g}')
    print(f'R(1 knockout | E > 80) / R(1 knockout) = {frac_80_1n:.5g}')
