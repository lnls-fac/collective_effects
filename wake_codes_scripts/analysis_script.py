#!/usr/bin/env python3
import sys
sys.path.append('/home/fac_files/code/collective_effects/wake_codes_scripts/')
import functions_for_analysis as funcs

newdir = ''
m = 1
bunlen = 0.5e-3
globdata = funcs.prepare_struct_for_load(newdir, m, bunlen)

# Load wakepotential result from referred software, rescale and save
#  txt-file on default format
globdata = funcs.load_wake(globdata)

# Calculates Impedance Spectrum from Wakepotential Results
globdata = funcs.calc_impedance(globdata)

# Calculates Loss Factor
if m == 0:
    globdata = funcs.calc_loss_factor(globdata)
elif m > 0:
    globdata = funcs.calc_kick_factor(globdata)

# Plot Results
funcs.plot_results(globdata,mostra=True)

# Export Results
funcs.save_results(globdata)
