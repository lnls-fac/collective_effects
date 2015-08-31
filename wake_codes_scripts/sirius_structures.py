#!/usr/bin/env python3

import numpy as np

class GlobData:
    def __init__(self):
        self.ringpar = Ringpar()
        self.simpar  = Simpar()
        self.results = Results()

class Ringpar:
    def __init__(self):
        self.circumf = 518.396       # Ring circumference [m]
        self.omega0 = 3.63361516e6   # Revolution frequency [rad/s]
        self.sigma = 2.5e-3          # Nominal Bunch Length[m]
        self.sigmamax = 15e-3        # Maximum Bunch Length of interest[m]
        self.Iavg = 500e-3           # Beam current [A]
        self.h = 864                 # Harmonic number
        self.beta = 0.9999999855     # Beam speed (fraction of c)
        self.E0 = 3e9                # Beam Energy [eV]
        self.Ee = 0.510998910e6      # Electron Rest Energy [eV]
        self.alpha = 1.74e-4         # Momentum Compaction Factor
        self.nux = 48.0977           # Horizontal Tune
        self.omegas = 2*np.pi*299792458/518.396*0.00393 # Synchrotron Frequency [rad/s]

class Simpar:
    def __init__(self):
        self.wakepath  = ''    # Path where the wake file (from any of the softwares) is
        self.targetdir = ''    # Path where all results will be saved
        self.datasource= ''    # CST3PO or GdfidL
        self.m         = 0     # 0=longipole transuadrupole trans
        self.offset    = 0     # offset for transverse analysis
        self.sym       = ''    # mirror symmetry for transverse?
        self.whichaxis = ''    # y or x
        self.units     = 0     # # of components in the ring
        self.bunlen    = 0e-3  # Bunch Length Used in simulation[m]
        self.cutoff    = 2     # multiple of the bunch frequency to calculate impedance

class PeakInfo:
    def __init__(self):
        self.omegar   = np.array([],dtype=float)     # Modes Center Frequencies
        self.BW       = np.array([],dtype=float)     # Modes Bandwidths (half height)
        self.Rshunt   = np.array([],dtype=float)     # Modes Peak Values
        self.Q        = np.array([],dtype=float)     # Modes Quality Factors
        self.ReZmodel = np.array([],dtype=float)     # Impedance Spectrum Modeled from Peaks Info
        self.wmodel   = np.array([],dtype=float)     # Frequency axis referred to ReZmodel

class Results:
    def __init__(self):
        self.s            = np.array([],dtype=float) # axis: distance from following to drive bunch [m]
        self.W            = np.array([],dtype=float) # Longitudinal or Transverse Wakepotential [V/pC or V/pC/m]
        self.freq         = np.array([],dtype=float) # axis: frequency obtained from FFT [GHz]
        self.naxis        = np.array([],dtype=float) # axis: omega/omega0
        self.interfreq    = np.array([],dtype=float) # interpolated frequency
        self.ReZlong      = np.array([],dtype=float) # Real Part of Longitudinal Impedance [Ohm]
        self.interReZlong = np.array([],dtype=float) # interpolated impedance
        self.ImZlong      = np.array([],dtype=float) # Imaginary Part of Longitudinal Impedance [Ohm]
        self.ImZoN        = np.array([],dtype=float) # Imaginary Part of Longitudinal Impedance over n [Ohm]
        self.interReZt    = np.array([],dtype=float) # interpolated impedance
        self.ImZt         = np.array([],dtype=float) # Imaginary Part of Vertical Dipole Impedance [KOhm/m]
        self.ReZt         = np.array([],dtype=float) # Real Part of Horizontal Dipole Impedance [KOhm/m]
        self.interImZt    = np.array([],dtype=float) # interpolated impedance
        self.peakinfo     = np.array([],dtype=float) # Omegar [rad/s], Rshunt [Ohm] and Q from ReZ
        self.klossW       = np.array([],dtype=float) # Single-bunch Loss Factor Calculated from Wlong [mV/pC]
        self.klossZ       = np.array([],dtype=float) # Single-bunch Loss Factor Calculated from ReZlong [mV/pC]
        self.kickW        = np.array([],dtype=float) # Vertical Kick Factor Calculated from Wy [V/pC/m]
        self.kickZ        = np.array([],dtype=float) # Vertical Kick Factor Calculated from ImZy [V/pC/m]
        self.sigmak       = np.array([],dtype=float) # axis: bunch length for kloss|kick integration [mm]
        self.Ploss        = np.array([],dtype=float) # Power loss from single-bunch loss factor [W]
        self.Plossvec     = np.array([],dtype=float) # Power loss vector for different sigmas [W]
        self.klossWM      = np.array([],dtype=float) # Multi-bunch Loss Factor Calculated from Wlong [mV/pC]
        self.klossZM      = np.array([],dtype=float) # Multi-bunch Loss Factor Calculated from ReZlong [mV/pC]
        self.PlossM       = np.array([],dtype=float) # Power loss from multi-bunch loss factor [W]
        self.ifast        = np.array([],dtype=float) # # of fastest CBM
        self.GRs          = np.array([],dtype=float) # Growth Rate value for each CBM
        self.GR_HOM       = np.array([],dtype=float) # Growth Rate value for each CBM accumulated through each HOM
        self.ReZsampl     = np.array([],dtype=float) # Impedance Spectrum Sampled by fastest CBM
        self.fsampl       = np.array([],dtype=float) # Frequency axis for sampled impedance
        self.peakinfo     = PeakInfo()
