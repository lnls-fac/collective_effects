import numpy as _np
from . import colleff as colefcts
from . import impedances as imp
import mathphys as _mp

c = _mp.constants.light_speed

def create_ring(phase=2):
    ring = colefcts.Ring()
    ring.version   = 'SI.v12.c02'
    ring.circ      = 518.396
    ring.T0        = ring.circ/c
    ring.f0        = 1/ring.T0  # revolution frequency [Hz]
    ring.w0        = 2*_np.pi*ring.f0 # revolution angular frequency [rad/s]
    ring.mom_cmpct = 1.7e-4    # momentum compaction factor
    ring.E         = 3e9       # energy [eV]
    ring.nuy       = 13.116    # vertical tune
    ring.nux       = 48.131    # horizontal tune
    ring.chromx    = 0.0       # horizontal chromaticity
    ring.chromy    = 0.0       # vertical chromaticity
    ring.harm_num  = 864       # harmonic Number
    ring.nbun      = 864       # number of bunches filled
    # ring.budget    = create_budget(phase)

    I = _np.linspace(0,4,num=40)
    ring.cur_bun   = I*1e-3
    if phase == 0: #commissioning
        ring.version    += '.Commisioning'
        ring.nom_cur     = 0.100       # total current [A]
        ring.nus         = 0.00435    # synchrotron tune
        ring.espread     = lambda x:7.64e-4 +0*x
        ring.sigma       = lambda x:3e-3    +0*x
        ring.emitx       = lambda x:271e-12 +0*x
        ring.emity       = lambda x:2.71e-12+0*x
        ring.damptx      = 17.1e-3
        ring.dampty      = 22.7e-3
        ring.dampte      = 13.6e-3
        ring.en_lost_rad = 456740.6 #eV
    elif phase == 1: #phase_1
        ring.version    += '.Phase1'
        ring.nom_cur     = 0.10         # total current [A]
        ring.nus         = 0.00435        # synchrotron tune
        ring.espread     = lambda x:1e-2*(9.4e-2+3.80e-2*x-1.83e-2*x**2+4.78e-3*x**3-4.73e-4*x**4)
        ring.sigma       = lambda x:3e-3    +0*x
        ring.emitx       = lambda x:1e-9*(2.3e-1+1.57e-1*x-6.36e-2*x**2+1.60e-2*x**3-1.55e-3*x**4)
        ring.emity       = lambda x:1e-12*(2.15 +1.87   *x-8.49e-1*x**2+2.25e-1*x**3-2.25e-3*x**4)
        ring.damptx      = 12.4e-3
        ring.dampty      = 15.1e-3
        ring.dampte      =  8.5e-3
        ring.en_lost_rad = 685374.1 #eV
    elif phase == 2: #phase_2
        ring.version    += '.Phase2'
        ring.nom_cur     = 0.35        # total current [A]
        ring.nus         = 0.00435    # synchrotron tune
        ring.espread     = lambda x:1e-2*(8.87e-2+1.58e-2*x-5.48e-3*x**2+1.25e-3*x**3-1.14e-4*x**4)
        ring.sigma       = lambda x:3e-3    +0*x
        ring.emitx       = lambda x:1e-9*(1.89e-1+5.61e-2*x-1.59e-2*x**2+3.44e-3*x**3-3.10e-4*x**4)
        ring.emity       = lambda x:1e-12*(1.6497+1.04220*x-5.15e-1*x**2+1.45e-1*x**3-1.51e-2*x**4)
        ring.damptx      = 10.6e-3
        ring.dampty      = 12.5e-3
        ring.dampte      =  6.9e-3
        ring.en_lost_rad = 829761.9 #eV
    elif phase == 3: #phase_2_HC
        ring.version    += '.Phase3'
        ring.nom_cur     = 0.5        # total current [A]
        ring.nus         = 0.00135    # synchrotron tune
        ring.espread     = lambda x:1e-2*(8.87e-2+1.58e-2*x-5.48e-3*x**2+1.25e-3*x**3-1.14e-4*x**4)
        ring.sigma       = lambda x:12e-3    +0*x
        ring.emitx       = lambda x:1e-9*(1.89e-1+5.68e-2*x-1.59e-2*x**2+3.45e-3*x**3-3.10e-4*x**4)
        ring.emity       = lambda x:1e-12*(1.6497+1.04220*x-5.15e-1*x**2+1.45e-1*x**3-1.51e-2*x**4)
        ring.damptx      = 10.6e-3
        ring.dampty      = 12.5e-3
        ring.dampte      =  6.9e-3
        ring.en_lost_rad = 829761.9 #eV
    return ring

def create_budget(phase = 2):
    return imp.Budget()
