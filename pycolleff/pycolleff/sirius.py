import numpy as _np
from . import colleff as colefcts
from . import impedances as imp
import mathphys as _mp

c = _mp.constants.light_speed

def create_ring(phase=2):
    ring = colefcts.Ring()
    ring.version   = 'SI.v20.01-s05.02'
    ring.circ      = 518.396
    ring.T0        = ring.circ/c
    ring.f0        = 1/ring.T0  # revolution frequency [Hz]
    ring.w0        = 2*_np.pi*ring.f0 # revolution angular frequency [rad/s]
    ring.mom_cmpct = 1.7e-4    # momentum compaction factor
    ring.E         = 3e9       # energy [eV]
    ring.nuy       = 14.165    # vertical tune
    ring.nux       = 49.110    # horizontal tune
    ring.chromx    = 2.5       # horizontal chromaticity
    ring.chromy    = 2.5       # vertical chromaticity
    ring.harm_num  = 864       # harmonic Number
    ring.nbun      = 864       # number of bunches filled
    # ring.budget    = create_budget(phase)

    I = _np.linspace(0,4,num=41)
    ring.cur_bun   = I*1e-3
    if phase == 0: #commissioning
        ring.version    += '.Commisioning'
        ring.nom_cur     = 0.100       # total current [A]
        ring.nus         = 0.0048      # synchrotron tune
        ring.espread     = lambda x:1e-2*(8.87e-2+ 1.58e-2*(x*1e3) - 5.48e-3*(x*1e3)**2 + 1.25e-3*(x*1e3)**3 - 1.14e-4*(x*1e3)**4)
        ring.bunlen      = lambda x:    (2.45e-3 + 5.80e-4*(x*1e3) - 1.51e-4*(x*1e3)**2 + 3.45e-5*(x*1e3)**3 - 3.15e-6*(x*1e3)**4)
        ring.emitx       = lambda x:1e-9*(2.44e-1+ 1.57e-1*(x*1e3) - 6.36e-2*(x*1e3)**2 + 1.60e-2*(x*1e3)**3 - 1.55e-3*(x*1e3)**4)
        ring.emity       = lambda x:1e-12*(2.15  + 1.87   *(x*1e3) - 8.49e-1*(x*1e3)**2 + 2.25e-1*(x*1e3)**3 - 2.25e-3*(x*1e3)**4)
        ring.damptx      = 17.1e-3
        ring.dampty      = 22.7e-3
        ring.dampte      = 13.6e-3
        ring.en_lost_rad = 475300 #eV
    elif phase == 1: #phase_1
        ring.version    += '.Phase1'
        ring.nom_cur     = 0.10         # total current [A]
        ring.nus         = 0.0048       # synchrotron tune
        ring.espread     = lambda x:1e-2*(8.87e-2+ 1.58e-2*(x*1e3) - 5.48e-3*(x*1e3)**2 + 1.25e-3*(x*1e3)**3 - 1.14e-4*(x*1e3)**4)
        ring.bunlen      = lambda x:    (2.45e-3 + 5.80e-4*(x*1e3) - 1.51e-4*(x*1e3)**2 + 3.45e-5*(x*1e3)**3 - 3.15e-6*(x*1e3)**4)
        ring.emitx       = lambda x:1e-9*(2.44e-1+ 1.57e-1*(x*1e3) - 6.36e-2*(x*1e3)**2 + 1.60e-2*(x*1e3)**3 - 1.55e-3*(x*1e3)**4)
        ring.emity       = lambda x:1e-12*(2.15  + 1.87   *(x*1e3) - 8.49e-1*(x*1e3)**2 + 2.25e-1*(x*1e3)**3 - 2.25e-3*(x*1e3)**4)
        ring.damptx      = 12.4e-3
        ring.dampty      = 15.1e-3
        ring.dampte      =  8.5e-3
        ring.en_lost_rad = 685374.1 #eV
    elif phase == 2: #phase_2
        ring.version    += '.Phase2'
        ring.nom_cur     = 0.35        # total current [A]
        ring.nus         = 0.0048      # synchrotron tune
        ring.espread     = lambda x:1e-2*(8.87e-2+ 1.58e-2*(x*1e3) - 5.48e-3*(x*1e3)**2 + 1.25e-3*(x*1e3)**3 - 1.14e-4*(x*1e3)**4)
        ring.bunlen      = lambda x:    (2.45e-3 + 5.80e-4*(x*1e3) - 1.51e-4*(x*1e3)**2 + 3.45e-5*(x*1e3)**3 - 3.15e-6*(x*1e3)**4)
        ring.emitx       = lambda x:1e-9*(2.44e-1+ 1.57e-1*(x*1e3) - 6.36e-2*(x*1e3)**2 + 1.60e-2*(x*1e3)**3 - 1.55e-3*(x*1e3)**4)
        ring.emity       = lambda x:1e-12*(2.15  + 1.87   *(x*1e3) - 8.49e-1*(x*1e3)**2 + 2.25e-1*(x*1e3)**3 - 2.25e-3*(x*1e3)**4)
        ring.damptx      = 10.6e-3
        ring.dampty      = 12.5e-3
        ring.dampte      =  6.9e-3
        ring.en_lost_rad = 829761.9 #eV
    elif phase == 3: #phase_2_HC
        ring.version    += '.Phase3'
        ring.nom_cur     = 0.5        # total current [A]
        ring.nus         = 0.00135    # synchrotron tune
        ring.espread     = lambda x:1e-2*(8.87e-2+1.58e-2*(x*1e3)-5.48e-3*(x*1e3)**2+1.25e-3*(x*1e3)**3-1.14e-4*(x*1e3)**4)
        ring.bunlen      = lambda x:12e-3    +0*(x*1e3)
        ring.emitx       = lambda x:1e-9*(1.89e-1+5.68e-2*(x*1e3)-1.59e-2*(x*1e3)**2+3.45e-3*(x*1e3)**3-3.10e-4*(x*1e3)**4)
        ring.emity       = lambda x:1e-12*(1.6497+1.04220*(x*1e3)-5.15e-1*(x*1e3)**2+1.45e-1*(x*1e3)**3-1.51e-2*(x*1e3)**4)
        ring.damptx      = 10.6e-3
        ring.dampty      = 12.5e-3
        ring.dampte      =  6.9e-3
        ring.en_lost_rad = 829761.9 #eV
    return ring

def create_budget(phase = 2):
    return imp.Budget()
