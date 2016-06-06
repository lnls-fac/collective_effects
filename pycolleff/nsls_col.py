import numpy as _np
from . import collective_effects as colefcts
from . import impedances as imp
import mathphys as _mp

c = _mp.constants.light_speed

def create_ring():
    ring = colefcts.Ring()
    ring.version   = 'NSLS-II'
    ring.circ      = 791.958
    ring.T0        = ring.circ/c
    ring.f0        = 1/ring.T0  # revolution frequency [Hz]
    ring.w0        = 2*_np.pi*ring.f0 # revolution angular frequency [rad/s]
    ring.mom_cmpct = 3.6e-4    # momentum compaction factor
    ring.E         = 3e9       # energy [eV]
    ring.nuy       = 16.26     # vertical tune
    ring.nux       = 33.22     # horizontal tune
    ring.chromx    = 0.0       # horizontal chromaticity
    ring.chromy    = 0.0       # vertical chromaticity
    ring.harm_num  = 1320      # harmonic Number
    ring.nbun      = 1320      # number of bunches filled
    # ring.budget    = create_budget(phase)

    I = _np.linspace(0,4,num=40)
    ring.cur_bun     = I*1e-3
    ring.version    += 'Commisioning'
    ring.nom_cur     = 0.100       # total current [A]
    ring.nus         = 0.0067     # synchrotron tune
    ring.espread     = lambda x:7.64e-4 +0*x
    ring.sigma       = lambda x:3.5e-3  +0*x
    ring.emitx       = lambda x:271e-12 +0*x
    ring.emity       = lambda x:2.71e-12+0*x
    ring.damptx      = 55.3e-3
    ring.dampty      = 55.3e-3
    ring.dampte      = 27.7e-3
    ring.en_lost_rad = 286400.1 #eV

    return ring
