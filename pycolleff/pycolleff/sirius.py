"""."""

from . import colleff as colefcts
from . import impedances as imp


def create_ring():
    ring = colefcts.Ring()
    ring.version = 'SI.v20.01-s05.02'
    ring.rf_freq = 499666600
    ring.mom_comp = 1.7e-4  # momentum compaction factor
    ring.energy = 3e9  # energy [eV]
    ring.tuney = 14.137  # vertical tune
    ring.tunex = 49.078  # horizontal tune
    ring.chromx = 2.5  # horizontal chromaticity
    ring.chromy = 2.5  # vertical chromaticity
    ring.harm_num = 864  # harmonic Number
    ring.num_bun = 864  # number of bunches filled

    ring.total_current = 0.10  # total current [A]
    ring.sync_tune = 0.00356  # synchrotron tune
    ring.espread = 8.87e-4
    ring.bunlen = 3.25e-3  # [m]
    ring.damptx = 16.9e-3  # [s]
    ring.dampty = 22.0e-3  # [s]
    ring.dampte = 12.9e-3  # [s]
    ring.en_lost_rad = 473e3  # [eV]
    return ring
