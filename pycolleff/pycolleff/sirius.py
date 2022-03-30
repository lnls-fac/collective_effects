"""."""

from . import colleff as colefcts
from . import impedances as imp


def create_ring(phase=2):
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

    if phase == 0:  # commissioning
        ring.version += '.Commissioning'
        ring.total_current = 0.030  # total current [A]
        ring.sync_tune = 0.00356  # synchrotron tune
        ring.espread = 8.87e-4
        ring.bunlen = 2.45e-3  # [m]
        ring.damptx = 16.9e-3  # [s]
        ring.dampty = 22.0e-3  # [s]
        ring.dampte = 12.9e-3  # [s]
        ring.en_lost_rad = 473e3  # [eV]
    elif phase == 1:  # phase_1
        ring.version += '.Phase1'
        ring.total_current = 0.10  # total current [A]
        ring.sync_tune = 0.0046  # synchrotron tune
        ring.espread = 8.87e-4
        ring.bunlen = 2.45e-3
        ring.damptx = 12.4e-3
        ring.dampty = 15.1e-3
        ring.dampte = 8.5e-3
        ring.en_lost_rad = 685374.1  # eV
    elif phase == 2:  # phase_2
        ring.version += '.Phase2'
        ring.total_current = 0.35  # total current [A]
        ring.sync_tune = 0.0046  # synchrotron tune
        ring.espread = 8.87e-4
        ring.bunlen = 2.45e-3
        ring.damptx = 10.6e-3
        ring.dampty = 12.5e-3
        ring.dampte = 6.9e-3
        ring.en_lost_rad = 829761.9  # [eV]
    elif phase == 3:  # phase_2_HC
        ring.version += '.Phase3'
        ring.total_current = 0.5  # total current [A]
        ring.sync_tune = 0.00135  # synchrotron tune
        ring.espread = 8.87e-4
        ring.bunlen = 12e-3
        ring.damptx = 10.6e-3
        ring.dampty = 12.5e-3
        ring.dampte = 6.9e-3
        ring.en_lost_rad = 829761.9  # [eV]
    return ring


def create_budget(phase=2):
    return imp.Budget()
