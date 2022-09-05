"""."""

from .colleff import Ring as _Ring


def create_ring():
    """."""
    ring = _Ring()
    ring.version = 'MAX-IV'
    ring.rf_freq = 99.931e6
    ring.mom_comp = 3.06e-4  # momentum compaction factor
    ring.energy = 3e9  # energy [eV]
    ring.tuney = 0  # vertical tune
    ring.tunex = 0  # horizontal tune
    ring.chromx = 0  # horizontal chromaticity
    ring.chromy = 0  # vertical chromaticity
    ring.harm_num = 176  # harmonic Number
    ring.num_bun = 176  # number of bunches filled

    ring.total_current = 0.300  # total current [A]
    # ring.total_current = 0.250  # total current [A]
    ring.sync_tune = 0.001638  # synchrotron tune
    ring.espread = 7.69e-4
    ring.bunlen = 35.672e-12 * 299792458   # [m]
    ring.damptx = 0  # [s]
    ring.dampty = 0  # [s]
    ring.dampte = 25.194e-3  # [s]
    ring.en_lost_rad = 363.8e3  # [eV]
    ring.gap_voltage = 1.397e6  # [V]
    # ring.gap_voltage = 1.251e6  # [V]
    return ring


def update_from_pymodels():
    """."""
    return None
