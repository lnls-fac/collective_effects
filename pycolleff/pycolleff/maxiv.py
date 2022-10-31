"""."""

from .colleff import Ring as _Ring


def create_ring(energy):
    """Create MAX-IV collective effects model ring.

    Args:
        energy (float): energy in [GeV]

    Raises:
        ValueError: if energy input is not 1.5 or 3 [GeV]

    Returns:
        ring (pycolleff.colleff.Ring): main parameters for collective effects

    """
    ring = _Ring()
    ring.energy = energy  # energy [eV]
    if energy == 3e9:
        ring.version = 'MAX-IV-3GeV'
        ring.rf_freq = 99.931e6
        ring.mom_comp = 3.06e-4  # momentum compaction factor

        ring.tuney = 42.20  # vertical tune
        ring.tunex = 16.28  # horizontal tune
        ring.chromx = 1.0  # horizontal chromaticity
        ring.chromy = 1.0  # vertical chromaticity
        ring.harm_num = 176  # harmonic Number
        ring.num_bun = 176  # number of bunches filled

        ring.total_current = 0.300  # total current [A]
        # ring.total_current = 0.250  # total current [A]
        ring.sync_tune = 0.001638  # synchrotron tune
        ring.espread = 7.69e-4
        ring.bunlen = 35.672e-12 * 299792458   # [m]
        ring.en_lost_rad = 363.8e3  # [eV]
        ring.gap_voltage = 1.397e6  # [V] for 300mA
        # ring.gap_voltage = 1.251e6  # [V] for 250mA

        alpha0 = ring.en_lost_rad/2/ring.rev_time/ring.energy
        jx = 1.847
        je = 3 - jx
        alphax = jx*alpha0
        alphae = je*alpha0
        ring.damptx = 1/alphax  # [s]
        ring.dampty = 1/alpha0  # [s]
        ring.dampte = 1/alphae  # [s]
    elif energy == 1.5e9:
        ring.version = 'MAX-IV-1.5GeV'
        ring.rf_freq = 99.931e6
        ring.mom_comp = 3.055e-3  # momentum compaction factor

        ring.tuney = 3.15  # vertical tune
        ring.tunex = 11.22  # horizontal tune
        ring.chromx = 1.0  # horizontal chromaticity
        ring.chromy = 1.0  # vertical chromaticity
        ring.harm_num = 32  # harmonic Number
        ring.num_bun = 32  # number of bunches filled

        ring.total_current = 0.500  # total current [A]
        ring.sync_tune = 0.002294   # synchrotron tune
        ring.espread = 7.41e-4
        ring.bunlen = 15.1e-3   # [m]
        ring.en_lost_rad = 114.4e3  # [eV]
        ring.gap_voltage = 520e3  # [V]

        alpha0 = ring.en_lost_rad/2/ring.rev_time/ring.energy
        jx = 1.46
        je = 3 - jx
        alphax = jx*alpha0
        alphae = je*alpha0
        ring.damptx = 1/alphax  # [s]
        ring.dampty = 1/alpha0  # [s]
        ring.dampte = 1/alphae  # [s]
    else:
        raise ValueError(
            'energy input for MAX-IV ring must be 1.5e9 or 3e9')
    return ring


def update_from_pymodels():
    """."""
    return None
