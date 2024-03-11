"""."""

import numpy as _np

from ..colleff import Ring as _Ring
from ..longitudinal_equilibrium import ImpedanceSource as _ImpSource


def create_ring():
    """Create HALF collective effects model ring.

    Based on data provided in [1].

    Returns:
        ring (pycolleff.colleff.Ring): main parameters for collective effects

    References:
    [1] 
    """
    ring = _Ring()
    ring.version = "HALF"
    # Circumference: 480 m
    ring.rf_freq = 499_654_096.67  # [Hz]
    ring.mom_comp = 8.1e-5  # momentum compaction factor
    ring.energy = 2.2e9  # energy [eV]
    ring.harm_num = 800  # harmonic Number
    ring.num_bun = 800  # number of bunches filled

    ring.total_current = 350e-3  # total current [A]
    ring.sync_tune = 0  # synchrotron tune
    ring.espread = 6.45e-4
    # Bunch duration: 6.76 ps
    ring.bunlen = 2.023e-3  # [m]
    ring.damptx = 0.0  # [s]
    ring.dampty = 0.0  # [s]
    ring.dampte = 22.7e-3  # [s]
    ring.en_lost_rad = 198.8e3  # [eV]
    ring.gap_voltage = 0.85e6  # [V]
    return ring


def create_harmonic_cavity():
    """Harmonic cavities parameters provided in [1]."""
    cav = _ImpSource()

    # High R/Q
    cav.Q = 5e5
    roverq = 90  # ohms
    cav.shunt_impedance = roverq * cav.Q  # ohms
    cav.harm_rf = 3
    freq_rf = 499_654_096.67
    cav.ang_freq_rf = 2 * _np.pi * freq_rf
    detune_f = 200e3
    cav.ang_freq = 2 * _np.pi * (cav.harm_rf * freq_rf + detune_f)
    return cav
