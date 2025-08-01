"""."""

import numpy as _np

from ..colleff import Ring as _Ring
from ..longitudinal_equilibrium import ImpedanceSource as _ImpSource


def create_ring(unstable_params=True):
    """Create ALS-U collective effects model ring.

    Based on data provided in [1-3].

    Returns:
        ring (pycolleff.colleff.Ring): main parameters for collective effects

    References:
    [1] Warnock, R. (2021). Equilibrium of an arbitrary bunch train with
    cavity resonators and short range wake: Enhanced iterative solution
    with Anderson acceleration. Physical Review Accelerators and Beams,
    24 (10), 104402. https://doi.org/10.1103/physrevaccelbeams.24.104402

    [2] Warnock, R. (2021). Equilibrium of an arbitrary bunch train in
    the presence of multiple resonator wakefields. Physical Review
    Accelerators and Beams, 24 (2), 1–27.
    https://doi.org/10.1103/PhysRevAccelBeams.24.024401

    [3] M. Venturini, "Passive higher-harmonic rf cavities with general
    settings and multibunch instabilities in electron storage rings"
    Phys. Rev. Accel. Beams 21, 114404 (2018)
    """
    ring = _Ring()
    ring.version = "ALS-U"
    if unstable_params:
        # From Ref. [3]
        ring.energy = 2e9  # energy [eV]
        ring.total_current = 500e-3  # total current [A]
        ring.mom_comp = 2.11e-4  # momentum compaction factor
        ring.espread = 0.943e-3
        ring.en_lost_rad = 0.217e6  # [eV]
        ring.harm_num = 328  # harmonic Number
        ring.num_bun = 328  # number of bunches filled
        ring.rf_freq = 500.417e6
        ring.gap_voltage = 0.6e6  # [V]
        ring.bunlen = 3.54e-3  # [m]
        ring.sync_tune = 1.75e-3  # synchrotron tune
        ring.damptx = 0.0  # [s]
        ring.dampty = 0.0  # [s]
        ring.dampte = 14e-3  # [s]

    else:
        ring.rf_freq = 500.390e6
        ring.mom_comp = 2.025e-4  # momentum compaction factor
        ring.energy = 2e9  # energy [eV]
        ring.harm_num = 328  # harmonic Number
        ring.num_bun = 328  # number of bunches filled

        ring.total_current = 500e-3  # total current [A]
        ring.sync_tune = 0  # synchrotron tune
        ring.espread = 1.02e-3
        ring.bunlen = 3.9e-3  # [m]
        ring.damptx = 0.0  # [s]
        ring.dampty = 0.0  # [s]
        ring.dampte = 0.0  # [s]
        ring.en_lost_rad = 315e3  # [eV]
        ring.gap_voltage = 600e3  # [V]
    return ring


def create_harmonic_cavity():
    """Harmonic cavities parameters provided in [1, 3]."""
    cav = _ImpSource()

    # Two ALS HHCs from Ref. [3]
    cav.shunt_impedance = 3.4e6
    cav.Q = 21_000
    cav.harm_rf = 3
    freq_rf = 500.417e6
    cav.ang_freq_rf = 2 * _np.pi * freq_rf
    detune_f = 584e3
    cav.ang_freq = 2 * _np.pi * (3 * freq_rf + detune_f)

    # # High R/Q
    # cav.shunt_impedance = 1.9e6
    # cav.Q = 2.4e4
    # cav.harm_rf = 3
    # freq_rf = 500.390e6
    # cav.ang_freq_rf = 2 * _np.pi * freq_rf
    # detune_f = 317.80e3
    # cav.ang_freq = 2 * _np.pi * (3 * freq_rf + detune_f)

    # # Low R/Q
    # cav.shunt_impedance = 1.4e6
    # cav.Q = 3.4e4
    # cav.harm_rf = 3
    # freq_rf = 500.390e6
    # cav.ang_freq_rf = 2*_np.pi*freq_rf
    # detune_f = 164.74e3
    # cav.ang_freq = 2*_np.pi*(3*freq_rf + detune_f)
    return cav


def create_main_cavity():
    """."""
    # Main cavities parameters provided in Ref. [2]

    # cav = _ImpSource()
    # cav.shunt_impedance = 0.8259e6
    # cav.Q = 3486
    # cav.harm_rf = 1
    # freq_rf = 500.390e6
    # cav.ang_freq_rf = 2*_np.pi*freq_rf
    # detune_f = -82.54e3
    # cav.ang_freq = 2*_np.pi*(freq_rf + detune_f)
    # cav.loop_ctrl_ang_freq = 2*_np.pi*freq_rf

    # Main cavities parameters provided in [1]
    cav = _ImpSource()
    beta = 9.983  # optimum
    # beta = 3.1  # ALS heritage
    Rs0 = 9.8e6
    cav.shunt_impedance = Rs0 / (1 + beta)
    Q0 = 3.6e4
    cav.Q = Q0 / (1 + beta)
    cav.harm_rf = 1
    freq_rf = 500.390e6
    cav.ang_freq_rf = 2 * _np.pi * freq_rf
    detune_f = -94.729e3
    cav.ang_freq = 2 * _np.pi * (freq_rf + detune_f)
    cav.loop_ctrl_ang_freq = 2 * _np.pi * freq_rf
    return cav
