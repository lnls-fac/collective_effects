"""."""

from .colleff import Ring as _Ring


def create_ring():
    """Create SIRIUS collective effects model ring.

    Returns:
        ring (pycolleff.colleff.Ring): main parameters for collective effects

    """
    ring = _Ring()
    ring.version = 'SI.v25.01-s05.02'
    ring.rf_freq = 499666600
    ring.mom_comp = 1.63e-4  # momentum compaction factor
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
    ring.gap_voltage = 1.75e6  # [V]
    return ring


def update_from_pymodels(ring, pyaccel_model=None):
    """Update pycolleff model with pymodels model equilibrium parameters.

    The gap voltage of the pycolleff model will be used in the process of
    getting the equilibrium parameters from the pyaccel model.

    The following attributes of the pycolleff model will be updated:
    damptx,
    dampty,
    dampte,
    mom_comp,
    bunlen,
    en_lost_rad,
    espread,
    sync_tune,
    tunex,   # Only the fractional part will be updated
    tuney    # Only the fractional part will be updated

    Args:
        ring (pycolleff.colleff.Ring): Collective effects model for sirius.
        pyaccel_model (pyaccel.accelerator.Accelerator, optional): Pyaccel
            model for Sirius. If None, the default model will be created
            internally. Defaults to None.

    """
    import pyaccel

    if pyaccel_model is None:
        from pymodels import si
        pyaccel_model = si.create_accelerator()

    idx = pyaccel.lattice.find_indices(
        pyaccel_model, 'frequency', 0, comparison=lambda x, y: x > y)[0]

    vgap0 = pyaccel_model[idx].voltage
    pyaccel_model[idx].voltage = ring.gap_voltage
    eqpar = pyaccel.optics.EqParamsFromBeamEnvelope(pyaccel_model)
    pyaccel_model[idx].voltage = vgap0

    ring.damptx = eqpar.tau1
    ring.dampty = eqpar.tau2
    ring.dampte = eqpar.tau3
    ring.mom_comp = -eqpar.etac
    ring.bunlen = eqpar.bunlen
    ring.en_lost_rad = eqpar.U0
    ring.espread = eqpar.espread0
    ring.sync_tune = eqpar.tune3
    ring.tunex = int(ring.tunex) + eqpar.tune1
    ring.tuney = int(ring.tuney) + eqpar.tune2
