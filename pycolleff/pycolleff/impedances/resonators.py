"""."""

import numpy as _np
import mathphys as _mp

_LSPEED = _mp.constants.light_speed


def longitudinal_resonator(Rs, Q, wr, w):
    """Return the longitudinal resonator impedance for w.

    Inputs:
    Rs    = Shunt impedance [Ohm]
    Q     = Quality Factor
    wr    = Angular resonance frequency [rad/s]
    w     = numpy array of frequencies to calculate the impedance [rad/s]

    Rs, Q and wr may be a float, list, tuple or a numpy array.

    Outputs:
    Zl    = Longitudinal Impedance (real and imaginary) [Ohm]

    Sign convention followed:
        Chao, A., Physics of Collective Beam Instabilities in High Energy
        Accelerators, Wiley 1993.

    Examples:
    >>> w = _np.linspace(-5e9,5e9,1000)
    >>> Rs, Q, wr = 1000, 1, 2*_np.pi*1e9
    >>> Zl = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zl.shape
    (1000,)
    >>> Rs, Q, wr = [1000,2000], [1,10], [2*_np.pi*1e9,2*_np.pi*5e8]
    >>> Zl = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zl.shape
    (1000,)
    >>> Rs, Q, wr = _np.array(Rs), _np.array(Q), _np.array(wr)
    >>> Zl = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zl.shape
    (1000,)

    """
    # I am using broadcasting
    isarr = isinstance(w, _np.ndarray)
    w = _np.array(w, ndmin=1)

    ndim = w.ndim + 1
    Rs = _np.moveaxis(_np.array(Rs, ndmin=ndim, dtype=float), -1, 0)
    Q = _np.moveaxis(_np.array(Q, ndmin=ndim, dtype=float), -1, 0)
    wr = _np.moveaxis(_np.array(wr, ndmin=ndim, dtype=float), -1, 0)
    w = _np.array(w, ndmin=ndim)
    # Zl = w*Rs / (w+1j*Q*(wr - w**2/wr))
    Zl = Rs/(1 + 1j*Q*(wr/w - w/wr))

    Zl = _np.squeeze(Zl.sum(0))
    Zl = _np.array(Zl, ndmin=1)
    if not isarr:
        Zl = Zl[0]
    return Zl


def transverse_resonator(Rs, Q, wr, w):
    """Return the transverse resonator impedance for w.

    Inputs:
    Rs    = Shunt impedance [Ohm]
    Q     = Quality Factor
    wr    = Angular resonance frequency [rad/s]
    w     = numpy array of frequencies to calculate the impedance [rad/s]

    Rs, Q and wr may be a float, list, tuple or a numpy array.

    Outputs:
    Zl    = Longitudinal Impedance (real and imaginary) [Ohm]

    Sign convention followed:
        Chao, A., Physics of Collective Beam Instabilities in High Energy
        Accelerators, Wiley 1993.

    Examples:
    >>> w = _np.linspace(-5e9,5e9,1000)
    >>> Rs, Q, wr = 1000, 1, 2*_np.pi*1e9
    >>> Zt = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zt.shape
    (1000,)
    >>> Rs, Q, wr = [1000,2000], [1,10], [2*_np.pi*1e9,2*_np.pi*5e8]
    >>> Zt = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zt.shape
    (1000,)
    >>> Rs, Q, wr = _np.array(Rs), _np.array(Q), _np.array(wr)
    >>> Zt = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zt.shape
    (1000,)

    """
    isarr = isinstance(w, _np.ndarray)
    w = _np.array(w, ndmin=1)

    ndim = w.ndim + 1
    Rs = _np.moveaxis(_np.array(Rs, ndmin=ndim, dtype=float), -1, 0)
    Q = _np.moveaxis(_np.array(Q, ndmin=ndim, dtype=float), -1, 0)
    wr = _np.moveaxis(_np.array(wr, ndmin=ndim, dtype=float), -1, 0)
    w = _np.array(w, ndmin=ndim)
    Zt = wr*Rs/(w + 1j*Q*(wr - w**2/wr))

    Zt = _np.squeeze(Zt.sum(0))
    Zt = _np.array(Zt, ndmin=1)
    if not isarr:
        Zt = Zt[0]
    return Zt


def wake_longitudinal_resonator(Rs, Q, wr, spos):
    """Return the longitudinal resonator wake-function for spos.

    Inputs:
    Rs    = Shunt impedance [Ohm]
    Q     = Quality Factor
    wr    = Angular resonance frequency [rad/s]
    spos  = numpy array of s positions to calculate the impedance [m]

    Rs, Q and wr may be a float, list, tuple or a numpy array.

    Outputs:
    Wl    = Longitudinal wake-function [Ohm]

    Sign convention followed:
        Chao, A., Physics of Collective Beam Instabilities in High Energy
        Accelerators, Wiley 1993.

    Examples:
    >>> spos = _np.linspace(0, 2, 1000)
    >>> Rs, Q, wr = 1000, 1, 2*_np.pi*1e9
    >>> Wl = wake_longitudinal_resonator(Rs, Q, wr, spos)
    >>> Wl.shape
    (1000,)
    >>> Rs, Q, wr = [1000, 2000], [1, 10], [2*_np.pi*1e9, 2*_np.pi*5e8]
    >>> Wl = wake_longitudinal_resonator(Rs, Q, wr, spos)
    >>> Wl.shape
    (1000,)
    >>> Rs, Q, wr = _np.array(Rs), _np.array(Q), _np.array(wr)
    >>> Wl = wake_longitudinal_resonator(Rs, Q, wr, spos)
    >>> Wl.shape
    (1000,)

    """
    isarr = isinstance(spos, _np.ndarray)
    spos = _np.array(spos, ndmin=1)

    ndim = spos.ndim + 1
    Rs = _np.moveaxis(_np.array(Rs, ndmin=ndim, dtype=float), -1, 0)
    Q = _np.moveaxis(_np.array(Q, ndmin=ndim, dtype=float), -1, 0)
    wr = _np.moveaxis(_np.array(wr, ndmin=ndim, dtype=float), -1, 0)

    alpha = wr / (2*Q)
    wrl = _np.sqrt(wr*wr - alpha*alpha)
    sel = spos >= 0.0
    phase = wrl*spos[sel]/_LSPEED

    wake = 2*alpha * Rs * _np.exp(-alpha*spos[sel]/_LSPEED)
    wake *= _np.cos(phase) - alpha/wrl*_np.sin(phase)

    Wl = _np.zeros(spos.shape)
    Wl[sel] = wake.sum(0).ravel()
    Wl[spos == 0.0] /= 2  # A particle sees half of its wake
    if not isarr:
        Wl = Wl[0]
    return Wl


def wake_transverse_resonator(Rs, Q, wr, spos):
    """Return the Transverse resonator wake-function for spos.

    Inputs:
    Rs    = Shunt impedance [Ohm]
    Q     = Quality Factor
    wr    = Angular resonance frequency [rad/s]
    spos  = numpy array of s positions to calculate the impedance [m]

    Rs, Q and wr may be a float, list, tuple or a numpy array.

    Outputs:
    Wl    = Longitudinal wake-function [Ohm]

    Sign convention followed:
        Chao, A., Physics of Collective Beam Instabilities in High Energy
        Accelerators, Wiley 1993.

    Examples:
    >>> spos = _np.linspace(0, 2, 1000)
    >>> Rs, Q, wr = 1000, 1, 2*_np.pi*1e9
    >>> Wt = wake_longitudinal_resonator(Rs, Q, wr, spos)
    >>> Wt.shape
    (1000,)
    >>> Rs, Q, wr = [1000, 2000], [1, 10], [2*_np.pi*1e9, 2*_np.pi*5e8]
    >>> Wt = wake_longitudinal_resonator(Rs, Q, wr, spos)
    >>> Wt.shape
    (1000,)
    >>> Rs, Q, wr = _np.array(Rs), _np.array(Q), _np.array(wr)
    >>> Wt = longitudinal_resonator(Rs, Q, wr, spos)
    >>> Wt.shape
    (1000,)

    """
    isarr = isinstance(spos, _np.ndarray)
    spos = _np.array(spos, ndmin=1)

    ndim = spos.ndim + 1
    Rs = _np.moveaxis(_np.array(Rs, ndmin=ndim, dtype=float), -1, 0)
    Q = _np.moveaxis(_np.array(Q, ndmin=ndim, dtype=float), -1, 0)
    wr = _np.moveaxis(_np.array(wr, ndmin=ndim, dtype=float), -1, 0)

    alpha = wr / (2*Q)
    wrl = _np.sqrt(wr*wr - alpha*alpha)
    sel = spos > 0.0

    wake = Rs * wr**2 / (Q*wrl) * _np.exp(-alpha*spos[sel]/_LSPEED)
    wake *= _np.sin(wrl*spos[sel]/_LSPEED)

    Wt = _np.zeros(spos.shape)
    Wt[sel] = wake.sum(0).flatten()
    if not isarr:
        Wt = Wt[0]
    return Wt
