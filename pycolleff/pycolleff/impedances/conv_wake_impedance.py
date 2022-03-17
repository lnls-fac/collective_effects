"""Implement some impedance models and Element and Budget classes."""

import numpy as _np
from scipy.interpolate import PchipInterpolator as _Pchip, \
    CubicSpline as _Spline

import mathphys as _mp

_Z0 = _mp.constants.vacuum_impedance
_LSPEED = _mp.constants.light_speed


def from_wake_to_impedance(z, wake, bunlen, cutoff=2):
    dt = (z[-1]-z[0]) / (z.shape[0]-1) / _LSPEED
    VHat = _np.fft.fft(wake) * dt
    w = 2 * _np.pi * _np.fft.fftfreq(len(z), d=dt)
    VHat = _np.fft.fftshift(VHat)   # shift the negative frequencies
    w = _np.fft.fftshift(w)         # to the center of the spectrum
    VHat *= _np.exp(-1j * w * z[0] / _LSPEED)

    wmax = cutoff / bunlen * _LSPEED
    indcs = _np.abs(w) <= wmax
    # Deconvolve the Transform with a gaussian bunch:
    Jwlist = _np.exp(-(w * bunlen / _LSPEED)**2 / 2)
    Z = VHat[indcs] / Jwlist[indcs]
    w = w[indcs]
    return w, Z.conj()


def from_impedance_to_wake(
        s, w, Z, plane='long', interp_type='spline', ret_interp=False, beta=1):
    """Calculate wake from given sampled impedance.

    Implementation follows the derivation presented in section 1.6.3 of:

    [1] Mounet, N. (2012). The LHC Transverse Coupled-Bunch Instability.
        École Polytechinique Fédérale de Lausanne.

    This implementation assumes the impedance (Z) is calculated for positive
    values of angular frequencies (w) only. Simmetry properties of the
    impedances are used to account for the Fourier integral in negative
    frequencies.
    It also assumes the sampling of the impedance is well performed such that
    the cubic Hermite interpolation is a good approximation for the impedance.
    The function gives an option of returning the interpolation object used in
    the calculations (`ret_interp=True`), so that the user can check
    afterwards whether this assumption is True or not.

    Args:
        s (numpy.ndarray, (N, )): positions where to evaluate wake [m].
        w (numpy.ndarray, (M, )): angular frequencies where impedance is
            defined [rad/s]. I assume this vector only contains posivite
            frequencies.
        Z (numpy.ndarray, (M, )): Longitudinal or transverse impedance
            [Ohm or Ohm/m].
        plane (str, optional): type of the impedance. May assume values:
            'long' or 'trans'. Defaults to 'long'.
        interp_type (str, optional): Defines algorithm to be used for the
            calculations of dZ/dw for the integration. May assume values:
            'spline' or 'monotone' (see ref.[1] and scipy.interpolate help
            for more infos) . Defaults to 'spline'.
        ret_interp (bool, optional): Whether or not to also return the
            interpolation object used to calculate the derivatives of the
            impedance. Defaults to False.
        beta (float, optional): The particle's velocity in units of the speed
            of light. Detaults to 1.

    Raises:
        ValueError: when any angular frequency is negative.

    Returns:
        wake (numpy.ndarray, (N, )): wake at `s` positions behind the source.
        interp (CubicSpline or PchipInterpolator from scipy.interpolate,
            optional): interpolation object used to calculate derivatives of
            the impedance. Only returned when `ret_interp=True`.

    """
    if w.min() < 0:
        raise ValueError('All angular frequencies must be positive.')

    # I neet to take the conjugate of the impedance here, because the
    # convention adopted by ref[1] is different than the one we consider here.
    Z = _np.array(Z.conj(), dtype=_np.complex256)

    if interp_type.lower().startswith('spline'):
        interp = _Spline(w, Z, bc_type='not-a-knot')
    else:
        interp = _Pchip(w, Z)
    dZdw = interp.derivative()(w)

    t = s[:, None]/(_LSPEED*beta)
    deltai = _np.diff(w)[None, :]

    arg = deltai * t
    phip, psip = _integral_funcs(arg)
    phin, psin = _integral_funcs(-arg)
    expdi = _np.exp(1j*arg)

    fi = Z[None, :-1]
    fip1 = Z[None, 1:]
    wi = w[None, :-1]
    di = dZdw[None, :-1]
    dip1 = dZdw[None, 1:]
    expi = _np.exp(1j*wi*t)

    # This numerical implementation has a grouping different from the one
    # presented in eq. 1.239, because this form is more robust against
    # numerical noise. This idea was taken from the IW2D implementation.
    du1 = deltai * (dip1*psip - di*expdi*psin)
    du2 = expdi*fi*phin + fip1*phip
    integ = expi*deltai*(du1+du2)
    integ = integ.sum(axis=1)

    # consider the integral up to infinity given by eq. E.133 (or 1.232)
    t = t.ravel()
    integ += _np.exp(1j*w[-1]*t) * 1j*Z[-1] / t

    if plane.lower().startswith('long'):
        integ = _np.array(integ.real, dtype=float)
    else:
        integ = _np.array(integ.imag, dtype=float)
    integ /= _np.pi
    if ret_interp:
        return integ, interp
    return integ


def _integral_funcs(x):
    """Calculate Fourier integrals of Hermite interpolation polynomials.

    Derivation follows ref:

    [1] Mounet, N. (2012). The LHC Transverse Coupled-Bunch Instability.
        École Polytechinique Fédérale de Lausanne.

    Args:
        x (numpy.ndarray): values where to calculate the integral

    Returns:
        phi: Phi function of equation 1.239 of ref. [1]
        psi: Psi function of equation 1.239 of ref. [1]

    """
    # I decided to use the exact formula for indices where x > 0.01
    # and Taylor series otherwise
    absx = _np.abs(x)
    iser = absx < 0.01
    ifor = ~iser
    x_ser = x[iser]  # series
    x_for = x[ifor]  # formula

    phi = _np.zeros(x.shape, dtype=_np.complex256)
    psi = _np.zeros(x.shape, dtype=_np.complex256)

    # The exact formula is taken from ref. [1], eqs. E.142 and E.143
    exp = _np.exp(1j*x_for)
    x2 = x_for*x_for
    x3 = x2*x_for
    x4 = x3*x_for
    phi[ifor] = -1j*exp/x_for - 6j*(exp+1)/x3 + 12*(exp-1)/x4
    psi[ifor] = exp/x2 + 2j*(2*exp+1)/x3 - 6*(exp-1)/x4

    # The series implementation follows Appendix E.3.3 of ref. [1]:
    # I noticed the accuracy estimative of eq. E.146 for Psi is more extrict
    # than the one of E.147 for Phi. Which means we could use only the
    # the first to control the convergence. Besides I dropped the term
    # (p + 5) in the denominator, relaxing it, so it would be easy to calculate
    # iteratively.
    # Lastly I took the logarithm of the errors to avoid divergence:
    absx = absx[iser].max()  # the slowest x is the one with maximum abs(x)
    log_conv = _np.log(2/5) + absx
    # For small x, the imaginary part of Phi and Psi goes linearly with |x|,
    # so we define the error as a fraction of the smallest x:
    log_err = _np.log(_np.min(absx)/1000)
    for i in range(1000):
        log_conv += _np.log(absx/(i+1))
        if not i:
            term = _np.ones(x_ser.shape, dtype=_np.complex256)
        else:
            term *= 1j*x_ser/i
        add = term / (i+3) / (i+4)
        psi[iser] -= add
        phi[iser] += add * (i+6)
        if log_conv < log_err:
            break
    return phi, psi
