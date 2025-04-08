"""."""

import time as _time
from functools import partial as _partial

import mathphys as _mp
import mpmath as _mpmath
import numpy as _np
import scipy.integrate as _scyint
from mpmath import besseli as _miv, besselk as _mkv, cosh as _cosh, \
    exp as _exp, eye as _eye, fabs as _fabs, matrix as _matrix, mpc as _mpc, \
    mpf as _mpf, sinh as _sinh, sqrt as _sqrt, zeros as _zeros
from scipy.special import iv as _iv, ive as _ive, kv as _kv, kve as _kve

_LSPEED = _mp.constants.light_speed
_ep0 = _mp.constants.vacuum_permitticity
_Z0 = _mp.constants.vacuum_impedance  # Z0 = 1/ep0/c = u0*c
_E0 = _mp.constants.electron_rest_energy * _mp.units.joule_2_eV

_pos = _np.logspace(-7, 0, 500)
_neg = _np.logspace(-7, -3, 100)
DEFAULT_SPOS_RW = _np.r_[-_np.flipud(_neg), _pos]


def get_default_reswall_w(radius=12e-3, energy=3e9):
    """Return a suggested vector of angular frequencies.

    The maximum suggested angular frequency is based on the discussion at the
    end of the page 47 of ref. [1].

    Refs:
        [1] Mounet, N. (2012). The LHC Transverse Coupled-Bunch Instability.
            École Polytechinique Fédérale de Lausanne.

    Args:
        radius (float, optional): inner radius or half gap in [m].
            Defaults to 12e-3.
        energy (float, optional): energy of the beam in [eV]. Defaults to 3e9.

    Returns:
        (numpy.ndarray, (401, )): angular frequencies [rad/s].

    """
    gamma = energy / _E0
    beta = _np.sqrt(1 - 1 / (gamma * gamma))
    cut = 10 * beta * _LSPEED * gamma / radius
    return _np.logspace(-1, _np.log10(cut), 401) * 2 * _np.pi


def get_impedance_for_negative_w(Z, w=None, impedance_type="ll"):
    """Get the impedance for negative frequencies from positive one.

    Uses the symmetry relations of impedances to get impedance at negative
    frequencies from values at positive frequencies.

    The longitudinal impedance must satisfy:
        Zl(-w) = conj(Zl(w))

    and the transverse impedance must satisfy:
        Zt(-w) = -conj(Zt(w))

    where conj(.) means the complex conjugate.

    Args:
        Z (numpy.ndarray, (N, )): Impedance for positive frequencies.
        w (numpy.ndarray, (N, ), optional): angular frequencies. If given, a
            vector containing negative and positive frequencies will also be
            returned. It is always assumed, without checking, that the
            frequencies are all positive. Defaults to None.
        impedance_type (str, optional): type of impedance. May assume values
            {'ll', 'long', 'trans', 't', 'dx', 'dy', 'qx', 'qy' }.
            Defaults to 'll'.

    Returns:
        Z (numpy.ndarray, (2*N, )): Impedance for negative and positive
            frequencies, respectively. If the input vector was given for
            frequencies in ascending order, this vector will maintain this
            property.
        w (numpy.ndarray, (2*N, )): Negative and positive angular frequencies.
            Only returned when input `w` is given.  If the input vector was in
            ascending order, this vector will maintain this property.

    """
    sig = 1 if impedance_type.startswith("l") else -1
    Z = _np.r_[sig * _np.flipud(Z.conj()), Z]
    if w is not None:
        w = _np.r_[-_np.flipud(w), w]
        return Z, w
    return Z


def yokoya_factors(impedance_type="ll", geometry="flat"):
    """Return the Yokoya factors for all planes in relation to 'round'.

    Refs:
        [1] Yokoya, K. (1993). Resistive Wall Impedance of Beam Pipes of
            General Cross Section. Particle Accelerators, 41, 221-248.

    Args:
        impedance_type (str, optional): type of geometry. Assume values in
        {'flat', 'square', 'round'}. Defaults to 'flat'.
        geometry (str, optional): type of impedance. Assume values in
        {'ll', 'dy', 'dx', 'qy', 'qx'}. Defaults to 'll'.

    Raises:
        ValueError: when geometry or impedance_type cannot be interpreted.

    Returns:
        float: the Yokoya factor for the selected inputs.

    """
    if geometry.lower().startswith("flat"):
        base = _np.pi**2 / 12
        yok = dict(ll=1, dy=base, dx=base / 2, qx=-base / 2, qy=base / 2)
    elif geometry.lower().startswith("square"):
        base = 0.85  # approximately from Figure 8.
        yok = dict(ll=1, dy=base, dx=base, qx=0, qy=0)
    elif geometry.lower().startswith("round"):
        yok = dict(ll=1, dy=1, dx=1, qx=0, qy=0)
    else:
        raise ValueError(
            "geometry not identified. Possible options: "
            "'flat', 'square', 'round'"
        )
    if impedance_type not in yok:
        raise ValueError(
            "impedance_type not identified. Possible options: "
            "'ll', 'dy', 'dx', 'qy', 'qx' "
        )
    return yok[impedance_type]


def prepare_inputs_epr_mur(w, epb, mub, ange, angm, sigmadc, tau):
    """Calculate frequency dependent electromagnetic properties.

    The complex relative electric permittivity (`epr`) and magnetic
    permeability (`mur`) are calculated from the materials properties
    described below for all the M layers. Please, consult eqs. 3.7 and 3.8 of
    ref. [1] and ref. [2] and the discussion that follows for more information
    on this model.

    Refs:
        [1] Mounet, N., & Métral, E. (2010). Electromagnetic fields and beam
            coupling impedances in a multilayer flat chamber. Geneva,
            Switzerland.
        [2] Mounet, N., & Métral, E. (2009). Electromagnetic field created by
            a macroparticle in an infinitel long and axisymmetric multilayer
            beam pipe. Geneva, Switzerland.

    Args:
        w (numpy.ndarray, (N, )): angular frequency where to calculate the
            properties in [rad/s].
        epb (numpy.ndarray, (M, )): relative dielectric permittivity/constant.
        mub (numpy.ndarray, (M, )): real part of the relative complex
            permeability.
        ange (numpy.ndarray, (M, )): dielectric loss angle (used to calculate
            the loss tangent).
        angm (numpy.ndarray, (M, )): magnetic loss angle (used to calculate
            the loss tangent).
        sigmadc (numpy.ndarray, (M, )): Conductivity at DC frequency in [S].
        tau (numpy.ndarray, (M, )): Material relaxation time according to
            Drude model.

    Returns:
        epr (numpy.ndarray, (M, N)): relative complex permittivity.
        mur (numpy.ndarray, (M, N)): relative complex permeability.

    """
    epr = _np.zeros((len(epb), len(w)), dtype=complex)
    mur = _np.zeros((len(epb), len(w)), dtype=complex)
    for j in range(len(epb)):
        epr[j] = epb[j] * (1 - 1j * _np.sign(w) * _np.tan(ange[j]))
        mur[j] = mub[j] * (1 - 1j * _np.sign(w) * _np.tan(angm[j]))
        epr[j] += sigmadc[j] / (1 + 1j * w * tau[j]) / (1j * w * _ep0)
    return epr, mur


def multilayer_round_chamber(
    w,
    L,
    E,
    epr,
    mur,
    b,
    precision=70,
    wmax_arb_prec=0.0,
    arb_prec_incl_long=False,
    print_progress=True,
):
    """Calculate multilayer resistive wall impedance for round cross-section.

    This method implements the calculation of the impedance theory developed
    in ref. [1] for a round chamber with an arbitrary number of layers of
    arbitrary linear materials (see figure 1 of ref. [1]), characterized by
    the complex relative eletric permittivity (`epr`) and the complex relative
    permeability (`mur`). The first layer is the region where the beam is
    located and is generally modelled as vacuum. The other layers inner
    surface are at a distance `b[i]` from the center of the chamber and the
    last layer extends to infinity.

    The values and frequency dependency of `epr` or `mur` is arbitrary,
    however a possible model can be created with the method
    `prepare_inputs_epr_mur`.

    This implementation only calculates the transverse dominant terms for the
    impedance as function of the radial coordinate. So it will return only the
    monopole longitudinal impedance, `Zl`, and the horizontal and vertical
    dipole impedances, `Zdx` and `Zdy`, respectively. The detuning impedances
    are not returned because for a ultra-relativistic beam, according to
    ref. [1], they vanish. However, since they are proportional to the
    monopole longitudinal impedance, one can easly obtain them.
    Please look at eqs 8.42-8.46 of ref. [1] for more details.

    IMPORTANT NOTE: Even though we follow ref. [1], the impedance definition
        used in this work is different from the one of ref.[3], which we
        follow in other functions and modules of this package. For this
        reason, we return the complex conjugate of the impedances, so they
        would agree with the impedance definition of ref.[3].

    Since the calculation of these formulas may be subjected for numerical
    noise when simple double-precision arithmetic is used, whe also give the
    option for the user to use arbitraty precision arithmetic througth the
    variable `wmax_arb_prec`. Since the numerical noises tend to happen for
    lower frequencies, this input variable controls above which frequency the
    regular double-precision is used. Besides that, the variable
    `arb_prec_incl_long` determines whether the arbitrary precision arithmetic
    should be used for the longitudinal impedance as well. This is useful
    because the longitunal impedance calculation suffer much less from this
    numerical noise.

    Please, have in mind that the arbitrary precision calculation is much
    slower. For this reason we strongly recommend the user to run this method
    wish `wmax_arb_prec=0` first to check whether or not arbitraty precision
    is needed.

    Finally, the wall impedance is non-negligible in a large frequency range,
    going from mHz to ~10^15 Hz, where it starts decaying to zero very fast.
    For this reason we recommend the user to use a log scale impedance
    sampling to evaluate the impedance. Experience with different number of
    layers showed so far that ~500 points 0.1Hz to 3e15Hz is good enough to
    sample most impedances.

    Refs:
        [1] Mounet, N., & Métral, E. (2009). Electromagnetic field created by
            a macroparticle in an infinitel long and axisymmetric multilayer
            beam pipe. Geneva, Switzerland.
        [2] Mounet, N. (2012). The LHC Transverse Coupled-Bunch Instability.
            École Polytechinique Fédérale de Lausanne.
        [3] Chao, A. W. (1993). Physics of Collective Beam Instabilities in
            High Energy Accelerators (1st ed.). New York: John Wiley & Sons.

    Args:
        w (numpy.ndarray, (N, )): Angular frequencies where to calculate
            impedance. Must be given in units of [rad/s].
        L (float): Length of the vacuum chamber in [m]. The theory in which
            this implementation was based consider an infinitely long
            structure, so this value is only a multiplicative factor.
        E (float): Energy of the beam in [eV].
        epr (numpy.ndarray, (M, N)): Frequency dependent relative complex
            permittivity for each layer. See `prepare_inputs_epr_mur` to
            understand how to model this parameter.
        mur (numpy.ndarray, (M, N)): Frequency dependent relative complex
            permeability for each layer. See `prepare_inputs_epr_mur` to
            understand how to model this parameter.
        b (numpy.ndarray, (M-1, )): Internal radius of each layer, starting
            from the second layer. The first layer is the one where the beam
            is, which is generally considered to be vacuum.
        precision (int, optional): Number of decimal places to consider in
            arbitrary-precision calculations. Defaults to 70.
        wmax_arb_prec (float, optional): For angular frequencies above this
            value [rad/s], arbitrary-precision calculation will be performed
            to evaluate impedances. Defaults to 0.0.
        arb_prec_incl_long (bool, optional): Whether or not to also calculate
            the longitudinal impedance using arbitrary-precision. Generally
            not needed because the longitudinal impedance expressions are not
            so demanding. Defaults to False.
        print_progress (bool, optional): Whether or not to print progress of
            the calculation. Defaults to True.

    Returns:
        Zl (numpy.ndarray, (N, )): Longitudinal impedance in [Ohm].
        Zdx (numpy.ndarray, (N, )): Horizontal dipolar impedance in [Ohm/m].
        Zdy (numpy.ndarray, (N, )): Vertical dipolar impedance in [Ohm/m].

    """
    _mpmath.mp.dps = precision
    Zl = _np.zeros(w.shape, dtype=complex)
    Zv = _np.zeros(w.shape, dtype=complex)

    regp = _np.abs(w) > wmax_arb_prec
    arbp = ~regp

    if arb_prec_incl_long:
        Zl[regp] = _round_chamber_reg_prec(
            w[regp], E, epr[:, regp], mur[:, regp], b, m=0
        )
        Zl[arbp] = _round_chamber_arb_prec(
            w[arbp],
            E,
            epr[:, arbp],
            mur[:, arbp],
            b,
            m=0,
            print_=print_progress,
        )
    else:
        Zl = _round_chamber_reg_prec(w, E, epr, mur, b, m=0)

    Zv[regp] = _round_chamber_reg_prec(
        w[regp], E, epr[:, regp], mur[:, regp], b, m=1
    )
    Zv[arbp] = _round_chamber_arb_prec(
        w[arbp], E, epr[:, arbp], mur[:, arbp], b, m=1, print_=print_progress
    )

    gamma = E / _E0
    beta = _np.sqrt(1 - 1 / gamma**2)
    k = w / (beta * _LSPEED)
    fac = 1j * k * _Z0 * L / (2 * _np.pi * beta * gamma**2)
    Zl *= fac
    Zv *= fac * k / gamma**2 / 2

    Zh = Zv.copy()
    # Need to take the conjugate because our convention of impedance follows
    # ref. [3], which is different from the ref. [1].
    return Zl.conj(), Zh.conj(), Zv.conj()


def _round_chamber_reg_prec(w, E, epr, mur, b, m=0):
    gamma = E / _E0
    beta = _np.sqrt(1 - 1 / gamma**2)
    k = w / (beta * _LSPEED)
    nu = _np.abs(k) * _np.sqrt(1 - beta**2 * epr * mur)
    return _round_alpha_tm(m, epr, mur, beta, nu, b)


def _round_chamber_arb_prec(w, E, epr, mur, b, m=0, print_=True):
    gamma = _mpf(E / _E0)
    beta = _sqrt(1 - 1 / gamma**2)

    Z = []
    b_ = [_mpc(b_) for b_ in b]
    for i, w_ in enumerate(w):
        t0_ = _time.time()
        k = w_ / (beta * _LSPEED)
        epr_ = [_mpc(x) for x in epr[:, i]]
        mur_ = [_mpc(x) for x in mur[:, i]]
        nu_ = [
            _fabs(k) * _sqrt(1 - beta**2 * e * m) for e, m in zip(epr_, mur_)
        ]

        Z_ = _round_alpha_tm_arb_prec(m, epr_, mur_, beta, nu_, b_)
        Z.append(complex(Z_))
        if print_:
            print(
                f"{i:04d}/{len(w):04d} -> freq = {w_/2/_np.pi:10.2g} "
                f" (ET: {_time.time()-t0_:.2f} s)"
            )

    return _np.array(Z)


def _round_alpha_tm(m, epr, mur, bet, nu, b):
    for i in range(len(b)):
        x = nu[i + 1] * b[i]
        y = nu[i] * b[i]
        Mt = _np.zeros((4, 4, nu.shape[1]), dtype=complex)
        D = _np.zeros((4, 4, nu.shape[1]), dtype=complex)

        kmx = _kve(m, x)
        kmy = _kve(m, y)
        imx = _ive(m, x)
        imy = _ive(m, y)
        km1x = _kve(m - 1, x)
        km1y = _kve(m - 1, y)
        im1x = _ive(m - 1, x)
        im1y = _ive(m - 1, y)

        if i < len(b) - 1:
            z = nu[i + 1] * b[i + 1]
            if not (z.real < 0).any():
                ind = z.real < 60

                A = _iv(m, z[ind])
                B = _kv(m, z[ind])
                C = _iv(m, x[ind])
                E = _kv(m, x[ind])

                D[0, 0] = 1
                D[2, 2] = 1
                D[1, 1, ind] = -B * C / (A * E)
                D[3, 3, ind] = D[1, 1, ind]
                D[1, 1, ~ind] = -_np.exp(-2 * (z[~ind] - x[~ind]))
                D[3, 3, ~ind] = D[1, 1, ~ind]
            else:
                print("z.real < 0")

        k_m_x = km1x / kmx + m / x
        k_m_y = km1y / kmy + m / y
        i_m_x = im1x / imx - m / x
        i_m_y = im1y / imy - m / y

        nu2_ep_b = nu[i + 1] ** 2 / epr[i + 1] * b[i]
        ep_nu_1 = epr[i + 1] / nu[i + 1]
        ep_nu = epr[i] / nu[i]
        Mt[0, 0] = nu2_ep_b * (ep_nu_1 * k_m_x + ep_nu * i_m_y)
        Mt[0, 1] = nu2_ep_b * (ep_nu_1 * k_m_x - ep_nu * k_m_y)
        Mt[1, 0] = -nu2_ep_b * (ep_nu_1 * i_m_x - ep_nu * i_m_y)
        Mt[1, 1] = -nu2_ep_b * (ep_nu_1 * i_m_x + ep_nu * k_m_y)

        nu2_mu_b = nu[i + 1] ** 2 / mur[i + 1] * b[i]
        mu_nu_1 = mur[i + 1] / nu[i + 1]
        mu_nu = mur[i] / nu[i]
        Mt[2, 2] = nu2_mu_b * (mu_nu_1 * k_m_x + mu_nu * i_m_y)
        Mt[2, 3] = nu2_mu_b * (mu_nu_1 * k_m_x - mu_nu * k_m_y)
        Mt[3, 2] = -nu2_mu_b * (mu_nu_1 * i_m_x - mu_nu * i_m_y)
        Mt[3, 3] = -nu2_mu_b * (mu_nu_1 * i_m_x + mu_nu * k_m_y)

        nur = nu[i + 1] / nu[i]
        nur *= nur
        Mt[0, 2] = (nur - 1) * m / (bet * epr[i + 1])
        Mt[0, 3] = Mt[0, 2]
        Mt[1, 2] = Mt[0, 2]
        Mt[1, 3] = Mt[0, 2]
        Mt[2, 0] = (nur - 1) * m / (bet * mur[i + 1])
        Mt[2, 1] = Mt[2, 0]
        Mt[3, 0] = Mt[2, 0]
        Mt[3, 1] = Mt[2, 0]

        if len(b) == 1:
            M = Mt
        else:
            if not i:
                M = _np.einsum("ijk,jlk->ilk", D, Mt)
            elif i < len(b) - 1:
                M = _np.einsum("ijk,jlk->ilk", Mt, M)
                M = _np.einsum("ijk,jlk->ilk", D, M)
            else:
                M = _np.einsum("ijk,jlk->ilk", Mt, M)

    B = M[0, 1] * M[2, 2] - M[2, 1] * M[0, 2]
    B /= M[0, 0] * M[2, 2] - M[0, 2] * M[2, 0]
    alphaTM = _kv(m, nu[0] * b[0]) / _iv(m, nu[0] * b[0]) * B
    return alphaTM


def _round_alpha_tm_arb_prec(m, epr, mur, bet, nu, b):
    for i in range(len(b)):
        x = nu[i + 1] * b[i]
        y = nu[i] * b[i]
        Mt = _matrix(4, 4)
        D = _eye(4)

        kmx = _mkv(m, x)
        kmy = _mkv(m, y)
        imx = _miv(m, x)
        imy = _miv(m, y)
        km1x = _mkv(m - 1, x)
        km1y = _mkv(m - 1, y)
        im1x = _miv(m - 1, x)
        im1y = _miv(m - 1, y)

        if i < len(b) - 1:
            z = nu[i + 1] * b[i + 1]
            if not (z.real < 0):
                D[1, 1] = -_mkv(m, z) * imx / (_miv(m, z) * kmx)
                D[3, 3] = D[1, 1]
            else:
                print("z.real < 0")

        k_m_x = km1x / kmx + m / x
        k_m_y = km1y / kmy + m / y
        i_m_x = im1x / imx - m / x
        i_m_y = im1y / imy - m / y

        nu2_ep_b = nu[i + 1] ** 2 / epr[i + 1] * b[i]
        ep_nu_1 = epr[i + 1] / nu[i + 1]
        ep_nu = epr[i] / nu[i]
        Mt[0, 0] = nu2_ep_b * (ep_nu_1 * k_m_x + ep_nu * i_m_y)
        Mt[0, 1] = nu2_ep_b * (ep_nu_1 * k_m_x - ep_nu * k_m_y)
        Mt[1, 0] = -nu2_ep_b * (ep_nu_1 * i_m_x - ep_nu * i_m_y)
        Mt[1, 1] = -nu2_ep_b * (ep_nu_1 * i_m_x + ep_nu * k_m_y)

        nu2_mu_b = nu[i + 1] ** 2 / mur[i + 1] * b[i]
        mu_nu_1 = mur[i + 1] / nu[i + 1]
        mu_nu = mur[i] / nu[i]
        Mt[2, 2] = nu2_mu_b * (mu_nu_1 * k_m_x + mu_nu * i_m_y)
        Mt[2, 3] = nu2_mu_b * (mu_nu_1 * k_m_x - mu_nu * k_m_y)
        Mt[3, 2] = -nu2_mu_b * (mu_nu_1 * i_m_x - mu_nu * i_m_y)
        Mt[3, 3] = -nu2_mu_b * (mu_nu_1 * i_m_x + mu_nu * k_m_y)

        nur = nu[i + 1] / nu[i]
        nur *= nur
        Mt[0, 2] = (nur - 1) * m / (bet * epr[i + 1])
        Mt[0, 3] = Mt[0, 2]
        Mt[1, 2] = Mt[0, 2]
        Mt[1, 3] = Mt[0, 2]
        Mt[2, 0] = (nur - 1) * m / (bet * mur[i + 1])
        Mt[2, 1] = Mt[2, 0]
        Mt[3, 0] = Mt[2, 0]
        Mt[3, 1] = Mt[2, 0]

        if not i:
            M = D * Mt
        elif i < len(b) - 1:
            M = D * Mt * M
        else:
            M = Mt * M

    B = M[0, 1] * M[2, 2] - M[2, 1] * M[0, 2]
    B /= M[0, 0] * M[2, 2] - M[0, 2] * M[2, 0]
    alphaTM = _mkv(m, nu[0] * b[0]) / _miv(m, nu[0] * b[0]) * B
    return alphaTM


def multilayer_flat_chamber(
    w,
    L,
    E,
    epr_up,
    mur_up,
    b_up,
    epr_dn=None,
    mur_dn=None,
    b_dn=None,
    precision=70,
    print_progress=True,
):
    """Calculate multilayer resistive wall impedance for flat cross-section.

    This method implements the calculation of the impedance theory developed
    in ref. [1] for a flat chamber with an arbitrary number (M) of layers of
    arbitrary linear materials (see figure 1 of ref. [1]), characterized by
    the complex relative eletric permittivity (`epr_up` and `epr_dn`) and the
    complex relative permeability (`mur_up` and `mur_dn`). The first upper and
    lower layer is the region where the beam is located and is generally
    modelled as vacuum. The other upper and lower layers inner surface are at
    a distance `b_up[i]` and `b_dn[i]` from the center of the chamber and the
    last layer extends to infinity.

    The values and frequency dependency of `epr` or `mur` is arbitrary,
    however a possible model can be created with the method
    `prepare_inputs_epr_mur`.

    If any of the bottom layer parameters `epr_dn`, `mur_dn` or `b_dn` is
    None, them top-bottom symmetry is assumed.

    IMPORTANTE NOTES:
    - Please note that `b_dn` must have negative values.
    - Even though we follow ref. [1], the impedance definition
        used in this work is different from the one of ref.[3], which we
        follow in other functions and modules of this package. For this
        reason, we return the complex conjugate of the impedances, so they
        would agree with the impedance definition of ref.[3].

    This implementation only calculates the transverse dominant terms for the
    impedance as function of the radial coordinate. So it will return only the
    monopole longitudinal impedance, `Zl`, the horizontal and vertical
    dipole impedances, `Zdx` and `Zdy`, and the horizontal and vertical
    detuning impedances, `Zqx` and `Zqy`. Please look at eqs 9.23-9.28 of
    ref. [1] for more details. Please, note that we do not calculate the
    monopolar vertical impedance given by equation 9.24, which is non-zero
    when there is assymetry between the upper and lower layers.

    Since the calculation of these formulas may be subjected for numerical
    noise when simple double-precision arithmetic is used, we employ an
    arbitrary precision library to perform all calculations, which makes this
    method rather slow.

    Finally, the wall impedance is non-negligible in a large frequency range,
    going from mHz to ~10^15 Hz, where it starts decaying to zero very fast.
    For this reason we recommend the user to use a log scale impedance
    sampling to evaluate the impedance. Experience with different number of
    layers showed so far that ~500 points 0.1Hz to 3e15Hz is good enough to
    sample most impedances.

    Refs:
        [1] Mounet, N., & Métral, E. (2010). Electromagnetic fields and beam
            coupling impedances in a multilayer flat chamber. Geneva,
            Switzerland.
        [2] Mounet, N. (2012). The LHC Transverse Coupled-Bunch Instability.
            École Polytechinique Fédérale de Lausanne.
        [3] Chao, A. W. (1993). Physics of Collective Beam Instabilities in
            High Energy Accelerators (1st ed.). New York: John Wiley & Sons.

    Args:
        w (numpy.ndarray, (N, )): Angular frequencies where to calculate
            impedance. Must be given in units of [rad/s].
        L (float): Length of the vacuum chamber in [m]. The theory in which
            this implementation was based consider an infinitely long
            structure, so this value is only a multiplicative factor.
        E (float): Energy of the beam in [eV].
        epr_up (numpy.ndarray, (M, N)): Frequency dependent relative complex
            permittivity for each upper layer. See `prepare_inputs_epr_mur` to
            understand how to model this parameter.
        mur_up (numpy.ndarray, (M, N)): Frequency dependent relative complex
            permeability for each upper layer. See `prepare_inputs_epr_mur` to
            understand how to model this parameter.
        b_up (numpy.ndarray, (M-1, )): Vertical position of internal surface
            of each upper layer, starting from the second layer. The first
            layer is where the beam is, which is generally considered to be
            vacuum.
        epr_dn (numpy.ndarray, (M, N), optional): Frequency dependent relative
            complex permittivity for each lower layer.
            See `prepare_inputs_epr_mur` to understand how to model this
            parameter. If is None, then top-bottom symmetry is assume.
            Defaults to None.
        mur_dn (numpy.ndarray, (M, N), optional): Frequency dependent relative
            complex permeability for each lower layer.
            See `prepare_inputs_epr_mur` to understand how to model this
            parameter. If is None, then top-bottom symmetry is assume.
            Defaults to None.
        b_dn (numpy.ndarray, (M-1, ), optional): Vertical position of internal
            surface of each lower layer, starting from the second layer. Note
            that these values must be negative! The first layer is where the
            beam is, which is generally considered to be vacuum. If is None,
            then top-bottom symmetry is assume. Defaults to None.
        precision (int, optional): Number of decimal places to consider in
            arbitrary-precision calculations. Defaults to 70.
        print_progress (bool, optional): Whether or not to print progress of
            the calculation. Defaults to True.

    Returns:
        Zl (numpy.ndarray, (N, )): Longitudinal impedance in [Ohm].
        Zdx (numpy.ndarray, (N, )): Horizontal dipolar impedance in [Ohm/m].
        Zdy (numpy.ndarray, (N, )): Vertical dipolar impedance in [Ohm/m].
        Zqx (numpy.ndarray, (N, )): Horizontal detuning impedance in [Ohm/m].
        Zqy (numpy.ndarray, (N, )): Vertical detuning impedance in [Ohm/m].

    Raises:
        ZeroDivisionError: May happen when `precision` is not large enough.

    """
    _mpmath.mp.dps = precision

    # Whether or not to make a variable transformation in integration:
    is_in_t = True

    symm = False
    if epr_dn is None or mur_dn is None or b_dn is None:
        symm = True
        epr_dn, mur_dn, b_dn = epr_up, mur_up, -b_up

    gamma = _mpf(E / _E0)
    beta = _sqrt(1 - 1 / gamma**2)
    k = [wi / (beta * _LSPEED) for wi in w]

    bu = [_mpf(x) for x in b_up]
    bd = [_mpf(x) for x in b_dn]

    epru, muru, eprd, murd = [], [], [], []
    for i in range(epr_up.shape[1]):
        epru.append([_mpc(x) for x in epr_up[:, i]])
        muru.append([_mpc(x) for x in mur_up[:, i]])
        eprd.append([_mpc(x) for x in epr_dn[:, i]])
        murd.append([_mpc(x) for x in mur_dn[:, i]])

    args = epru, muru, bu, eprd, murd, bd
    alpha00, alpha11, alpha02 = _flat_calc_alphas(
        k,
        gamma,
        beta,
        *args,
        is_in_t=is_in_t,
        symm=symm,
        print_=print_progress,
    )

    # Starting here everything is double-precision:
    gamma = E / _E0
    beta = _np.sqrt(1 - 1 / gamma**2)
    k = w / (beta * _LSPEED)

    fac = 1j * k * _Z0 * L / (2 * _np.pi * beta * gamma**2)
    Zll = fac * alpha00
    # This term is the monopolar vertical impedance from eq. 9.24 of ref. [1]:
    # Zmy = fac / gamma * alpha01

    fac *= k / gamma**2
    Zdy = fac * alpha11
    fac /= 2
    Zqy = fac * (alpha00 + alpha02)
    Zqx = fac * (alpha00 - alpha02)
    Zdx = -Zqx.copy()

    # Need to take the conjugate because our convention of impedance follows
    # ref. [3], which is different from the ref. [1].
    return Zll.conj(), Zdx.conj(), Zdy.conj(), Zqx.conj(), Zqy.conj()


def _flat_calc_alphas(
    k, gamma, beta, *args, is_in_t=True, symm=False, print_=True
):
    epru, muru, bu, eprd, murd, bd = args

    supl = _np.inf
    if is_in_t:
        supl = 1
    alpha00 = _np.zeros(len(k), dtype=complex)
    alpha11 = _np.zeros(len(k), dtype=complex)
    alpha02 = _np.zeros(len(k), dtype=complex)
    for i, ki in enumerate(k):
        t0_ = _time.time()

        nuu, nud = [], []
        for eu, mu, ed, md in zip(epru[i], muru[i], eprd[i], murd[i]):
            nuu.append(_fabs(ki) * _sqrt(1 - beta**2 * eu * mu))
            nud.append(_fabs(ki) * _sqrt(1 - beta**2 * ed * md))

        kovgamma = ki / gamma
        kwrgs = dict(
            is_in_t=is_in_t,
            symm=symm,
            kovgamma=kovgamma,
            beta=beta,
            epru=epru[i],
            muru=muru[i],
            nuu=nuu,
            bu=bu,
            eprd=eprd[i],
            murd=murd[i],
            nud=nud,
            bd=bd,
            prec=_mpmath.mp.dps,  # needed for child processes in quad_vec
        )

        # The output of this integration is in standard precision:
        alphas, err, info = _scyint.quad_vec(
            _partial(_flat_integrand_arb_prec, **kwrgs),
            0,
            supl,
            epsabs=1e-6,
            norm="max",
            epsrel=1e-4,
            full_output=True,
            workers=-1,
        )
        alpha00[i], alpha11[i], alpha02[i] = alphas[::2] + 1j * alphas[1::2]
        if not print_:
            continue
        print(
            f'{i:04d}/{len(k):04d} -> '
            f'freq = {float(ki*beta)*_LSPEED/2/_np.pi:8.2g} Hz,   '
            f'converged = {"no " if info.status else "yes":s},   '
            f'n_evals = {info.neval:04d} '
            f' (ET: {_time.time()-t0_:.2f} s)'
        )

    return alpha00, alpha11, alpha02


# NOTE: This function exists only for debuging purposes
def _debug_flat_calc_integrand(
    u, w, E, *args, is_in_t=False, symm=False, prec=70
):
    _mpmath.mp.dps = prec
    gamma = _mpf(E / _E0)
    beta = _sqrt(1 - 1 / gamma**2)
    k = w / (beta * _LSPEED)
    epru, muru, bu, eprd, murd, bd = args

    epru = [_mpc(epi) for epi in epru]
    muru = [_mpc(mui) for mui in muru]
    bu = [_mpf(bi) for bi in bu]
    eprd = [_mpc(epi) for epi in eprd]
    murd = [_mpc(mui) for mui in murd]
    bd = [_mpf(bi) for bi in bd]
    beta = _mpf(beta)

    kovgamma = k / gamma
    nuu, nud = [], []
    for i, _ in enumerate(epru):
        nuu.append(_fabs(k) * _sqrt(1 - beta**2 * epru[i] * muru[i]))
        nud.append(_fabs(k) * _sqrt(1 - beta**2 * eprd[i] * murd[i]))

    kwrgs = dict(
        kovgamma=kovgamma,
        beta=beta,
        is_in_t=is_in_t,
        symm=symm,
        epru=epru,
        muru=muru,
        nuu=nuu,
        bu=bu,
        eprd=eprd,
        murd=murd,
        nud=nud,
        bd=bd,
        prec=_mpmath.mp.dps,
    )

    alpha00 = _np.zeros(u.shape, dtype=complex)
    alpha11 = _np.zeros(u.shape, dtype=complex)
    alpha02 = _np.zeros(u.shape, dtype=complex)
    for i, ui in enumerate(u):
        alphas = _flat_integrand_arb_prec(ui, **kwrgs)
        alpha00[i] = alphas[0] + 1j * alphas[1]
        alpha11[i] = alphas[2] + 1j * alphas[3]
        alpha02[i] = alphas[4] + 1j * alphas[5]
    return alpha00, alpha11, alpha02


def _flat_integrand_arb_prec(u, kovgamma=1, is_in_t=False, **kwrgs):
    _mpmath.mp.dps = kwrgs.pop('prec')  # Needed for Mac and Windows systems.

    if is_in_t:
        t, u = u, (1 - u) / u
    u = _mpf(u)
    v = kovgamma * _sinh(u)
    xi1, xi2, eta1, eta2 = _flat_eta_xi_funs_arb_prec(v, **kwrgs)

    # General expression:
    #   itgmn = xi1 + (-1)**m*eta1 + (-1)**n*xi2 + (-1)**(m+n)*eta2
    #   itgmn *= _cosh(m*u)*_cosh(n*u)

    itg00 = xi1 + eta1 + xi2 + eta2
    itg02 = itg00 * _cosh(2 * u)

    itg11 = xi1 - eta1 - xi2 + eta2
    itg11 *= _cosh(u) ** 2

    itg00 = complex(itg00)
    itg11 = complex(itg11)
    itg02 = complex(itg02)
    itgs = _np.array(
        [
            itg00.real,
            itg00.imag,
            itg11.real,
            itg11.imag,
            itg02.real,
            itg02.imag,
        ]
    )
    if is_in_t:
        itgs /= t * t
    return itgs


def _flat_eta_xi_funs_arb_prec(
    v, beta, epru, muru, nuu, bu, eprd, murd, nud, bd, symm=False
):
    Mps = _flat_calc_matrix_arb_prec(v, beta, epru, muru, nuu, bu)
    if symm:
        Mns = _flat_symmetric_matrix(Mps)
    else:
        Mns = _flat_calc_matrix_arb_prec(v, beta, eprd, murd, nud, bd)

    Mp = _flat_mult_mats(Mps)
    Mn = _flat_mult_mats(Mns)
    matp = _zeros(4, 4)
    matp[:2, :] = Mp[::2, :]
    matp[2:, :] = Mn[1::2, :]

    # This inversion may raise ZeroDivisionError when precision is
    # not large enough:
    pinv = matp**-1

    xi1 = pinv[0, 0] * Mp[0, 1] + pinv[0, 1] * Mp[2, 1]
    xi2 = pinv[1, 0] * Mp[0, 1] + pinv[1, 1] * Mp[2, 1]
    eta1 = pinv[0, 2] * Mn[1, 0] + pinv[0, 3] * Mn[3, 0]
    eta2 = pinv[1, 2] * Mn[1, 0] + pinv[1, 3] * Mn[3, 0]
    return xi1, xi2, eta1, eta2


def _flat_calc_matrix_arb_prec(kx, beta, epr, mur, nu, b):
    Mts = []
    for p, b_ in enumerate(b):
        Mt = _matrix(4, 4)

        kyp = _sqrt(kx * kx + nu[p] * nu[p])
        kyp1 = _sqrt(kx * kx + nu[p + 1] * nu[p + 1])

        # This next four lines may raise ZeroDivisionError when precision is
        # not large enough:
        exp_nn = _exp(-(kyp + kyp1) * b_)
        exp_pn = _exp((kyp - kyp1) * b_)
        exp_pp = 1 / exp_nn
        exp_np = 1 / exp_pn

        nu_ri2 = nu[p + 1] / nu[p]
        nu_ri2 *= nu_ri2
        ky_r = kyp / kyp1
        mu_r = mur[p] / mur[p + 1]
        ep_r = epr[p] / epr[p + 1]

        pre = ky_r * nu_ri2 / 2
        fac = pre * ep_r
        Mt[0, 0] = (0.5 + fac) * exp_pn
        Mt[0, 1] = (0.5 - fac) * exp_nn
        Mt[1, 0] = (0.5 - fac) * exp_pp
        Mt[1, 1] = (0.5 + fac) * exp_np

        fac = kx * (nu_ri2 - 1) / (2 * beta * kyp1 * epr[p + 1])
        Mt[0, 2] = -fac * exp_pn
        Mt[0, 3] = -fac * exp_nn
        Mt[1, 2] = fac * exp_pp
        Mt[1, 3] = fac * exp_np

        fac = epr[p + 1] / mur[p + 1]
        Mt[2, 0] = fac * Mt[0, 2]
        Mt[2, 1] = fac * Mt[0, 3]
        Mt[3, 0] = fac * Mt[1, 2]
        Mt[3, 1] = fac * Mt[1, 3]

        fac = pre * mu_r
        Mt[2, 2] = (0.5 + fac) * exp_pn
        Mt[2, 3] = (0.5 - fac) * exp_nn
        Mt[3, 2] = (0.5 - fac) * exp_pp
        Mt[3, 3] = (0.5 + fac) * exp_np

        Mts.append(Mt)
    return Mts


def _flat_symmetric_matrix(Mps):
    Mns = []
    for Mp in Mps:
        Mn = _zeros(4, 4)
        Mn[0, 0] = Mp[1, 1]
        Mn[0, 1] = Mp[1, 0]
        Mn[1, 0] = Mp[0, 1]
        Mn[1, 1] = Mp[0, 0]

        Mn[0, 2] = -Mp[1, 3]
        Mn[0, 3] = -Mp[1, 2]
        Mn[1, 2] = -Mp[0, 3]
        Mn[1, 3] = -Mp[0, 2]

        Mn[2, 0] = -Mp[3, 1]
        Mn[2, 1] = -Mp[3, 0]
        Mn[3, 0] = -Mp[2, 1]
        Mn[3, 1] = -Mp[2, 0]

        Mn[2, 2] = Mp[3, 3]
        Mn[2, 3] = Mp[3, 2]
        Mn[3, 2] = Mp[2, 3]
        Mn[3, 3] = Mp[2, 2]
        Mns.append(Mn)
    return Mns


def _flat_mult_mats(Ms):
    M = Ms[0]
    for mm in Ms[1:]:
        M = mm * M
    return M
