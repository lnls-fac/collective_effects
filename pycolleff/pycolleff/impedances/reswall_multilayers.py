"""."""

import time as _time
from functools import partial as _partial

import numpy as _np
from scipy.special import kve as _kve, ive as _ive, kv as _kv, iv as _iv
import scipy.integrate as _scyint
import mpmath as _mpmath
from mpmath import besseli as _miv, besselk as _mkv, mpc as _mpc, mpf as _mpf

import mathphys as _mp

_LSPEED = _mp.constants.light_speed
_ep0 = _mp.constants.vacuum_permitticity
_Z0 = _mp.constants.vacuum_impedance
_E0 = _mp.constants.electron_rest_energy * _mp.units.joule_2_eV

RW_default_w = _np.logspace(-1, 13, 1000)


def yokoya_factors(plane='ll'):
    if plane == 'll':
        return 1
    elif plane == 'dy':
        return _np.pi**2/12
    elif plane in {'qy', 'dx'}:
        return _np.pi**2/24
    elif plane == 'qx':
        return -_np.pi**2/24
    else:
        raise Exception(
            "Plane not identified. Possible options: "
            "'ll', 'dy', 'dx', 'qy', 'qx' ")


def prepare_input_epr_mur(w, epb, mub, ange, angm, sigmadc, tau):
    epr = _np.zeros((len(epb), len(w)), dtype=_np.complex)
    mur = _np.zeros((len(epb), len(w)), dtype=_np.complex)
    for j in range(len(epb)):
        epr[j] = epb[j]*(1 - 1j*_np.sign(w)*_np.tan(ange[j]))
        mur[j] = mub[j]*(1 - 1j*_np.sign(w)*_np.tan(angm[j]))

        epr[j] += sigmadc[j]/(1 + 1j*w*tau[j]) / (1j*w*_ep0)
    return epr, mur


def resistive_multilayer_round_pipe(w, epr, mur, b, L, E):

    def Mtil(m, epr, mur, bet, nu, b):
        for i in range(len(b)):  # length(b) = # de camadas - 1
            x = nu[i+1] * b[i]
            y = nu[i] * b[i]
            Mt = _np.zeros((4, 4, w.shape[0]), dtype=_np.complex)

            if i < len(b)-1:
                D = _np.zeros((4, 4, nu.shape[1]), dtype=_np.complex)
                z = nu[i+1] * b[i+1]
                if not (z.real < 0).any():
                    ind = z.real < 60

                    A = _iv(m, z[ind])
                    B = _kv(m, z[ind])
                    C = _iv(m, x[ind])
                    E = _kv(m, x[ind])

                    D[0, 0] = 1
                    D[2, 2] = 1
                    D[1, 1, ind] = - B*C / (A*E)
                    D[3, 3, ind] = - B*C / (A*E)
                    D[1, 1, ~ind] = - _np.exp(-2*(z[~ind]-x[~ind]))
                    D[3, 3, ~ind] = - _np.exp(-2*(z[~ind]-x[~ind]))

            Mt[0, 0] = -nu[i+1]**2*b[i]/epr[i+1]*(
                epr[i+1]/nu[i+1]*(-_kve(m-1, x)/_kve(m, x) - m/x) -
                epr[i]/nu[i]*(_ive(m-1, y)/_ive(m, y) - m/y))
            Mt[0, 1] = -nu[i+1]**2*b[i]/epr[i+1]*(
                epr[i+1]/nu[i+1]*(-_kve(m-1, x)/_kve(m, x) - m/x) -
                epr[i]/nu[i]*(-_kve(m-1, y)/_kve(m, y) - m/y))
            Mt[1, 0] = -nu[i+1]**2*b[i]/epr[i+1]*(
                epr[i+1]/nu[i+1] * (_ive(m-1, x)/_ive(m, x) - m/x) -
                epr[i]/nu[i] * (_ive(m-1, y)/_ive(m, y) - m/y))
            Mt[1, 1] = -nu[i+1]**2*b[i]/epr[i+1]*(
                epr[i+1]/nu[i+1] * (_ive(m-1, x)/_ive(m, x) - m/x) -
                epr[i]/nu[i] * (-_kve(m-1, y)/_kve(m, y) - m/y))

            Mt[0, 2] = (nu[i+1]**2/nu[i]**2 - 1)*m/(bet*epr[i+1])
            Mt[0, 3] = Mt[0, 2]
            Mt[1, 2] = Mt[0, 2]
            Mt[1, 3] = Mt[0, 2]
            Mt[2, 2] = -nu[i+1]**2*b[i]/mur[i+1]*(
                mur[i+1]/nu[i+1] * (-_kve(m-1, x)/_kve(m, x) - m/x) -
                mur[i]/nu[i] * (_ive(m-1, y)/_ive(m, y) - m/y))
            Mt[2, 3] = -nu[i+1]**2*b[i]/mur[i+1]*(
                mur[i+1]/nu[i+1] * (-_kve(m-1, x)/_kve(m, x) - m/x) -
                mur[i]/nu[i] * (-_kve(m-1, y)/_kve(m, y) - m/y))
            Mt[3, 2] = -nu[i+1]**2*b[i]/mur[i+1]*(
                mur[i+1]/nu[i+1]*(_ive(m-1, x)/_ive(m, x) - m/x) -
                mur[i]/nu[i] * (_ive(m-1, y)/_ive(m, y) - m/y))
            Mt[3, 3] = -nu[i+1]**2*b[i]/mur[i+1]*(
                mur[i+1]/nu[i+1]*(_ive(m-1, x)/_ive(m, x) - m/x) -
                mur[i]/nu[i] * (-_kve(m-1, y)/_kve(m, y) - m/y))
            Mt[2, 0] = (nu[i+1]**2/nu[i]**2 - 1) * m/(bet*mur[i+1])
            Mt[2, 1] = Mt[2, 0]
            Mt[3, 0] = Mt[2, 0]
            Mt[3, 1] = Mt[2, 0]

            if len(b) == 1:
                M = Mt
            else:
                if not i:
                    M = _np.einsum('ijk,jlk->ilk', D, Mt)
                elif i < len(b)-1:
                    M = _np.einsum('ijk,jlk->ilk', Mt, M)
                    M = _np.einsum('ijk,jlk->ilk', D, M)
                else:
                    M = _np.einsum('ijk,jlk->ilk', Mt, M)
        return M

    def alphaTM(m, epr, mur, bet, nu, b):
        M = Mtil(m, epr, mur, bet, nu, b)

        B = M[0, 1]*M[2, 2] - M[2, 1]*M[0, 2]
        B /= M[0, 0]*M[2, 2] - M[0, 2]*M[2, 0]
        alphaTM = _kv(m, nu[0]*b[0])/_iv(m, nu[0]*b[0]) * B
        return alphaTM

    ####################
    gamma = E/_E0
    beta = _np.sqrt(1 - 1/gamma**2)
    nu = _np.abs(w/_LSPEED)*_np.sqrt(1 - beta**2*epr*mur)

    Zl = 1j*L*w / (2*_np.pi*_ep0 * (beta*_LSPEED)**2*gamma**2)
    Zl *= alphaTM(0, epr, mur, beta, nu, b)

    Zv = 1j*L*w**2 / (4*_np.pi*_ep0*_LSPEED**2*(beta*_LSPEED)*gamma**4)
    Zv *= alphaTM(1, epr, mur, beta, nu, b)

    # The code cant handle w = 0;
    ind0 = (w == 0).nonzero()[0]
    if ind0:
        Zl[ind0] = (Zl[ind0+1] + Zl[ind0-1]).real/2
        Zv[ind0] = 0 + 1j*Zv[ind0+1].imag

    Zh = Zv.copy()

    return Zl.conj(), Zv.conj(), Zh.conj()


def resistive_multilayer_round_pipe_multiprecision(
        w, epr, mur, b, L, E, prec=70, print_progress=True):

    def Mtil(m, epr, mur, bet, nu, b):
        def produto(A, B):
            C = _mpmath.matrix(A.rows, B.cols)
            for i in range(C.rows):
                for j in range(C.cols):
                    for k in range(A.cols):
                        C[i, j] += A[i, k]*B[k, j]
            return C

        for i in range(len(b)):  # length(b) = # de camadas - 1
            x = nu[i+1] * b[i]
            y = nu[i] * b[i]
            Mt = _mpmath.matrix(4, 4)

            if i < len(b)-1:
                D = _mpmath.matrix(4, 4)
                z = nu[i+1] * b[i+1]
                if not (z.real < 0):
                    D[0, 0] = 1
                    D[2, 2] = 1
                    if z.real < 60:
                        A = _miv(m, z)
                        B = _mkv(m, z)
                        C = _miv(m, x)
                        E = _mkv(m, x)
                        D[1, 1] = - B*C / (A*E)
                        D[3, 3] = - B*C / (A*E)
                    else:
                        D[1, 1] = - _mpmath.exp(-2*(z - x))
                        D[3, 3] = - _mpmath.exp(-2*(z - x))

            Mt[0, 0] = -nu[i+1]**2*b[i]/epr[i+1]*(
                    epr[i+1]/nu[i+1]*(-_mkv(m-1, x)/_mkv(m, x) - m/x) -
                    epr[i]/nu[i]*(_miv(m-1, y)/_miv(m, y) - m/y))
            Mt[0, 1] = -nu[i+1]**2*b[i]/epr[i+1]*(
                    epr[i+1]/nu[i+1]*(-_mkv(m-1, x)/_mkv(m, x) - m/x) -
                    epr[i]/nu[i]*(-_mkv(m-1, y)/_mkv(m, y) - m/y))
            Mt[1, 0] = -nu[i+1]**2*b[i]/epr[i+1]*(
                    epr[i+1]/nu[i+1] * (_miv(m-1, x)/_miv(m, x) - m/x) -
                    epr[i]/nu[i] * (_miv(m-1, y)/_miv(m, y) - m/y))
            Mt[1, 1] = -nu[i+1]**2*b[i]/epr[i+1]*(
                    epr[i+1]/nu[i+1] * (_miv(m-1, x)/_miv(m, x) - m/x) -
                    epr[i]/nu[i] * (-_mkv(m-1, y)/_mkv(m, y) - m/y))

            Mt[0, 2] = (nu[i+1]**2/nu[i]**2 - 1)*m/(bet*epr[i+1])
            Mt[0, 3] = Mt[0, 2]
            Mt[1, 2] = Mt[0, 2]
            Mt[1, 3] = Mt[0, 2]
            Mt[2, 2] = -nu[i+1]**2*b[i]/mur[i+1]*(
                    mur[i+1]/nu[i+1] * (-_mkv(m-1, x)/_mkv(m, x) - m/x) -
                    mur[i]/nu[i] * (_miv(m-1, y)/_miv(m, y) - m/y))
            Mt[2, 3] = -nu[i+1]**2*b[i]/mur[i+1]*(
                    mur[i+1]/nu[i+1] * (-_mkv(m-1, x)/_mkv(m, x) - m/x) -
                    mur[i]/nu[i] * (-_mkv(m-1, y)/_mkv(m, y) - m/y))
            Mt[3, 2] = -nu[i+1]**2*b[i]/mur[i+1]*(
                    mur[i+1]/nu[i+1]*(_miv(m-1, x)/_miv(m, x) - m/x) -
                    mur[i]/nu[i] * (_miv(m-1, y)/_miv(m, y) - m/y))
            Mt[3, 3] = -nu[i+1]**2*b[i]/mur[i+1]*(
                    mur[i+1]/nu[i+1]*(_miv(m-1, x)/_miv(m, x) - m/x) -
                    mur[i]/nu[i] * (-_mkv(m-1, y)/_mkv(m, y) - m/y))
            Mt[2, 0] = (nu[i+1]**2/nu[i]**2 - 1) * m/(bet*mur[i+1])
            Mt[2, 1] = Mt[2, 0]
            Mt[3, 0] = Mt[2, 0]
            Mt[3, 1] = Mt[2, 0]

            if len(b) == 1:
                M = Mt
            else:
                if not i:
                    M = produto(D, Mt)
                elif i < len(b)-1:
                    M = produto(Mt, M)
                    M = produto(D, M)
                else:
                    M = produto(Mt, M)
        return M

    def alphaTM(m, epr, mur, bet, nu, b):
        M = Mtil(m, epr, mur, bet, nu, b)

        B = M[0, 1]*M[2, 2] - M[2, 1]*M[0, 2]
        B /= M[0, 0]*M[2, 2] - M[0, 2]*M[2, 0]
        alphaTM = _mkv(m, nu[0]*b[0])/_miv(m, nu[0]*b[0]) * B
        return alphaTM

    ####################
    gamma = E/_E0
    beta = _np.sqrt(1 - 1/gamma**2)
    nu = _np.ones(
        (epr.shape[0], 1))*abs(w/_LSPEED)*_np.sqrt(1 - beta**2*epr*mur)

    _mpmath.mp.dps = prec

    Zl = []
    Zv = []
    beta_ = _mpf(beta)
    b_ = [_mpc(b_) for b_ in b]
    for i, w_ in enumerate(w):
        t0_ = _time.time()
        epr_ = [_mpc(epr[j, i]) for j in range(epr.shape[0])]
        mur_ = [_mpc(mur[j, i]) for j in range(mur.shape[0])]
        nu_ = [_mpc(nu[j, i]) for j in range(nu.shape[0])]
        Zl_ = complex(alphaTM(0, epr_, mur_, beta_, nu_, b_))
        Zl_ *= 1j*L*w_ / (2*_np.pi*_ep0 * (beta*_LSPEED)**2*gamma**2)
        Zl.append(Zl_)

        Zv_ = complex(alphaTM(1, epr_, mur_, beta_, nu_, b_))
        Zv_ *= 1j*L*w_**2 / (4*_np.pi*_ep0*_LSPEED**2*(beta*_LSPEED)*gamma**4)
        Zv.append(Zv_)
        if print_progress:
            print(
                f'{i:04d}/{len(w):04d} -> freq = {w_/2/_np.pi:10.2g} '
                f' (ET: {_time.time()-t0_:.2f} s)')
    Zl = _np.array(Zl)
    Zv = _np.array(Zv)

    # The code can't handle w = 0;
    ind0 = (w == 0).nonzero()[0]
    if ind0:
        Zl[ind0] = (Zl[ind0+1] + Zl[ind0-1]).real/2
        Zv[ind0] = 0 + 1j*Zv[ind0+1].imag

    Zh = Zv.copy()

    return Zl.conj(), Zv.conj(), Zh.conj()


def _calc_matrix_m_prec(kx, beta, epr, mur, nu, b):
    for p, b_ in enumerate(b):  # length(b) = # de camadas - 1
        Mt = _mpmath.matrix(4, 4)

        kyp = _mpmath.sqrt(kx*kx + nu[p]*nu[p])
        kyp1 = _mpmath.sqrt(kx*kx + nu[p+1]*nu[p+1])

        exp_nn = _mpmath.exp(-(kyp + kyp1)*b_)
        exp_pn = _mpmath.exp((kyp - kyp1)*b_)
        try:
            exp_pp = 1/exp_nn
            exp_np = 1/exp_pn
        except ZeroDivisionError:
            return Mt

        nu_ri2 = nu[p+1]/nu[p]
        nu_ri2 *= nu_ri2
        ky_r = kyp/kyp1
        mu_r = mur[p]/mur[p+1]
        ep_r = epr[p]/epr[p+1]

        pre = ky_r*nu_ri2/2
        fac = pre*ep_r
        Mt[0, 0] = (0.5+fac)*exp_pn
        Mt[0, 1] = (0.5-fac)*exp_nn
        Mt[1, 0] = (0.5-fac)*exp_pp
        Mt[1, 1] = (0.5+fac)*exp_np

        fac = kx*(nu_ri2 - 1)/(2*beta*kyp1*epr[p+1])
        Mt[0, 2] = -fac*exp_pn
        Mt[0, 3] = -fac*exp_nn
        Mt[1, 2] = fac*exp_pp
        Mt[1, 3] = fac*exp_np

        fac = epr[p+1]/mur[p+1]
        Mt[2, 0] = fac*Mt[0, 2]
        Mt[2, 1] = fac*Mt[0, 3]
        Mt[3, 0] = fac*Mt[1, 2]
        Mt[3, 1] = fac*Mt[1, 3]

        fac = pre*mu_r
        Mt[2, 2] = (0.5+fac)*exp_pn
        Mt[2, 3] = (0.5-fac)*exp_nn
        Mt[3, 2] = (0.5-fac)*exp_pp
        Mt[3, 3] = (0.5+fac)*exp_np

        if not p:
            M = Mt
        else:
            M = Mt*M
    return M


def _eta_xi_funcs_prec(v, beta, epru, muru, nuu, bu, eprd, murd, nud, bd):
    epru = [_mpmath.mpc(epi) for epi in epru]
    muru = [_mpmath.mpc(mui) for mui in muru]
    nuu = [_mpmath.mpc(nui) for nui in nuu]
    bu = [_mpmath.mpf(bi) for bi in bu]
    eprd = [_mpmath.mpc(epi) for epi in eprd]
    murd = [_mpmath.mpc(mui) for mui in murd]
    nud = [_mpmath.mpc(nui) for nui in nud]
    bd = [_mpmath.mpf(bi) for bi in bd]
    v = _mpmath.mpf(v)
    beta = _mpmath.mpf(beta)

    matmp = _calc_matrix_m_prec(v, beta, epru, muru, nuu, bu)
    matmn = _calc_matrix_m_prec(v, beta, eprd, murd, nud, bd)
    matp = _mpmath.zeros(4, 4)
    matp[:2, :] = matmp[::2, :]
    matp[2:, :] = matmn[1::2, :]
    try:
        pinv = matp**-1
    except ZeroDivisionError:
        print(v, 'zero')
        return [_mpmath.mpc(0) for _ in range(4)]

    xi1 = pinv[0, 0]*matmp[0, 1] + pinv[0, 1] * matmp[2, 1]
    xi2 = pinv[1, 0]*matmp[0, 1] + pinv[1, 1] * matmp[2, 1]
    eta1 = pinv[0, 2]*matmn[1, 0] + pinv[0, 3] * matmn[3, 0]
    eta2 = pinv[1, 2]*matmn[1, 0] + pinv[1, 3] * matmn[3, 0]
    return xi1, xi2, eta1, eta2


def _integrand(u, kovgamma=1, is_in_t=False, **kwrgs):
    if is_in_t:
        t, u = u, (1-u)/u
    u = _mpmath.mpf(u)
    v = kovgamma * _mpmath.sinh(u)
    xi1, xi2, eta1, eta2 = _eta_xi_funcs_prec(v, **kwrgs)

    # General expression:
    #   itgmn = xi1 + (-1)**m*eta1 + (-1)**n*xi2 + (-1)**(m+n)*eta2
    #   itgmn *= _mpmath.cosh(m*u)*_mpmath.cosh(n*u)

    itg00 = xi1 + eta1 + xi2 + eta2
    itg02 = itg00*_mpmath.cosh(2*u)

    itg11 = xi1 - eta1 - xi2 + eta2
    itg11 *= _mpmath.cosh(u)**2

    itg00 = complex(itg00)
    itg11 = complex(itg11)
    itg02 = complex(itg02)
    itgs = _np.array([
        itg00.real, itg00.imag,
        itg11.real, itg11.imag,
        itg02.real, itg02.imag])
    if is_in_t:
        itgs /= t*t
    return itgs


def calc_integrand(u, w, E, *args, is_in_t=False, prec=70):
    _mpmath.mp.dps = prec
    gamma = E/_E0
    beta = _np.sqrt(1 - 1/gamma**2)
    k = w/(beta*_LSPEED)

    epr_up, mur_up, b_up, epr_dn, mur_dn, b_dn = args
    kovgamma = k/gamma
    nu_up = _np.abs(k)*_np.sqrt(1 - beta**2*epr_up*mur_up)
    nu_dn = _np.abs(k)*_np.sqrt(1 - beta**2*epr_dn*mur_dn)

    kwrgs = dict(
        kovgamma=kovgamma, beta=beta, is_in_t=is_in_t,
        epru=epr_up, muru=mur_up, nuu=nu_up, bu=b_up,
        eprd=epr_dn, murd=mur_dn, nud=nu_dn, bd=b_dn)

    alpha00 = _np.zeros(u.shape, dtype=complex)
    alpha11 = _np.zeros(u.shape, dtype=complex)
    alpha02 = _np.zeros(u.shape, dtype=complex)
    for i, ui in enumerate(u):
        alphas = _integrand(ui, **kwrgs)
        alpha00[i] = alphas[0] + 1j*alphas[1]
        alpha11[i] = alphas[2] + 1j*alphas[3]
        alpha02[i] = alphas[4] + 1j*alphas[5]

    return alpha00, alpha11, alpha02


def calc_alphas(k, gamma, beta, *args, is_in_t=True):
    epr_up, mur_up, b_up, epr_dn, mur_dn, b_dn = args
    kovgamma = k/gamma
    nu_up = _np.abs(k)*_np.sqrt(1 - beta**2*epr_up*mur_up)
    nu_dn = _np.abs(k)*_np.sqrt(1 - beta**2*epr_dn*mur_dn)

    supl = _np.inf
    if is_in_t:
        supl = 1
    alpha00 = _np.zeros(k.shape, dtype=complex)
    alpha11 = _np.zeros(k.shape, dtype=complex)
    alpha02 = _np.zeros(k.shape, dtype=complex)
    for i, kovgi in enumerate(kovgamma):
        t0_ = _time.time()
        kwrgs = dict(
            is_in_t=is_in_t, kovgamma=kovgi, beta=beta,
            epru=epr_up[:, i], muru=mur_up[:, i], nuu=nu_up[:, i], bu=b_up,
            eprd=epr_dn[:, i], murd=mur_dn[:, i], nud=nu_dn[:, i], bd=b_dn)

        alphas, err, info = _scyint.quad_vec(
            _partial(_integrand, **kwrgs),
            0, supl, epsabs=1e-6, norm='max', epsrel=1e-4, full_output=True,
            workers=-1)
        alpha00[i], alpha11[i], alpha02[i] = alphas[::2] + 1j*alphas[1::2]
        print(
            f'{i:04d}/{len(k):04d} -> '
            f'freq = {k[i]*beta*_LSPEED/2/_np.pi:8.2g} Hz,   '
            f'converged = {"no " if info.status else "yes":s},   '
            f'n_evals = {info.neval:04d} '
            f' (ET: {_time.time()-t0_:.2f} s)')

    return alpha00, alpha11, alpha02


def resistive_multilayer_flat_pipe(
        w, L, E, epr_up, mur_up, b_up, epr_dn, mur_dn, b_dn,
        prec=70, is_in_t=True):
    _mpmath.mp.dps = prec

    gamma = E/_E0
    beta = _np.sqrt(1 - 1/gamma**2)
    k = w/(beta*_LSPEED)

    args = epr_up, mur_up, b_up, epr_dn, mur_dn, b_dn
    alpha00, alpha11, alpha02 = calc_alphas(
        k, gamma, beta, *args, is_in_t=is_in_t)

    fac = 1j*k * _Z0 * L / (2*_np.pi * beta*gamma**2)
    Zll = fac * alpha00
    # Zmy = fac / gamma * alpha01

    fac *= k / gamma**2
    Zdy = fac * alpha11
    fac /= 2
    Zqy = fac * (alpha00 + alpha02)
    Zqx = fac * (alpha00 - alpha02)
    Zdx = -Zqx.copy()

    return Zll.conj(), Zdx.conj(), Zdy.conj(), Zqx.conj(), Zqy.conj()
