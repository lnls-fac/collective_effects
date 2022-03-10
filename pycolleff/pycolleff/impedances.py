"""Implement some impedance models and Element and Budget classes."""

import os as _os
import time as _time

import numpy as _np
import matplotlib.pyplot as _plt
import scipy.signal as _scysig
from scipy.special import kve as _kve, ive as _ive, kv as _kv, iv as _iv, \
    hyp1f1 as _hyp1f1, gamma as _gamma, airy as _airy, airye as _airye
from scipy.interpolate import PchipInterpolator as _Pchip, \
    CubicSpline as _Spline
import mpmath as _mpmath
from mpmath import besseli as _miv, besselk as _mkv, mpc as _mpc, mpf as _mpf

import lnls as _lnls
import mathphys as _mp

_c = _mp.constants.light_speed
_mu0 = _mp.constants.vacuum_permeability
_ep0 = _mp.constants.vacuum_permitticity
_Z0 = _mp.constants.vacuum_impedance
E0 = _mp.constants.electron_rest_energy * _mp.units.joule_2_eV

_IMPS = {'Zll', 'Zdy', 'Zdx', 'Zqy', 'Zqx'}
_WAKES = {'Wll', 'Wdy', 'Wdx', 'Wqy', 'Wqx'}
_TITLE = {
    'Zll': 'Longitudinal Impedance',
    'Zdy': 'Driving Vertical Impedance',
    'Zdx': 'Driving Horizontal Impedance',
    'Zqy': 'Detuning Vertical Impedance',
    'Zqx': 'Detuning Horizontal Impedance'}
_FACTOR = {
    'Zll': 1e-3, 'Zdy': 1e-3, 'Zdx': 1e-3, 'Zqy': 1e-3, 'Zqx': 1e-3,
    'Wll': 1e-3, 'Wdy': 1e-6, 'Wdx': 1e-6, 'Wqy': 1e-6, 'Wqx': 1e-6}
_BETA = {
    'Zll': lambda x: 1,       'Wll': lambda x: 1,
    'Zdy': lambda x: x.betay, 'Wdy': lambda x: x.betay,
    'Zdx': lambda x: x.betax, 'Wdx': lambda x: x.betax,
    'Zqy': lambda x: x.betay, 'Wqy': lambda x: x.betay,
    'Zqx': lambda x: x.betax, 'Wqx': lambda x: x.betax,
    }

RW_default_w = _np.logspace(1, 13, 3000)
RW_default_w = _np.sort(_np.hstack([-RW_default_w, RW_default_w]))


def _plotlog(x, y, color=None, label=None, ax=None, linewidth=1.5):
    if ax is None:
        ax = _plt.gca()

    if any(y > 0):
        ax.loglog(x, y, color=color, label=label, linewidth=linewidth)
        ax.loglog(x, -y, '--', color=color, linewidth=linewidth)
    else:
        ax.loglog(x, -y, '--', color=color, linewidth=linewidth)
        ax.loglog(x, y, color=color, label=label, linewidth=linewidth)


def _prepare_props(props):
    if isinstance(props, str):
        if props.lower() == 'all':
            props = sorted(list(_IMPS))
        elif props in _IMPS:
            props = [props]
        else:
            raise AttributeError(props, ' not supported for plot.')
    elif isinstance(props, (list, tuple)):
        wrong = set(props) - _IMPS
        if wrong:
            raise AttributeError(wrong, ' not supported for plot.')
    else:
        raise TypeError("Type '"+type(props)+"' not supported for 'props'.")
    return props


class Element:

    _YLABEL = {
        'Zll': r'$Z_l [k\Omega]$',
        'Zdy': r'$Z_y^D [k\Omega/m]$',
        'Zdx': r'$Z_x^D [k\Omega/m]$',
        'Zqy': r'$Z_y^Q [k\Omega/m]$',
        'Zqx': r'$Z_x^Q [k\Omega/m]$'}

    def __init__(
            self, name=None, path=None, betax=None, betay=None, quantity=None):
        self.name = name or 'element'
        if path is not None:
            path = _os.path.abspath(path)
        self.path = path or _os.path.abspath('.')
        self.quantity = quantity or 0  # this field must only be used in Budget
        self.betax = betax or 0.0  # this field shall only be used in Budget
        self.betay = betay or 0.0  # this field shall only be used in Budget
        self.w = _np.array([], dtype=float)
        self.Zll = _np.array([], dtype=complex)
        self.Zdy = _np.array([], dtype=complex)
        self.Zdx = _np.array([], dtype=complex)
        self.Zqy = _np.array([], dtype=complex)
        self.Zqx = _np.array([], dtype=complex)
        self.s = _np.array([], dtype=float)
        self.Wll = _np.array([], dtype=float)
        self.Wdy = _np.array([], dtype=float)
        self.Wdx = _np.array([], dtype=float)
        self.Wqy = _np.array([], dtype=float)
        self.Wqx = _np.array([], dtype=float)

    def copy(self):
        other = Element(
            name=self.name, path=self.path, betax=self.betax,
            betay=self.betay, quantity=self.quantity)
        other.w = self.w.copy()
        other.Zll = self.Zll.copy()
        other.Zdy = self.Zdy.copy()
        other.Zdx = self.Zdx.copy()
        other.Zqy = self.Zqy.copy()
        other.Zqx = self.Zqx.copy()
        other.s = self.s.copy()
        other.Wll = self.Wll.copy()
        other.Wdy = self.Wdy.copy()
        other.Wdx = self.Wdx.copy()
        other.Wqy = self.Wqy.copy()
        other.Wqx = self.Wqx.copy()
        return other

    def save(self):
        name = self.name.replace(' ', '_').lower()
        _lnls.utils.save_pickle(
            _os.path.sep.join([self.path, name]), element=self)

    def load(self):
        name = self.name.replace(' ', '_').lower()
        data = _lnls.utils.load_pickle(_os.path.sep.join([self.path, name]))
        return data['element']

    def plot(
            self, props='all', logscale=True, show=True, save=False, name='',
            figsize=(8, 4)):

        if name:
            name = '_'+name
        props = _prepare_props(props)

        for prop in props:
            Imp2 = getattr(self, prop)
            if Imp2 is None or len(Imp2) == 0:
                continue
            _plt.figure(figsize=figsize)
            Imp = Imp2*_FACTOR[prop]
            w = self.w
            if logscale:
                _plotlog(w, Imp.real, color='b', label='Real')
                _plotlog(w, Imp.imag, color='r', label='Imag')
            else:
                _plt.plot(w, Imp.real, 'b', label='Real')
                _plt.plot(w, Imp.imag, 'r', label='Imag')
            _plt.legend(loc='best')
            _plt.grid(True)
            _plt.xlabel(r'$\omega [rad/s]$')
            _plt.ylabel(Element._YLABEL[prop])
            _plt.title(self.name+': '+_TITLE[prop])
            if save:
                _plt.savefig(
                    _os.path.sep.join((self.path, prop + name + '.svg')))
        if show:
            _plt.show()


class Budget(list):

    _YLABEL = {
        'Zll': r'$Z_l [k\Omega]$',
        'Zdy': r'$\beta_yZ_y^D [k\Omega]$',
        'Zdx': r'$\beta_xZ_x^D [k\Omega]$',
        'Zqy': r'$\beta_yZ_y^Q [k\Omega]$',
        'Zqx': r'$\beta_xZ_x^Q [k\Omega]$'}

    def __init__(self, lista=None, name=None, path=None):
        lista = lista or []
        if lista and not isinstance(lista[0], Element):
            assert 'Input must be a sequence of Element objects.'
        super().__init__(lista)
        self.name = name or 'Budget'
        if path is not None:
            path = _os.path.abspath(path)
        self.path = path or _os.path.abspath('.')

    def __setitem__(self, k, v):
        assert isinstance(v, Element)
        super().__setitem__(k, v)

    def __setattr__(self, name, value):
        if name in {'name','path'}:
            self.__dict__[name] = str(value)
        else:
            raise AttributeError('Attribute '+name+' is read only.')

    def __getattr__(self, name):
        if name not in _IMPS | _WAKES | {'w', 's'}:
            return [getattr(x, name) for x in self]

        w = _np.unique(_np.concatenate([getattr(x, 'w') for x in self]))
        if name == 'w':
            return w
        if name in _IMPS:
            temp = _np.zeros(w.shape, dtype=complex)
            for el in self:
                attr = getattr(el, name)
                if attr is None or len(attr) == 0:
                    continue
                tmp = _np.interp(w, el.w, attr.imag, left=0.0, right=0.0)*1j
                tmp += _np.interp(w, el.w, attr.real, left=0.0, right=0.0)
                tmp *= el.quantity*_BETA[name](el)
                temp += tmp
            return temp

        s = _np.unique(_np.concatenate([getattr(x, 's') for x in self]))
        if name == 's':
            return s
        if name in _WAKES:
            temp = _np.zeros(s.shape, dtype=float)
            for el in self:
                attr = getattr(el, name)
                if attr is None or len(attr) == 0:
                    continue
                tmp = _np.interp(s, el.s, attr, left=0.0, right=0.0)
                tmp *= el.quantity*_BETA[name](el)
                temp += tmp
            return temp
        raise AttributeError(
            "'" + self.__class__.__name__ + "' object has no attribute '" +
            name + "'")

    def __str__(self):
        string = '{0:^48s}\n'.format(self.name)
        string += '{0:^15s}: {1:^10s} {2:^10s} {3:^10s}\n'.format(
            'Element', 'Quantity', 'Betax', 'Betay')
        for el in self:
            string += '{0:<15s}: {1:^10d} {2:^10.1f} {3:^10.1f}\n'.format(
                el.name, el.quantity, el.betax, el.betay)
        string += '\n'
        return string

    def copy(self):
        other = Budget(name=self.name, path=self.path)
        for el in self:
            other.append(el.copy())
        return other

    def budget2element(self, name=None, path=None):
        ele = Element(name=name, path=path)
        for prop in _IMPS | _WAKES | {'w', 's'}:
            Imp2 = getattr(self, prop)
            if not _np.isclose(Imp2, 0).all():
                setattr(ele, prop, Imp2.copy())
        ele.betax = 1.0
        ele.betay = 1.0
        ele.quantity = 1
        return ele

    def save(self):
        name = self.name.replace(' ', '_').lower()
        _lnls.utils.save_pickle(
            _os.path.sep.join([self.path, name]), budget=self)

    def load(self):
        name = self.name.replace(' ', '_').lower()
        data = _lnls.utils.load_pickle(_os.path.sep.join([self.path, name]))
        return data['budget']

    def plot(
            self, props='all', logscale=True, show=True, save=False,
            name='', figsize=(8, 6), fontsize=14, linewidth=1.5):

        color_map = _plt.get_cmap('nipy_spectral')
        if name:
            name = '_'+name
        props = _prepare_props(props)

        for prop in props:
            a = True
            for el in self:
                Imp3 = getattr(el, prop)
                a &= Imp3 is None or len(Imp3) == 0
                if not a:
                    break
            if a:
                continue
            f, ax = _plt.subplots(2, 1, sharex=True, figsize=figsize)
            N = len(self)
            for i, el in enumerate(self):
                Imp2 = getattr(el, prop)
                if Imp2 is None or len(Imp2) == 0:
                    continue
                Imp = Imp2*_FACTOR[prop] * el.quantity * _BETA[prop](el)
                w = el.w
                cor = color_map(i/N)
                if logscale:
                    _plotlog(
                        w, Imp.real, color=cor, ax=ax[0], linewidth=linewidth)
                    _plotlog(
                        w, Imp.imag, color=cor, label=el.name, ax=ax[1],
                        linewidth=linewidth)
                else:
                    ax[0].plot(w, Imp.real, color=cor, linewidth=linewidth)
                    ax[1].plot(
                        w, Imp.imag, color=cor, label=el.name,
                        linewidth=linewidth)
            ax[1].legend(loc='best', fontsize=10)
            ax[0].grid(True)
            ax[0].tick_params(labelsize=fontsize)
            ax[1].grid(True)
            ax[1].tick_params(labelsize=fontsize)
            ax[1].set_xlabel(r'$\omega [rad/s]$', fontsize=fontsize)
            ax[0].set_ylabel(r'Re'+Budget._YLABEL[prop], fontsize=fontsize)
            ax[1].set_ylabel(r'Im'+Budget._YLABEL[prop], fontsize=fontsize)
            ax[0].set_title(self.name+': '+_TITLE[prop], fontsize=fontsize)
            if save:
                f.savefig(_os.path.sep.join((self.path, prop + name + '.svg')))
        if show:
            _plt.show()


def load_budget(fname):
    data = _lnls.utils.load_pickle(fname)
    return data['budget']


def load_element(fname):
    data = _lnls.utils.load_pickle(fname)
    return data['element']


def longitudinal_resonator(Rs, Q, wr, w):
    """Returns the longitudinal resonator impedance for w.

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
    Zl = w*Rs / (w+1j*Q*(wr - w**2/wr))

    Zl = _np.squeeze(Zl.sum(0))
    Zl = _np.array(Zl, ndmin=1)
    if not isarr:
        Zl = Zl[0]
    return Zl


def transverse_resonator(Rs, Q, wr, w):
    """Returns the transverse resonator impedance for w.

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
    phase = wrl*spos[sel]/_c

    wake = 2*alpha * Rs * _np.exp(-alpha*spos[sel]/_c)
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

    wake = Rs * wr**2 / (Q*wrl) * _np.exp(-alpha*spos[sel]/_c)
    wake *= _np.sin(wrl*spos[sel]/_c)

    Wt = _np.zeros(spos.shape)
    Wt[sel] = wake.sum(0).flatten()
    if not isarr:
        Wt = Wt[0]
    return Wt


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
    gamma = E/E0
    beta = _np.sqrt(1 - 1/gamma**2)
    nu = _np.abs(w/_c)*_np.sqrt(1 - beta**2*epr*mur)

    Zl = 1j*L*w / (2*_np.pi*_ep0 * (beta*_c)**2*gamma**2)
    Zl *= alphaTM(0, epr, mur, beta, nu, b)

    Zv = 1j*L*w**2 / (4*_np.pi*_ep0*_c**2*(beta*_c)*gamma**4)
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
    gamma = E/E0
    beta = _np.sqrt(1 - 1/gamma**2)
    nu = _np.ones((epr.shape[0], 1))*abs(w/_c)*_np.sqrt(1 - beta**2*epr*mur)

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
        Zl_ *= 1j*L*w_ / (2*_np.pi*_ep0 * (beta*_c)**2*gamma**2)
        Zl.append(Zl_)

        Zv_ = complex(alphaTM(1, epr_, mur_, beta_, nu_, b_))
        Zv_ *= 1j*L*w_**2 / (4*_np.pi*_ep0*_c**2*(beta*_c)*gamma**4)
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


def kicker_coupled_flux(w, h, W, t, L, mur, Zg):
    # Calculates Impedances for a ferrite kicker:
    #   - For the Coupled Flux, it uses Davino-Hahn model.
    #
    #   DAVINO-HAHN MODEL:
    #
    #  #######################################    |
    #  ###############FERRITE#################    t
    #  #######################################    |
    #  ###**                             **###  |
    #  ###**  VACUUM                     **###  |      ______|  |_________|
    #  ###**                             **###  |            |  |         |
    #  ###**             +   .           **###  w            |  |         |
    #  ###**                             **###  |            )||(         \
    #  ###**             |_D_|           **###  |      Zk  L1)||(L2     Zg/
    #  ###**                             **###  |            )||(         \
    #  #######################################               | M|         /
    #  #######################################               |  |         |
    #  #######################################         ______|  |_________|
    #      |______________h______________|
    #
    # Bibliografias:
    #
    # - Davino_D Hahn_H - Improved Analytical Model of the transverse coupling
    #   impedance of ferrite kicker magnets - Phys. Rev. ST-AB v6 012001 2003
    #
    # - Nassibian G Sacherer F - Methods for measuring tranverse coupling
    #   impedances in circular Accelerators,
    #   Nucl Inst and Meth. 159 21-27 1979

    # Equivalent Circuit model.
    D = 0.5e-3
    M = L*D*_mu0/W
    L2 = L*h*_mu0/W
    # L2 = L*h*_mu0/W*(mur*t/(mur*t+h*(h/W+1)))

    Zk = (M/L2)**2 * Zg*L2*1j/(1j*w*L2 + Zg) * w
    Zx = (M/L2)**2 * Zg*L2*1j/(1j*w*L2 + Zg) * _c/D**2

    return Zk.conj(), Zx.conj()  # conjugate to adapt impedance convention


def kicker_tsutsui_model(w, epr, mur, a, b, d, L, n):
    '''
       - For the Uncoupled Flux, we can choose between three models:

           TSUTSUI MODEL:

      ******************PEC**********************
      *******************************************
      **#######################################**      |
      **################FERRITE################**      |
      **#######################################**      |
      **                                       **  |   d
      **                                       **  b   |
      **                                       **  |   |
      **     VACUUM        .                   **  |   |
      **                                       **
      **                                       **
      **                                       **
      **#######################################**
      **#######################################**
      **#######################################**
      *******************************************
      *******************************************
                           |__________a________|

     Inputs:

     w   = vector of angular frequencies to evaluate impedances [rad/s]
     epr = vector with real and imaginary electric permeability of ferrite for
           the frequency range of interest
     mur = vector with real and imaginary magnetic permissivity of ferrite for
           the frequency range of interest
     n   = max order of terms to sum
     L   = length of the structure [m]

     Outputs:

     Zl = Impedancia Longitudinal [Ohm]
     Zh = Impedancia Horizontal [Ohm/m]
     Zv = Impedancia Vertical   [Ohm/m]

     Bibliografias:

     - Tsutsui_H - Some Simplified Models of Ferrite Kicker Magnet for
       Calculation of longitudinal Coupling Impedance - CERN-SL-2000-004

     - Tsutsui_H - Transverse Coupling Impedance of a Simplified Ferrite
       Kicker Magnet Model - LHC Project Note 234 - 2000

     - Salvant, B. et al - Quadrupolar Impedance of Simple
     Models of Kickers, Proceedings of IPAC 2010, pp. 2054-2057
     '''

    # Valores do artigo do Wang et al para testar a implementacao das formulas
    # do modelo do Tsutui.
    # a = 103e-3
    # b = 67e-3
    # d = 122e-3
    # L = 1.0
    # Valores do artigo do Tsutsui para testar a implementacao das formulas
    # a = 20e-3
    # b = 16e-3
    # d = 66e-3
    # L = 0.6

    # Terms for the infinite sum:
    n = _np.arange(0, n+1)[:, None]

    k = _np.ones(n.shape)*w/_c
    epr = _np.ones(n.shape)*epr
    mur = _np.ones(n.shape)*mur

    kxn = _np.repeat((2*n + 1)*_np.pi/2/a, w.shape[0], axis=1)
    kyn = _np.sqrt((epr*mur-1)*k**2 - kxn**2)
    sh = _np.sinh(kxn*b)
    ch = _np.cosh(kxn*b)
    tn = _np.tan(kyn*(b-d))
    ct = 1/_np.tan(kyn*(b-d))

    Zl = 1j*_Z0*L/2/a / (
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(
            mur*sh**2*tn - epr*ch**2*ct)
         )/(epr*mur-1) - k/kxn*sh*ch)
    Zq = -Zl*kxn*kxn/k
    Zq = Zq.sum(0)
    Zl = Zl.sum(0)

    Zv = 1j*_Z0*L/2/a * kxn**2/k/(
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(
            mur*ch**2*tn - epr*sh**2*ct)
         )/(epr*mur-1) - k/kxn*sh*ch)
    Zv = Zv.sum(0)

    kxn = _np.repeat(2*(n + 1)*_np.pi/2/a, w.shape[0], axis=1)
    kyn = _np.sqrt((epr*mur-1)*k**2 - kxn**2)
    sh = _np.sinh(kxn*b)
    ch = _np.cosh(kxn*b)
    tn = _np.tan(kyn*(b-d))
    ct = 1/_np.tan(kyn*(b-d))

    Zh = 1j*_Z0*L/2/a * kxn**2/k/(
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(
            mur*sh**2*tn - epr*ch**2*ct)
         )/(epr*mur-1) - k/kxn*sh*ch)
    Zh = Zh.sum(0)

    return Zl.conj(), Zh.conj(), Zv.conj(), Zq.conj()  # impedance convention


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


def taper(w, r1, r2, t, wid=0, geom='round'):
    '''
                   L2
    _____|- - - - - - - - - - - - |_____
          \                      /    :
           \                    /     :
            \                  /      :
             \_______L1_______/       : R2
                :                     :
                : R1                  :
                :                     :
    - - - - - - - - - - - - - - - - - - -
    '''
    diff = _np.abs(r2-r1)
    summ = r2 + r1
    prod = r2*r1
    ums = _np.ones(w.shape)
    if geom == 'round':
        Zll = 2 * -1j*w*_Z0/4/_np.pi/_c * diff/t
        Zdx = 2 * -1j*_Z0/2/_np.pi * diff/t/prod * ums
        Zdy = 1*Zdx
        Zqx = 0*w
    else:
        Zll = 2 * -1j*0.43*w*_Z0/_np.pi/_c * diff/t
        Zdx = 2 * -1j*_Z0/4/_np.pi * diff/t/prod * ums
        Zdy = 2 * -1j*_Z0/2 * wid/t*summ*diff/prod/prod * ums
        Zqx = -1*Zdx
    return Zll, Zdx, Zdy, Zqx


class CSRElement:
    _X = _np.linspace(-900, 900, 10001)
    _Y = None

    def __init__(
            self, rho=17.2, h=12e-3, bl=2.5e-3, nus=4.6e-3,
            espread=8.5e-4, rev_time=1.73e-6):
        self.bl = bl
        self.rho = rho
        self.h = h
        self.nus = nus
        self.espread = espread
        self.rev_time = rev_time

    @property
    def shielding(self):
        return self.rho**(1/2) * self.bl / self.h**(3/2)

    @property
    def threshold(self):
        S = self.calc_normalized_strength(1e-3) / 1e-3
        Ith = (0.5 + 0.12*self.shielding)/S
        return Ith

    def calc_normalized_current(self, I0):
        return (120/4*_c*I0*self.rev_time /
                (2*_np.pi*self.nus*3e9*self.espread))

    def calc_normalized_strength(self, I0):
        curr = self.calc_normalized_current(I0)
        return curr*self.rho**(1/3)/self.bl**(4/3)

    def calc_formation_length(self, w):
        return (24*self.rho**2/w/_c)**(1/3)

    def wake(self, z, maxi=25, L=None, bunlen=80e-6, convolved=True):
        gaus_to_si = 120*_c/4  # _Z0 _c / (4 pi)
        L = L or (2 * _np.pi * self.rho)

        # Free Space Term
        inds = z < 0
        W0 = _np.zeros(len(z))
        W0[inds] = (-2/3**(4/3)/self.rho**(2/3) /
                    _np.power((-z[inds]), 4/3)*gaus_to_si*L)

        # Shielding Term
        W1 = _np.zeros(len(z))
        zshield = self.shielding / self.bl * z
        for i in range(1, maxi):
            uai = self._getY(3*zshield/i**(3/2))
            W1 += 8*_np.pi*(-1)**(i+1)/i/i*uai*(3 - uai)/(1 + uai)**3
        W1 *= -1/self.h**2 / (2 * _np.pi) * gaus_to_si * L

        # If want the convolved values
        if convolved:
            # For the free space used analytical formulas of convolution
            bl = bunlen
            C = _Z0*_c/2**(13/6)/_np.pi**(3/2)/(3*self.rho**2*bl**10)**(1/3)*L
            W0 = C*(2**(1/2)*_gamma(5/6)*(
                        bl**2*_hyp1f1(-1/3, 1/2, -z*z/2/bl**2) -
                        z**2*_hyp1f1(2/3, 3/2, -z*z/2/bl**2)) +
                    z*bl*_gamma(4/3)*(
                        3*_hyp1f1(1/6, 1/2, -z*z/2/bl**2) -
                        2*_hyp1f1(1/6, 3/2, -z*z/2/bl**2)))
            # For shielding perform convolution numerically
            bunch = _np.exp(-(z*z/bl**2)/2)/_np.sqrt(2*_np.pi)/bl  # gaussian
            W1 = _scysig.fftconvolve(W1, bunch, mode='same') * (z[1]-z[0])
        return W0, W1

    def impedance(self, w, L=None, imax=3, free=False):
        L = L or (2*_np.pi * self.rho)
        Lfrac = L / (2*_np.pi * self.rho)
        n = w/_c * self.rho

        if free:
            Z = (120/2 * 1.354/3**(1/3)*_np.exp(1j*_np.pi/6) *
                 (w/_c/self.rho**2)**(1/3) * L)
            return Z

        u0 = _np.pi**2/2**(2/3) * (n * (2*self.h/self.rho)**(3/2))**(-4/3)
        Z = (1 + 0*1j)*_Z0 * 16 * u0

        F = _np.zeros(len(w), dtype=complex)
        for p in range(0, imax):
            up = u0*(2*p + 1)**2
            Ai, Ail, Bi, Bil = _airy(up)
            Ri = Ail*Ail + up * Ai*Ai
            Ai, Ail, Bi, Bil = _airye(up)
            Im = Ail*Bil + up * Ai*Bi
            F += Ri - 1j * Im
        Z *= F
        return Z * n * self.h / self.rho * Lfrac

    def _getY(self, values):
        if self._Y is None:
            self._Y = _np.zeros(self._X.shape)
            for i, x in enumerate(self._X):
                sol = _np.roots([1, 0, 0, -x, -3])
                sol2 = []
                for s in sol:
                    if _np.isclose(s.imag, 0):
                        sol2.append(s.real)
                if len(sol2) != 2:
                    print('war')
                self._Y[i] = min(sol2)
            self._Y = self._Y*self._Y*self._Y*self._Y
        return _np.interp(values, self._X, self._Y, left=None, right=None)


def from_wake_to_impedance(z, wake, bunlen, cutoff=2):
    dt = (z[-1]-z[0]) / (z.shape[0]-1) / _c
    VHat = _np.fft.fft(wake) * dt
    w = 2 * _np.pi * _np.fft.fftfreq(len(z), d=dt)
    VHat = _np.fft.fftshift(VHat)   # shift the negative frequencies
    w = _np.fft.fftshift(w)         # to the center of the spectrum
    VHat *= _np.exp(-1j * w * z[0] / _c)

    wmax = cutoff / bunlen * _c
    indcs = _np.abs(w) <= wmax
    # Deconvolve the Transform with a gaussian bunch:
    Jwlist = _np.exp(-(w * bunlen / _c)**2 / 2)
    Z = VHat[indcs] / Jwlist[indcs]
    w = w[indcs]
    return w, Z.conj()
