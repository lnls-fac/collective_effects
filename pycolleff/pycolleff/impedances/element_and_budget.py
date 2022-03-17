import os as _os

import numpy as _np
import matplotlib.pyplot as _plt

import lnls as _lnls


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
        if name in {'name', 'path'}:
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
