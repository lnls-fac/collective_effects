"""."""
import os as _os

import numpy as _np
import matplotlib.pyplot as _plt

from mathphys.functions import save_pickle as _save_pickle, \
    load_pickle as _load_pickle


_IMPS = {'Zll', 'Zdy', 'Zdx', 'Zqy', 'Zqx'}
_WAKES = {'Wll', 'Wdy', 'Wdx', 'Wqy', 'Wqx'}
_TITLE = {
    'Zll': 'Longitudinal Impedance',
    'Zdy': 'Driving Vertical Impedance',
    'Zdx': 'Driving Horizontal Impedance',
    'Zqy': 'Detuning Vertical Impedance',
    'Zqx': 'Detuning Horizontal Impedance',
    'Wll': 'Longitudinal Wake-function',
    'Wdy': 'Driving Vertical Wake-function',
    'Wdx': 'Driving Horizontal Wake-function',
    'Wqy': 'Detuning Vertical Wake-function',
    'Wqx': 'Detuning Horizontal Wake-function',
    }
_FACTOR = {
    'Zll': 1e-3, 'Zdy': 1e-3, 'Zdx': 1e-3, 'Zqy': 1e-3, 'Zqx': 1e-3,
    'Wll': 1e-12, 'Wdy': 1e-12, 'Wdx': 1e-12, 'Wqy': 1e-12, 'Wqx': 1e-12}
_BETA = {
    'Zll': lambda x: 1,       'Wll': lambda x: 1,
    'Zdy': lambda x: x.betay, 'Wdy': lambda x: x.betay,
    'Zdx': lambda x: x.betax, 'Wdx': lambda x: x.betax,
    'Zqy': lambda x: x.betay, 'Wqy': lambda x: x.betay,
    'Zqx': lambda x: x.betax, 'Wqx': lambda x: x.betax,
    }


def _plotlogy(x, y, color=None, label=None, ax=None, linewidth=1.5):
    if ax is None:
        ax = _plt.gca()

    # NOTE: I tried to use `symlog` of matplotlib, but the result is too ugly.
    ax.plot(x, y, color=color, label=label, linewidth=linewidth)
    if (y < 0).any():
        if isinstance(label, str):
            label = '(-1)*'+label
        ax.plot(x, -y, '--', color=color, label=label, linewidth=linewidth)
    ax.set_yscale('log')


def _prepare_props(props):
    allp = _IMPS | _WAKES
    if isinstance(props, str):
        if props.lower() == 'all':
            props = sorted(list(allp))
        elif props in allp:
            props = [props]
        else:
            raise AttributeError(props, ' not supported for plot.')
    elif isinstance(props, (list, tuple)):
        wrong = set(props) - allp
        if wrong:
            raise AttributeError(wrong, ' not supported for plot.')
    else:
        raise TypeError("Type '"+type(props)+"' not supported for 'props'.")
    return props


class Element:
    """."""

    _YLABEL = {
        'Zll': r'$Z_l [k\Omega]$',
        'Zdy': r'$Z_y^D [k\Omega/m]$',
        'Zdx': r'$Z_x^D [k\Omega/m]$',
        'Zqy': r'$Z_y^Q [k\Omega/m]$',
        'Zqx': r'$Z_x^Q [k\Omega/m]$',
        'Wll': r'$W_l [V/pC]$',
        'Wdy': r'$W_y^D [V/pC/m]$',
        'Wdx': r'$W_x^D [V/pC/m]$',
        'Wqy': r'$W_y^Q [V/pC/m]$',
        'Wqx': r'$W_x^Q [V/pC/m]$',
        }

    def __init__(
            self, name='element', betax=0.0, betay=0.0, quantity=1):
        """."""
        self.name = name

        self.quantity = quantity  # this field must only be used in Budget
        self.betax = betax  # this field shall only be used in Budget
        self.betay = betay  # this field shall only be used in Budget

        self.ang_freq = _np.array([], dtype=float)
        self.Zll = _np.array([], dtype=complex)
        self.Zdy = _np.array([], dtype=complex)
        self.Zdx = _np.array([], dtype=complex)
        self.Zqy = _np.array([], dtype=complex)
        self.Zqx = _np.array([], dtype=complex)

        self.pos = _np.array([], dtype=float)
        self.Wll = _np.array([], dtype=float)
        self.Wdy = _np.array([], dtype=float)
        self.Wdx = _np.array([], dtype=float)
        self.Wqy = _np.array([], dtype=float)
        self.Wqx = _np.array([], dtype=float)

    def copy(self):
        """."""
        other = Element(
            name=self.name, betax=self.betax, betay=self.betay,
            quantity=self.quantity)
        other.ang_freq = self.ang_freq.copy()
        other.Zll = self.Zll.copy()
        other.Zdy = self.Zdy.copy()
        other.Zdx = self.Zdx.copy()
        other.Zqy = self.Zqy.copy()
        other.Zqx = self.Zqx.copy()
        other.pos = self.pos.copy()
        other.Wll = self.Wll.copy()
        other.Wdy = self.Wdy.copy()
        other.Wdx = self.Wdx.copy()
        other.Wqy = self.Wqy.copy()
        other.Wqx = self.Wqx.copy()
        return other

    def to_dict(self):
        """."""
        return {
            'name': self.name,
            'betax': self.betax,
            'betay': self.betay,
            'quantity': self.quantity,
            'ang_freq': self.ang_freq.copy(),
            'Zll': self.Zll.copy(),
            'Zdy': self.Zdy.copy(),
            'Zdx': self.Zdx.copy(),
            'Zqy': self.Zqy.copy(),
            'Zqx': self.Zqx.copy(),
            'pos': self.pos.copy(),
            'Wll': self.Wll.copy(),
            'Wdy': self.Wdy.copy(),
            'Wdx': self.Wdx.copy(),
            'Wqy': self.Wqy.copy(),
            'Wqx': self.Wqx.copy(),
            }

    def from_dict(self, dic):
        """."""
        self.name = dic.get('name', self.name)
        self.betax = dic.get('betax', self.betax)
        self.betay = dic.get('betay', self.betay)
        self.quantity = dic.get('quantity', self.quantity)
        self.ang_freq = dic.get('ang_freq', self.ang_freq).copy()
        self.Zll = dic.get('Zll', self.Zll).copy()
        self.Zdy = dic.get('Zdy', self.Zdy).copy()
        self.Zdx = dic.get('Zdx', self.Zdx).copy()
        self.Zqy = dic.get('Zqy', self.Zqy).copy()
        self.Zqx = dic.get('Zqx', self.Zqx).copy()
        self.pos = dic.get('pos', self.pos).copy()
        self.Wll = dic.get('Wll', self.Wll).copy()
        self.Wdy = dic.get('Wdy', self.Wdy).copy()
        self.Wdx = dic.get('Wdx', self.Wdx).copy()
        self.Wqy = dic.get('Wqy', self.Wqy).copy()
        self.Wqx = dic.get('Wqx', self.Wqx).copy()

    def save(self, overwrite=False):
        """."""
        name = self.name.replace(' ', '_').lower()
        dic = self.to_dict()
        _save_pickle(dic, name, overwrite=overwrite)

    def load(self):
        """."""
        name = self.name.replace(' ', '_').lower()
        data = _load_pickle(name)
        self.from_dict(data)

    def plot(
            self, props='all', logx=True, logy=True, show=True, save=False,
            figname='', figsize=(8, 4)):
        """Create plots for properties.

        Args:
            props (str or list, optional): List of properties to plot. Value
                must assume {'Zll', 'Zdx', 'Zdy', 'Zqx', 'Zqy', 'Wll', 'Wdx',
                'Wdy', 'Wqx', 'Wqy'}. Can also assume the string 'all', to
                plot all properties. Defaults to 'all'.
            logx (bool, optional): Whether or not the horizontal scale will be
                in log-scale. Negative values of x will not be displayed.
                Defaults to True.
            logy (bool, optional): Whether or not the vertical scale will be
                in log-scale. Negative values of y will be displayed with '--'
                linestyle. Defaults to True.
            show (bool, optional): Whether or not to show the figures.
                Defaults to True.
            save (bool, optional): If True the figures will be saved to files.
                Defaults to False.
            figname (str, optional): suffix for file name. Defaults to ''.
            figsize (tuple, optional): figure size. Defaults to (8, 6).
            fontsize (int, optional): fontsize of the plots. Defaults to 14.
            linewidth (int, optional): linewidth of the lines. Defaults to 2.

        """
        if figname:
            figname = '_' + figname
        props = _prepare_props(props)

        func = _plotlogy if logy else _plt.plot

        for prop in props:
            is_imp = prop.startswith('Z')
            imp2 = getattr(self, prop)
            if imp2 is None or len(imp2) == 0:
                continue
            fig, ax = _plt.subplots(1, 1, figsize=figsize)
            imp = imp2*_FACTOR[prop]
            if is_imp:
                w = self.ang_freq
                func(w, imp.real, color='b', label='Real')
                func(w, imp.imag, color='r', label='Imag')
                ax.set_xlabel(r'$\omega [rad/s]$')
            else:
                pos = self.pos
                func(pos, imp, color='b')
                ax.set_xlabel('Distance to Source [m]')

            if logx:
                ax.set_xscale('log')
            ax.legend(loc='best')
            ax.grid(True)
            ax.set_ylabel(Element._YLABEL[prop])
            ax.set_title(self.name+': '+_TITLE[prop])
            fig.tight_layout()
            if save:
                fig.savefig(prop + figname + '.svg')
            if show:
                fig.show()


class Budget(list):
    """Collection of impedance source elements of the ring."""

    _YLABEL = {
        'Zll': r'$Z_l [k\Omega]$',
        'Zdy': r'$\beta_y \times Z_y^D [k\Omega]$',
        'Zdx': r'$\beta_x \times Z_x^D [k\Omega]$',
        'Zqy': r'$\beta_y \times Z_y^Q [k\Omega]$',
        'Zqx': r'$\beta_x \times Z_x^Q [k\Omega]$',
        'Wll': r'$W_l [V/pC]$',
        'Wdy': r'$\beta_y \times W_y^D [V/pC]$',
        'Wdx': r'$\beta_x \times W_x^D [V/pC]$',
        'Wqy': r'$\beta_y \times W_y^Q [V/pC]$',
        'Wqx': r'$\beta_x \times W_x^Q [V/pC]$',
        }

    def __init__(self, lista=None, name='budget'):
        """Create a budget element.

        Args:
            lista (list or Budget, optional): list of elements.
                Defaults to None.
            name (str, optional): name of the budget. Defaults to 'budget'.

        Raises:
            ValueError: When items in lista are not of type Element.

        """
        lista = lista or []
        if lista and not isinstance(lista[0], Element):
            raise ValueError('Input must be a sequence of Element objects.')
        super().__init__(lista)
        self._name = name

    def __str__(self):
        """."""
        string = '{0:^48s}\n'.format(self.name)
        string += '{0:^15s}: {1:^10s} {2:^10s} {3:^10s}\n'.format(
            'Element', 'Quantity', 'Betax', 'Betay')
        for el in self:
            string += '{0:<15s}: {1:^10d} {2:^10.1f} {3:^10.1f}\n'.format(
                el.name, el.quantity, el.betax, el.betay)
        string += '\n'
        return string

    def __setitem__(self, k, v):
        """."""
        if not isinstance(v, Element):
            raise ValueError('Item must be instance of Element.')
        super().__setitem__(k, v)

    @property
    def name(self):
        """."""
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def ang_freq(self):
        """."""
        return _np.unique(_np.r_[[getattr(x, 'ang_freq') for x in self]])

    @property
    def Zll(self):
        """."""
        return self._get_impedance('Zll')

    @property
    def Zdx(self):
        """."""
        return self._get_impedance('Zdx')

    @property
    def Zdy(self):
        """."""
        return self._get_impedance('Zdy')

    @property
    def Zqx(self):
        """."""
        return self._get_impedance('Zqx')

    @property
    def Zqy(self):
        """."""
        return self._get_impedance('Zqy')

    @property
    def pos(self):
        """."""
        return _np.unique(_np.r_[[getattr(x, 'pos') for x in self]])

    @property
    def Wll(self):
        """."""
        return self._get_wake('Wll')

    @property
    def Wdx(self):
        """."""
        return self._get_wake('Wdx')

    @property
    def Wdy(self):
        """."""
        return self._get_wake('Wdy')

    @property
    def Wqx(self):
        """."""
        return self._get_wake('Wqx')

    @property
    def Wqy(self):
        """."""
        return self._get_wake('Wqy')

    def copy(self):
        """."""
        other = Budget(name=self.name)
        for el in self:
            other.append(el.copy())
        return other

    def budget2element(self, name=''):
        """Transform a Budget object into an Element object.

        Args:
            name (str, optional): name of the element.
                Defaults to 'element_'+self.name.

        Returns:
            Element: element object.

        """
        name = name or 'element_'+self.name
        ele = Element(name=name)
        for prop in _IMPS | _WAKES | {'ang_freq', 'pos'}:
            Imp2 = getattr(self, prop)
            if not _np.isclose(Imp2, 0).all():
                setattr(ele, prop, Imp2.copy())
        ele.betax = 1.0
        ele.betay = 1.0
        ele.quantity = 1
        return ele

    def to_dict(self):
        """."""
        dic = {'name': self.name}
        dic['elements'] = [el.to_dict() for el in self]
        return dic

    def from_dict(self, dic):
        """."""
        self.name = dic.get('name', self.name)
        if 'elements' not in dic:
            return
        for el_data in dic['elements']:
            el = Element()
            el.from_dict(el_data)
            self.append(el)

    def save(self, overwrite=False):
        """."""
        name = self.name.replace(' ', '_').lower()
        dic = self.to_dict()
        _save_pickle(dic, name, overwrite=overwrite)

    def load(self):
        """."""
        name = self.name.replace(' ', '_').lower()
        data = _load_pickle(name)
        self.from_dict(data)

    def plot_impedances(
            self, props='all', logx=True, logy=True, show=True, save=False,
            figname='', figsize=(8, 6), fontsize=14, linewidth=2):
        """Create plots for impedances.

        Args:
            props (str or list, optional): List of properties to plot. Value
                must assume {'Zll', 'Zdx', 'Zdy', 'Zqx', 'Zqy'}. Can also
                assume the string 'all', to plot all properties.
                Defaults to 'all'.
            logx (bool, optional): Whether or not the horizontal scale will be
                in log-scale. Negative values of x will not be displayed.
                Defaults to True.
            logy (bool, optional): Whether or not the vertical scale will be
                in log-scale. Negative values of y will be displayed with '--'
                linestyle. Defaults to True.
            show (bool, optional): Whether or not to show the figures.
                Defaults to True.
            save (bool, optional): If True the figures will be saved to files.
                Defaults to False.
            figname (str, optional): suffix for file name. Defaults to ''.
            figsize (tuple, optional): figure size. Defaults to (8, 6).
            fontsize (int, optional): fontsize of the plots. Defaults to 14.
            linewidth (int, optional): linewidth of the lines. Defaults to 2.

        """
        color_map = _plt.get_cmap('nipy_spectral')
        if figname:
            figname = '_' + figname
        props = _prepare_props(props)

        fun = _plotlogy if logy else _plt.plot

        for prop in props:
            a = True
            for el in self:
                Imp3 = getattr(el, prop)
                a &= Imp3 is None or len(Imp3) == 0
                if not a:
                    break
            if a:
                continue
            fig, axs = _plt.subplots(2, 1, sharex=True, figsize=figsize)
            size = len(self)
            for i, el in enumerate(self):
                imp2 = getattr(el, prop)
                if imp2 is None or len(imp2) == 0:
                    continue
                imp = imp2*_FACTOR[prop] * el.quantity * _BETA[prop](el)
                w = el.ang_freq
                cor = color_map(i/size)
                _plt.sca(axs[0])
                fun(w, imp.real, color=cor, linewidth=linewidth)
                _plt.sca(axs[1])
                fun(w, imp.imag, color=cor, label=el.name, linewidth=linewidth)

            if logx:
                axs[0].set_xscale('log')
            axs[1].legend(loc='best', fontsize=10)
            axs[0].grid(True)
            axs[0].tick_params(labelsize=fontsize)
            axs[1].grid(True)
            axs[1].tick_params(labelsize=fontsize)
            axs[1].set_xlabel(r'$\omega [rad/s]$', fontsize=fontsize)
            axs[0].set_ylabel(r'Re'+Budget._YLABEL[prop], fontsize=fontsize)
            axs[1].set_ylabel(r'Im'+Budget._YLABEL[prop], fontsize=fontsize)
            axs[0].set_title(self.name+': '+_TITLE[prop], fontsize=fontsize)
            fig.tight_layout()
            if save:
                fig.savefig(prop + figname + '.svg')
            if show:
                fig.show()

    def plot_wakes(
            self, props='all', logx=True, logy=True, show=True, save=False,
            figname='', figsize=(8, 6), fontsize=14, linewidth=2):
        """Create plots for wakes.

        Args:
            props (str or list, optional): List of properties to plot. Value
                must assume {'Wll', 'Wdx', 'Wdy', 'Wqx', 'Wqy'}. Can also
                assume the string 'all', to plot all properties.
                Defaults to 'all'.
            logx (bool, optional): Whether or not the horizontal scale will be
                in log-scale. Negative values of x will not be displayed.
                Defaults to True.
            logy (bool, optional): Whether or not the vertical scale will be
                in log-scale. Negative values of y will be displayed with '--'
                linestyle. Defaults to True.
            show (bool, optional): Whether or not to show the figures.
                Defaults to True.
            save (bool, optional): If True the figures will be saved to files.
                Defaults to False.
            figname (str, optional): suffix for file name. Defaults to ''.
            figsize (tuple, optional): figure size. Defaults to (8, 6).
            fontsize (int, optional): fontsize of the plots. Defaults to 14.
            linewidth (int, optional): linewidth of the lines. Defaults to 2.

        """
        color_map = _plt.get_cmap('nipy_spectral')
        if figname:
            figname = '_' + figname
        props = _prepare_props(props)

        fun = _plotlogy if logy else _plt.plot

        for prop in props:
            a = True
            for el in self:
                imp3 = getattr(el, prop)
                a &= imp3 is None or len(imp3) == 0
                if not a:
                    break
            if a:
                continue
            fig, ax = _plt.subplots(1, 1, figsize=figsize)
            size = len(self)
            for i, el in enumerate(self):
                wake = getattr(el, prop)
                if wake is None or len(wake) == 0:
                    continue
                wake = wake*_FACTOR[prop] * el.quantity * _BETA[prop](el)
                pos = el.pos
                cor = color_map(i/size)
                fun(pos, wake, color=cor, linewidth=linewidth, label=el.name)

            if logx:
                ax.set_xscale('log')
            ax.legend(loc='best', fontsize=10)
            ax.grid(True)
            ax.tick_params(labelsize=fontsize)
            ax.set_xlabel('Distance to Source [m]', fontsize=fontsize)
            ax.set_ylabel(Budget._YLABEL[prop], fontsize=fontsize)
            ax.set_title(self.name+': '+_TITLE[prop], fontsize=fontsize)
            fig.tight_layout()
            if save:
                fig.savefig(prop + figname + '.svg')
            if show:
                fig.show()

    def _get_impedance(self, name):
        ang_freq = self.ang_freq

        temp = _np.zeros(ang_freq.shape, dtype=complex)
        kws = {'left': 0.0, 'right': 0.0}
        for el in self:
            attr = getattr(el, name)
            if attr is None or len(attr) == 0:
                continue
            tmp = _np.interp(ang_freq, el.ang_freq, attr.imag, **kws)*1j
            tmp += _np.interp(ang_freq, el.ang_freq, attr.real, **kws)
            tmp *= el.quantity*_BETA[name](el)
            temp += tmp
        return temp

    def _get_wake(self, name):
        pos = self.pos
        temp = _np.zeros(pos.shape, dtype=float)
        for el in self:
            attr = getattr(el, name)
            if attr is None or len(attr) == 0:
                continue
            tmp = _np.interp(pos, el.pos, attr, left=0.0, right=0.0)
            tmp *= el.quantity*_BETA[name](el)
            temp += tmp
        return temp


def load_budget(fname):
    """Load impedance budget from file.

    Args:
        fname (str): path and name of the file to be loaded.

    Returns:
        Budget: budget object.

    """
    data = _load_pickle(fname)
    bud = Budget()
    bud.from_dict(data)
    return bud


def load_element(fname):
    """Load impedance element from file.

    Args:
        fname (str): path and name of the file to be loaded.

    Returns:
        Element: element object.

    """
    data = _load_pickle(fname)
    ele = Element()
    ele.from_dict(data)
    return ele
