"""Analysis of Wake simulation codes results.

Troughout the code I am assuming:
    s positive means particle behind source -->  Wl, Wt = 0 s < 0

Definitions:
    - Longitudinal Wake:
        Wl(s) = -c/Q * int El(ct-s,t) dt

    - Dipolar wake:
        Wdx(s) = - 1/xs * int_-inf^s dWl(s)/dx ds'
    where xs is source particle displacement

    - Quadrupolar wake:
        Wqx(s) = - int_-inf^s dWl(s)/dx^2 ds'

    - Longitudinal impedance:
        Zl = int exp(i*w*s) Wl(s) ds

    - Dipolar and Quadrupolar impedances:
        Zx = i*int exp(i*w*s) Wx(s) ds

"""

import stat as _stat
from os import mkdir as _mkdir, chmod as _chmod, lstat as _lstat
from os.path import join as _jnpth, abspath as _abspth, isdir as _isdir
import gzip as _gzip
import pickle as _pickle
import logging as _log

import numpy as _np
import matplotlib.pyplot as _plt
from matplotlib import rc as _rc

from mathphys.constants import light_speed as c

from pycolleff import sirius as _sirius

try:
    from pyaccel import naff as _naff
    bool_pyaccel = True
except Exception:
    bool_pyaccel = False


ANALYSIS_TYPES = {
    'dx',  # horizontal impedance
    'dy',  # vertical impedance
    'db',  # both planes are symmetric
    'll'   # longitudinal and transverse quadrupolar impedances
    }


class EMSimulData:
    """Class that gather processed simulated results."""

    DEFAULT_FNAME_SAVE = 'SimulData.pickle'
    PLANES = ('ll', 'dx', 'dy', 'qx', 'qy')
    TITLES = {
        'll': 'Longitudinal',
        'dx': 'Dipolar Horizontal',
        'dy': 'Dipolar Vertical',
        'qx': 'Quadrupolar Horizontal',
        'qy': 'Quadrupolar Vertical'}
    WAKE_YLABELS = {
        'll': r'$W_l$ [V/pC]',
        'dx': r'$W_{{D_x}}$ [V/pC/m]',
        'dy': r'$W_{{D_y}}$ [V/pC/m]',
        'qx': r'$W_{{Q_x}}$ [V/pC/m]',
        'qy': r'$W_{{Q_y}}$ [V/pC/m]'}
    IMPS_YLABELS = {
        'll': r'$Z_l$ [$\Omega$]',
        'dx': r'$Z_{{D_x}}$ [$\Omega$/m]',
        'dy': r'$Z_{{D_y}}$ [$\Omega$/m]',
        'qx': r'$Z_{{Q_x}}$ [$\Omega$/m]',
        'qy': r'$Z_{{Q_y}}$ [$\Omega$/m]'}

    def __init__(self, code=None):
        """."""
        # Sirius model to perform calculations
        self._si = _sirius.create_ring()
        # CST, ACE3P, GdfidL, ECHOz1 ECHOz2, ECHOzR, ECHO2D, ECHO3D...
        self.code = code
        # Bunch Length Used in simulation [m]
        self.bunlen = 0.0
        # positions where bunch is defined [m]
        self.sbun = _np.array([], dtype=float)
        # bunch profile used in the simulation [As/m]
        self.bun = _np.array([], dtype=float)
        # axis: distance from driving bunch  [m]
        # (positive when particle is behind driving bunch)
        self.s = _np.array([], dtype=float)
        # Longitudinal Wakepotential [V/C]
        self.Wll = _np.array([], dtype=float)
        # Dipolar Horizontal Wakepotential [V/C/m]
        self.Wdx = _np.array([], dtype=float)
        # Dipolar Vertical Wakepotential [V/C/m]
        self.Wdy = _np.array([], dtype=float)
        # Quadrupolar Horizontal Wakepotential [V/C/m]
        self.Wqx = _np.array([], dtype=float)
        # Quadrupolar Vertical Wakepotential [V/C/m]
        self.Wqy = _np.array([], dtype=float)
        # axis: frequency obtained from FFT [GHz]
        self.freq = _np.array([], dtype=float)
        # Longitudinal Impedance [Ohm]
        self.Zll = _np.array([], dtype=complex)
        # Dipolar Horizontal Impedance [Ohm/m]
        self.Zdx = _np.array([], dtype=complex)
        # Dipolar Vertical Impedance [Ohm/m]
        self.Zdy = _np.array([], dtype=complex)
        # Quadrupolar Horizontal Impedance [Ohm/m]
        self.Zqx = _np.array([], dtype=complex)
        # Quadrupolar Vertical Impedance [Ohm/m]
        self.Zqy = _np.array([], dtype=complex)
        self._klossW = None
        self._kckdxW = None
        self._kckdyW = None
        self._kckqxW = None
        self._kckqyW = None

    def copy(self):
        """."""
        other = EMSimulData()
        other.code = self.code
        other.bunlen = self.bunlen
        other.sbun = self.sbun.copy()
        other.bun = self.bun.copy()
        other.s = self.s.copy()
        other.Wll = self.Wll.copy()
        other.Wdx = self.Wdx.copy()
        other.Wdy = self.Wdy.copy()
        other.Wqx = self.Wqx.copy()
        other.Wqy = self.Wqy.copy()
        other.freq = self.freq.copy()
        other.Zll = self.Zll.copy()
        other.Zdx = self.Zdx.copy()
        other.Zdy = self.Zdy.copy()
        other.Zqx = self.Zqx.copy()
        other.Zqy = self.Zqy.copy()
        other._klossW = self._klossW
        other._kckdxW = self._kckdxW
        other._kckdyW = self._kckdyW
        other._kckqxW = self._kckqxW
        other._kckqyW = self._kckqyW
        return other

    def get_klossW(self):
        """."""
        if self._klossW:
            return self._klossW
        sigs, spos = self.bunlen, self.s
        wake = self.Wll
        if wake is None or _np.allclose(wake, 0, atol=0):
            return None
        rhos = (1/(sigs*_np.sqrt(2*_np.pi)))*_np.exp(-spos**2/(2*sigs**2))
        kW = _np.trapz(wake*rhos, x=spos)
        self._klossW = kW
        return kW

    def get_PlossW(self, rev_time=None, harm_num=None, Iavg=500e-3):
        """."""
        rev_time = rev_time or self._si.T0
        harm_num = harm_num or self._si.harm_num
        kW = self.klossW()
        Ploss = kW * Iavg**2 * rev_time * 1e12 / harm_num
        return Ploss

    def get_kick_factorW(self, pl='dy'):
        """."""
        kick = getattr(self, '_kck'+pl+'W')
        if kick:
            return kick
        sigs, spos = self.bunlen, self.s
        wake = getattr(self, 'W'+pl)
        if wake is None or _np.all(wake == 0):
            return None
        rhos = (1/(sigs*_np.sqrt(2*_np.pi)))*_np.exp(-spos**2/(2*sigs**2))
        kW = _np.trapz(wake*rhos, x=spos)
        setattr(self, '_kck'+pl+'W', kW)
        return kW

    def get_klossZ(self, bunlen=2.65e-3, nbunches=1):
        """."""
        self._si.nbun = nbunches
        klossZ, *_ = self._si.loss_factor(
            w=self.freq*2*_np.pi, Zl=self.Zll, bunlen=bunlen)
        return klossZ

    def get_kick_factorZ(self, pl='dy', bunlen=2.65e-3, nbunches=1):
        """."""
        self._si.nbun = nbunches
        Z = getattr(self, 'Z'+pl)
        if Z is None or _np.allclose(Z, 0, atol=0):
            return None
        kckZ, *_ = self._si.kick_factor(
            w=self.freq*2*_np.pi, Z=Z, bunlen=bunlen)
        return kckZ

    def get_PlossZ(self, bunlen=2.65e-3, nbunches=1):
        """."""
        self._si.nbun = nbunches
        _, PlossZ, *_ = self._si.loss_factor(
            w=self.freq*2*_np.pi, Zl=self.Zll, bunlen=bunlen)
        return PlossZ

    def calc_impedance(
            self, use_win='phase', pl=None, cutoff=2, s_min=None,
            s_max=None, silent=False):
        """."""
        _log.basicConfig(level=_log.CRITICAL if silent else _log.INFO)
        _log.info('#'*60 + '\n' + 'Calculating Impedances')

        planes = self.PLANES
        if pl is not None:
            planes = [pl]

        # Extracts Needed Variables
        sigt = self.bunlen / c  # bunch time-length
        spos = self.s.copy()

        if s_min is None:
            s_min = spos[0]
        if s_max is None:
            s_max = spos[-1]
        inds = (spos >= s_min) & (spos <= s_max)
        spos = spos[inds]
        # numpy's fft algorithm is slow for large primes:
        p, n = self._choose_fft_length(spos.shape[0])
        spos = spos[:n]

        if p > 0:
            _log.info(
                'Last {0:d} point{1:s} of wake {2:s} '.format(
                    p, *(('s', 'were') if p > 1 else ('', 'was'))) +
                'not considered to gain performance in FFT.')

        if use_win is True:
            _log.info('Using Half-Hanning Window')
            # Half Hanning window to zero the end of the signal
            window = _np.hanning(2*spos.shape[0])[spos.shape[0]:]
        elif isinstance(use_win, str) and use_win.lower().startswith('phase'):
            _log.info('Using Half-Hanning Window to correct the phases')
            # Half Hanning window to smooth the final of the signal
            window = _np.hanning(2*spos.shape[0])[spos.shape[0]:]
        else:
            _log.info('Not using Window')
            window = _np.ones(spos.shape[0])

        _log.info(f'Cutoff frequency w = {cutoff:d}/sigmat')

        for pl in planes:
            _log.info(f'Performing FFT on W{pl:s}: ')
            Wpl = getattr(self, 'W'+pl).copy()
            if Wpl is None or _np.all(Wpl == 0):
                _log.info('No Data found.')
                continue
            _log.info('Data found. ')
            Wpl = Wpl[inds]
            Wpl = Wpl[:n]

            self.freq, Zpl = self._get_impedance(
                spos, Wpl*window, sigt, cutoff)
            if isinstance(use_win, str) and use_win.lower().startswith('phase'):
                _, Zpl2 = self._get_impedance(spos, Wpl, sigt, cutoff)
                Zpl = _np.abs(Zpl2)*_np.exp(1j*_np.angle(Zpl))

            # I have to take the conjugate of the fft because:
            # fftt == \int exp(-i*2pi*f*t/n) G(t) dt
            # while impedance, according to Chao and Ng, is given by:
            # Z == \int exp(i*2pi*f*t/n) G(t) dt
            Zpl = Zpl.conj()

            if pl == 'll':
                self.Zll = Zpl
            else:
                # Transverse impedance, according to Chao and Ng, is given by:
                # Z == i\int exp(i*2pi*f*t/n) G(t) dt
                setattr(self, 'Z'+pl, 1j*Zpl)
            _log.info('Impedance Calculated.')

        _log.info('#'*60 + '\n')

    def calc_impedance_naff(
            self, pl='ll', s_min=None, s_max=None, win=1, nr_ff=20):
        """."""
        if pl not in self.PLANES:
            raise Exception(
                'Value of variable pl not accepted. Must be one of these: ' +
                ', '.join(self.PLANES))

        # Extracts Needed Variables
        sigt = self.bunlen / c  # bunch time-length
        spos = self.s.copy()
        W = getattr(self, 'W' + pl).copy()
        sizeW = len(W)
        if W is None or not len(W) or _np.all(W == 0):
            raise Exception('No Data found.')

        if s_min is None:
            s_min = spos[0]
        if s_max is None:
            s_max = spos[-1]

        inds = _np.logical_and(spos >= s_min, spos <= s_max)
        spos = spos[inds]
        W = W[inds]
        dt = (spos[1]-spos[0])/c
        leng = len(W) - (len(W)-1) % 6
        spos = spos[-leng:]
        W = W[-leng:]
        if 0.49 < win < 0.51:
            W *= _np.hanning(2*spos.shape[0])[spos.shape[0]:]
            tu, a = _naff.naff_general(
                W, use_win=0, is_real=False, nr_ff=nr_ff)
        elif isinstance(win, int):
            tu, a = _naff.naff_general(
                W, use_win=win, is_real=False, nr_ff=nr_ff)
        else:
            raise Exception(
                'Win must be 1/2 for half-hanning window or an integer '
                'for other windows(0 --> no window).')

        freq = tu/dt
        w = 2*_np.pi*freq
        # Longitudinal position shift to match center of the bunch with zero z:
        a *= _np.exp(-1j*w*(spos[0])/c)

        # Reconstruct the signal
        S = _np.zeros(leng, dtype=complex)
        for wi, ai in zip(w, a):
            S += ai*_np.exp(1j*wi*spos/c)
        S = S.real

        # Deconvolve the Transform with a gaussian bunch:
        a /= _np.exp(-(w*sigt)**2/2)

        # Must multiply by the vector length due to difference in the meaning
        # of the amplitune in the NAFF transform and the fourier transform
        Z = a*dt*sizeW

        # I have to take the conjugate of the fft because:
        # fftt == \int exp(-i*2pi*f*t/n) G(t) dt
        # while impedance, according to Chao and Ng, is given by:
        # Z == \int exp(i*2pi*f*t/n) G(t) dt
        Z = Z.conj()

        if pl != 'll':
            # the Transverse impedance, according to Chao and Ng, is given by:
            # Z == i\int exp(i*2pi*f*t/n) G(t) dt
            Z = 1j*Z
        return freq, Z, leng, S

    def save_processed_data(self, silent=False, pth2sv=None):
        """."""
        _log.basicConfig(level=_log.CRITICAL if silent else _log.INFO)
        _log.info('#'*60 + '\nSaving Processed data:')
        spos = self.s
        freq = self.freq

        if pth2sv is None:
            _log.info('Saving in the same folder of the raw data')
            pth2sv = _abspth('.')
        elif type(pth2sv) is str:
            _log.info('Saving to subfolder: ' + pth2sv)
            if not _isdir(pth2sv):
                _log.info('Folder does not exist. Creating it...')
                _mkdir(pth2sv)
        else:
            msg = 'pth2sv must be a string or None object'
            _log.info(msg)
            raise Exception(msg)

        # Save wakes
        for par in self.PLANES:
            unit = 'V/C' if par == 'll' else 'V/C/m'
            header = '{0:30s} {1:30s}'.format(
                's [m]', 'W{0:s} [{1:s}]'.format(par, unit))
            fname = _jnpth(pth2sv, 'W'+par+'.gz')
            wake = getattr(self, 'W'+par)
            if wake is None or _np.all(wake == 0):
                continue
            _log.info('Saving W' + par + ' data to .gz file')
            _np.savetxt(
                fname, _np.array([spos, wake]).transpose(),
                fmt=['%30.16g', '%30.16g'], header=header)

        # Save Impedances
        for par in self.PLANES:
            unit = 'Ohm' if par == 'll' else 'Ohm/m'
            header = '{0:30s} {1:30s} {2:30s}'.format(
                'Frequency [GHz]',
                'ReZ{0:s} [{1:s}]'.format(par, unit),
                'ImZ{0:s} [{1:s}]'.format(par, unit))
            fname = _jnpth(pth2sv, 'Z'+par+'.gz')
            Z = getattr(self, 'Z'+par)
            if Z is None or _np.allclose(Z, 0, atol=0):
                continue
            _log.info('Saving Z' + par + ' data to .gz file')
            _np.savetxt(
                fname, _np.array([freq/1e9, Z.real, Z.imag]).transpose(),
                fmt=['%30.16g', '%30.16g', '%30.16g'], header=header)

        _log.info(
            'Saving the Complete EMSimulData structure to a .pickle file.')
        with _gzip.open(_jnpth((pth2sv, self.DEFAULT_FNAME_SAVE)), 'wb') as f:
            _pickle.dump(self, f, _pickle.HIGHEST_PROTOCOL)

        _log.info('All Data Saved\n' + '#'*60)

    @staticmethod
    def load_processed_data(filename):
        """."""
        with _gzip.open(filename, 'rb') as fh:
            simul_data = _pickle.load(fh)
        return simul_data

    def plot_wakes(self, save_figs=False, pth2sv=None, show=False, pls=None):
        """."""
        sbun = self.sbun
        sigs = self.bunlen
        spos = self.s
        if pls is None:
            pls = self.PLANES
        for pl in pls:
            wake = getattr(self, 'W'+pl)*1e-12  # V/C -> V/pC
            if wake is None or _np.allclose(wake, 0, atol=0):
                continue
            max_wake = wake[_np.abs(wake).argmax()]
            bunchshape = self.bun * (max_wake/self.bun.max())

            f, axs = _plt.subplots(
                nrows=1, ncols=2, sharey=True, figsize=(14, 6))
            ax = axs[0]
            b = ax.get_position()
            b.x0, b.x1 = 0.05, 0.45
            ax.set_position(b)
            ax.plot(
                sbun*1000, bunchshape, 'b', linewidth=2, label='Bunch Shape')
            ax.plot(spos*1000, wake, 'r', linewidth=2, label='Wake Potential')
            ax.grid(True)
            ax.set_ylabel(self.WAKE_YLABELS[pl], fontsize=13)
            ax.set_xlim([spos[0]*1000, 8000*sigs])
            ax.set_ylim([wake.min()*1.1, wake.max()*1.1])
            ax.legend(loc='best')

            ax = axs[1]
            b = ax.get_position()
            b.x0, b.x1 = 0.45, 0.95
            ax.set_position(b)
            ax.plot(spos*1000, wake, 'r', linewidth=2)
            ax.grid(True)
            tit = ax.set_title(self.TITLES[pl], fontsize=13)
            tit.set_x(0.1)
            xl = ax.set_xlabel('s [mm]', fontsize=13)
            xl.set_x(0.1)
            ax.set_xlim([8000*sigs, spos[-1]*1000])
            ax.set_ylim([wake.min()*1.1, wake.max()*1.1])

            if save_figs:
                f.savefig(_jnpth((pth2sv, 'W'+pl+'.svg')))
        if show:
            _plt.show()

    def plot_impedances(
            self, save_figs=False, pth2sv=None, show=False, pls=None):
        """."""
        freq = self.freq
        if pls is None:
            pls = self.PLANES
        for pl in pls:
            Z = getattr(self, 'Z'+pl)
            if Z is None or _np.allclose(Z, 0, atol=0):
                continue

            _plt.figure()
            _plt.plot(freq/1e9, Z.real, 'r', linewidth=2, label='Re')
            _plt.plot(freq/1e9, Z.imag, 'b--', linewidth=2, label='Im')
            _plt.xlabel('Frequency [GHz]', fontsize=13)
            _plt.grid(True)
            _plt.title(self.TITLES[pl], fontsize=13)
            _plt.ylabel(self.IMPS_YLABELS[pl], fontsize=13)
            _plt.legend(loc='best')
            _plt.xlim(freq[[0, -1]]/1e9)
            if save_figs:
                _plt.savefig(_jnpth((pth2sv, 'Z'+pl+'.svg')))
        if show:
            _plt.show()

    def plot_losskick_factors(
            self, save_figs=False, pth2sv=None, show=False, pls=None):
        """."""
        # Extracts and Initialize Needed Variables:
        self._si.nom_cur = 500e-3
        bunlen = self.bunlen
        sigi = _np.linspace(bunlen, 18e-3, num=50)
        fill_pat = _np.array([1, 864, 864/2, 864/4], dtype=int)

        if pls is None:
            pls = self.PLANES
        pls2 = []
        for pl in pls:
            W = getattr(self, 'W'+pl)
            if W is None or _np.allclose(W, 0, atol=0):
                continue
            Z = getattr(self, 'Z'+pl)
            if Z is None or _np.allclose(Z, 0, atol=0):
                continue
            pls2.append(pl)

        for pl in pls2:
            if pl == 'll':
                fname = 'Loss_factor'
                fig, _ = self._get_figure_loss_factor(sigi, fill_pat, bunlen)
            else:
                fname = 'Kck'+pl+'_factor'
                fig, _ = self._get_figure_kick_factor(
                    sigi, fill_pat, bunlen, pl)
            if save_figs:
                fig.savefig(_jnpth((pth2sv, fname+'.svg')))
        if show:
            _plt.show()

    @classmethod
    def create_make_fig_file(cls, path=None):
        """."""
        if path is None:
            path = _abspth('.')
        fname = _jnpth(path, 'create_figs.py')
        analysis = '#!/usr/bin/env python3\n\n'
        analysis += 'from os import path\n'
        analysis += 'import matplotlib.pyplot as mplt\n'
        analysis += 'import pycolleff.process_wakes as ems\n\n'
        analysis += 'opts = dict(save_figs=False, show=False)\n'
        analysis += 'pth = path.abspath(__file__).rpartition(path.sep)[0]\n'
        analysis += "file_name = path.sep.join([path,'{0:s}'])\n".format(
            cls.DEFAULT_FNAME_SAVE)
        analysis += 'simul_data = ems.load_processed_data(file_name)\n'
        analysis += 'ems.plot_wakes(simul_data, **opts)\n'
        analysis += 'ems.plot_impedances(simul_data, **opts)\n'
        analysis += 'ems.plot_losskick_factors(simul_data, **opts)\n'
        analysis += 'mplt.show()\n'
        with open(fname, 'w') as f:
            f.writelines(analysis)

        permissions = _stat.S_IMODE(_lstat(fname).st_mode)
        permissions |= _stat.S_IXUSR | _stat.S_IXGRP | _stat.S_IXOTH
        _chmod(path, permissions)

    # ---------------------- Auxiliary Methods ----------------------

    @staticmethod
    def _choose_fft_length(n, max_factor=1000):
        for p in range(n):
            n2, i = n-p, 2
            while (i * i <= n2 and i < max_factor):
                if n2 % i:
                    i += 1
                else:
                    n2 //= i
                if n2 < max_factor:
                    return p, n-p

    @staticmethod
    def _get_impedance(spos, wake, sigt, cutoff):
        # frequency scale (Hz):
        dt = (spos[-1]-spos[0]) / (spos.shape[0]-1) / c
        # fft == \int exp(-i*2pi*f*t/n) G(t) dt:
        VHat = _np.fft.fft(wake, wake.shape[0]) * dt
        freq = _np.fft.fftfreq(wake.shape[0], d=dt)
        VHat = _np.fft.fftshift(VHat)  # shift the negative frequencies
        freq = _np.fft.fftshift(freq)  # to the center of the spectrum
        # Longitudinal position shift to match center of the bunch with zero z:
        w = 2*_np.pi*freq
        VHat *= _np.exp(-1j*w*(spos[0])/c)
        # find the maximum useable frequency
        wmax = cutoff/sigt
        indcs = _np.abs(w) <= wmax
        # Deconvolve the Transform with a gaussian bunch:
        Jwlist = _np.exp(-(w*sigt)**2/2)
        Z = VHat[indcs]/Jwlist[indcs]
        return freq[indcs], Z

    def _get_figure_kick_factor(self, sigi, fill_pat, bunlen, pl):
        fig = _plt.figure(figsize=(6, 6))
        ax = _plt.axes()
        for i in range(fill_pat.shape[0]):
            kZi = _np.zeros(sigi.shape[0])
            for j in range(sigi.shape[0]):
                kZi[j] = self.get_kick_factorZ(
                    pl=pl, bunlen=sigi[j],
                    nbunches=fill_pat[i]) * 1e-12  # V/pC/m
            ax.plot(
                sigi*1e3, kZi, 'o', markersize=4,
                label=r'n = {0:03d}'.format(fill_pat[i]))

        # Calculates kickW:
        kickW = self.get_kick_factorW(pl=pl) * 1e-12
        # Print loss factor calculated in both ways
        ax.plot(
            bunlen*1e3, kickW, '*', markersize=7, color=[1, 0, 0],
            label=r'$K_{{{0:s}_{1:s}}}^W$'.format(
                pl[0].upper(), pl[1]))
        ax.set_title('Kick Factor for $n$ equally spaced bunches.')
        ax.set_xlabel(r'$\sigma$ [mm]', fontsize=13)
        ax.set_ylabel(
            r'$K_{{{0:s}_{1:s}}}$ [V/pC/m]'.format(
                pl[0].upper(), pl[1]),
            fontsize=13)
        ax.legend(loc='best')
        ax.grid(True)
        stri = r'$K_{{{0:s}_{1:s}}}^W = {2:5.2f}$ V/pC/m'.format(
            pl[0].upper(), pl[1], kickW)
        ax.annotate(stri, xy=(bunlen*1.1e3, kickW), fontsize=13)
        return fig, [ax, ]

    def _get_figure_loss_factor(self, sigi, fill_pat, bunlen):
        # bunch length scenarios:
        sigvec = _np.array([2.65, 5, 8, 10, 15], dtype=float)*1e-3
        # current scenarios
        Ivec = _np.linspace(10e-3, self._si.nom_cur, num=50)

        fig, axs = _plt.subplots(
            nrows=1, ncols=2, figsize=(12, 6),
            gridspec_kw=dict(left=0.08, right=0.97))
        ax = axs[0]
        for i in range(fill_pat.shape[0]):
            kZi = _np.zeros(sigi.shape[0])
            for j in range(sigi.shape[0]):
                # V/pC:
                kZi[j] = self.get_klossZ(
                    bunlen=sigi[j], nbunches=fill_pat[i]) * 1e-12
            ax.semilogy(
                sigi * 1e3, kZi * 1e3, 'o', markersize=4,
                label=r'$n = {0:03d}$'.format(fill_pat[i]))
        # Calculates klossW
        kW = self.get_klossW() * 1e-12
        # Print loss factor calculated in both ways
        ax.semilogy(
            bunlen * 1e3, kW * 1e3, '*', markersize=7, color=[1, 0, 0],
            label=r'$K_L^W$')
        ax.set_title('Loss Factor for $n$ equally spaced bunches.')
        ax.set_xlabel(r'$\sigma$ [mm]')
        ax.set_ylabel(r'$K_L$ [mV/pC]')
        ax.legend(loc='best')
        ax.grid(True)
        ax.annotate(
            r'$K_L^W = {0:5.2f}$ mV/pC'.format(kW*1e3),
            xy=(bunlen*1.1e3, kW*1e3), fontsize=12)
        ax = axs[1]
        kZvec = _np.zeros(sigvec.shape[0])
        labels = []
        for i in range(sigvec.shape[0]):
            kZvec[i] = self.get_klossZ(bunlen=sigvec[i])  # V/C
            labels.append(r'$\sigma = {0:05.2f}$ mm'.format(
                sigvec[i]*1e3))
        Plossvec = kZvec[None, :] * Ivec[:, None]**2
        Plossvec *= self._si.T0 / self._si.harm_num

        ax.semilogy(Ivec*1e3, Plossvec, markersize=4)
        ax.set_title(
            'Power Loss for ${0:d}$ equally spaced bunches.'.format(
                self._si.harm_num))
        ax.set_xlabel(r'$I_{{avg}}$ [mA]')
        ax.set_ylabel(r'Power [W]')
        ax.legend(labels, loc='best')
        ax.grid(True)
        return fig, axs
